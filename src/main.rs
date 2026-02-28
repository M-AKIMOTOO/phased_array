mod acf;
mod args;
mod cor;
mod geom;
mod ifile;
mod plot;
mod utils;
mod xcf;
mod xml;

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::{mpsc, Arc};
use std::thread;

use clap::Parser;
use num_complex::Complex;
use rayon::prelude::*;

use acf::finalize_auto_spectrum;
use args::{parse_levels, parse_shuffle, resolve_per_antenna_config, resolve_weight, DEFAULT_SHUFFLE_IN};
use cor::{epoch_to_yyyydddhhmmss, CorHeaderConfig, CorStation, CorWriter};
use plot::{plot_multi_series_f64_x, plot_series_f64_x, plot_series_with_x, BLUE, GREEN, RED};
use utils::{apply_delay_and_rate_regular_bins, apply_integer_sample_shift_zerofill, build_decode_plan, decode_block_into_with_plan, quantise_frame, DynError, FftHelper};
use xcf::finalize_cross_spectrum;

const DEFAULT_COARSE_DELAY_S: f64 = 1.6e-6;

fn carrier_phase_from_delay(f_hz: f64, tau_s: f64) -> Complex<f64> {
    if f_hz == 0.0 {
        Complex::new(1.0, 0.0)
    } else {
        Complex::from_polar(1.0, -2.0 * std::f64::consts::PI * f_hz * tau_s)
    }
}

fn delay_seconds_at_time(d_s: f64, r_sps: f64, a_sps2: f64, t_s: f64) -> f64 {
    d_s + r_sps * t_s + 0.5 * a_sps2 * t_s.powi(2)
}

fn split_delay_to_integer_and_fractional(delay_seconds: f64, fs_hz: f64) -> (i64, f64) {
    let integer_samples = (delay_seconds * fs_hz).round() as i64;
    let fractional_seconds = delay_seconds - (integer_samples as f64 / fs_hz);
    (integer_samples, fractional_seconds)
}

#[derive(Clone, Copy, Debug)]
enum OutputGrid { Ant1, Ant2 }

const VSREC_LEVELS_2BIT: &str = "-1.5,-0.5,0.5,1.5";
const VSREC_SHUFFLE_EXTERNAL: &str =
    "24,25,26,27,28,29,30,31,16,17,18,19,20,21,22,23,8,9,10,11,12,13,14,15,0,1,2,3,4,5,6,7";

#[derive(Clone, Copy, Debug)]
struct BandAlignment { shift_bins: isize, a1s: usize, a1e: usize, a2s: usize }

fn compute_band_alignment(fft: usize, fs: f64, c1: f64, c2: f64, bw1: f64, bw2: f64) -> Result<BandAlignment, DynError> {
    let df = fs / fft as f64; let h = fft / 2 + 1;
    let lim = if fft % 2 == 0 { h - 1 } else { h };
    let v1 = (((bw1 * 1e6) / df).floor() as usize).saturating_add(1).min(lim);
    let v2 = (((bw2 * 1e6) / df).floor() as usize).saturating_add(1).min(lim);
    let sr = (c1 * 1e6 - c2 * 1e6) / df; let sb = sr.round() as isize;
    let a1s = 0isize.max(-sb) as usize; let a1e = (v1 as isize).min(v2 as isize - sb) as usize;
    if a1e <= a1s { return Err("No band overlap".into()); }
    Ok(BandAlignment { shift_bins: sb, a1s, a1e, a2s: (a1s as isize + sb) as usize })
}

fn format_bit_codes(b: usize) -> String { if b == 0 { "n/a".into() } else { (0..(1<<b).min(1024)).map(|c| format!("{c:0b$b}", b=b)).collect::<Vec<_>>().join(" ") } }
fn format_level_map(b: usize, l: &[f64]) -> String { if b == 0 { "n/a".into() } else { (0..(1<<b).min(l.len())).map(|c| format!("{c:0b$b}->{:.6}", l[c], b=b)).collect::<Vec<_>>().join(", ") } }

fn normalize_level_args(level_args: &[String], bit1: usize, bit2: usize) -> Result<Vec<String>, DynError> {
    if level_args.is_empty() {
        return Ok(level_args.to_vec());
    }
    if level_args
        .iter()
        .any(|v| v.contains("ant1:") || v.contains("ant2:"))
    {
        return Ok(level_args.to_vec());
    }
    if bit1 != bit2 {
        return Ok(level_args.to_vec());
    }
    let expected = 1usize << bit1;
    let all_numeric = level_args.iter().all(|v| v.parse::<f64>().is_ok());
    if all_numeric {
        if level_args.len() != expected {
            return Err(format!(
                "--level without ant1:/ant2: must provide exactly {expected} values for --bit {bit1}"
            )
            .into());
        }
        return Ok(vec![level_args.join(",")]);
    }
    Ok(level_args.to_vec())
}

fn normalize_shuffle_args(shuffle_args: &[String]) -> Result<Vec<String>, DynError> {
    if shuffle_args.is_empty() {
        return Ok(shuffle_args.to_vec());
    }
    if shuffle_args
        .iter()
        .any(|v| v.contains("ant1:") || v.contains("ant2:"))
    {
        return Ok(shuffle_args.to_vec());
    }
    let all_numeric = shuffle_args.iter().all(|v| v.parse::<usize>().is_ok());
    if all_numeric {
        if shuffle_args.len() != 32 {
            return Err("--shuffle without ant1:/ant2: must provide exactly 32 values".into());
        }
        return Ok(vec![shuffle_args.join(",")]);
    }
    Ok(shuffle_args.to_vec())
}

fn seek_forward_samples(r: &mut BufReader<File>, s: u64, b: usize, len: u64, sought: &mut isize, label: &str, reason: &str) -> Result<(), DynError> {
    let bytes = (s * b as u64 + 7) / 8;
    if bytes >= len { return Err(format!("Seek for {} exceeds file", label).into()); }
    r.seek(SeekFrom::Current(bytes as i64))?;
    *sought += ((bytes * 8) / b as u64) as isize;
    println!("[info] Seeking {} forward by {} bytes for {}", label, bytes, reason);
    Ok(())
}

fn read_with_padding(reader: &mut BufReader<File>, buf: &mut [u8]) -> Result<usize, DynError> {
    let mut total = 0;
    while total < buf.len() {
        match reader.read(&mut buf[total..]) {
            Ok(0) => { buf[total..].fill(0); return Ok(total); }
            Ok(n) => total += n,
            Err(ref e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            Err(e) => return Err(e.into()),
        }
    }
    Ok(total)
}

fn apply_integer_delay_with_forward_seek(ds: i64, b1: usize, b2: usize, res: &str, r1: &mut BufReader<File>, r2: &mut BufReader<File>, l1: u64, l2: u64, s1: &mut isize, s2: &mut isize) -> Result<(), DynError> {
    if ds == 0 { return Ok(()); }
    if ds > 0 { seek_forward_samples(r2, ds as u64, b2, l2, s2, "ant2", res)?; }
    else { seek_forward_samples(r1, (-ds) as u64, b1, l1, s1, "ant1", res)?; }
    Ok(())
}

fn resolve_input_paths(args: &args::Args, pe: &str, meta: Option<&ifile::IFileData>) -> Result<(PathBuf, PathBuf, String, String), DynError> {
    if let (Some(a1), Some(a2)) = (&args.ant1, &args.ant2) { let (_, tag) = epoch_to_yyyydddhhmmss(pe)?; return Ok((a1.clone(), a2.clone(), pe.to_string(), tag)); }
    let data_dir = args.data.clone().unwrap_or_else(|| PathBuf::from("."));
    let (_, tag) = epoch_to_yyyydddhhmmss(pe)?;
    let mut candidates = vec![("YAMAGU32", "YAMAGU34")];
    if let Some(m) = meta { if let (Some(n1), Some(n2)) = (&m.ant1_station_name, &m.ant2_station_name) { candidates.insert(0, (n1, n2)); } }
    for (p1, p2) in candidates {
        let a1 = data_dir.join(format!("{}_{}.raw", p1, tag)); let a2 = data_dir.join(format!("{}_{}.raw", p2, tag));
        if a1.exists() && a2.exists() { println!("[info] Auto-resolved inputs: {} / {}", a1.display(), a2.display()); return Ok((a1, a2, pe.to_string(), tag)); }
    }
    Err("Input files not found".into())
}

fn resolve_output_layout(args: &args::Args, tag: &str) -> Result<(PathBuf, PathBuf), DynError> {
    let stem = format!("YAMAGU66_{tag}");
    let dir = args.output.as_ref().and_then(|p| p.parent()).map(|p| p.to_path_buf()).unwrap_or_else(|| PathBuf::from("phased_array").join(&stem));
    std::fs::create_dir_all(&dir)?;
    let path = args.output.clone().unwrap_or_else(|| dir.join(format!("{stem}.raw")));
    Ok((dir, path))
}

fn main() -> Result<(), DynError> {
    let args = args::Args::parse();
    if args.cpu == 0 {
        return Err("--cpu must be >= 1".into());
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.cpu)
        .build_global()
        .map_err(|e| format!("failed to configure rayon thread pool: {e}"))?;
    let if_d = if let Some(p) = &args.ifile { Some(ifile::parse_ifile(p)?) } else { None };
    let fft_len = if args.fft == args::DEFAULT_FFT { if_d.as_ref().and_then(|d| d.fft).unwrap_or(args.fft) } else { args.fft };
    
    let bit_a = if args.bit.is_empty() { vec!["2".into()] } else { args.bit.clone() };
    let (bit1, bit2) = resolve_per_antenna_config(&bit_a, if_d.as_ref().and_then(|d| d.ant1_bit).unwrap_or(2), |s| Ok(s.parse()?))?;
    if args.vsrec && (bit1 != 2 || bit2 != 2) {
        return Err("--vsrec requires 2-bit input on both antennas".into());
    }
    let bit_out = bit1.max(bit2);
    let level_src = if args.vsrec {
        vec![VSREC_LEVELS_2BIT.to_string()]
    } else {
        args.level.clone()
    };
    let level_args = normalize_level_args(&level_src, bit1, bit2)?;
    let (lv1_s, lv2_s) = resolve_per_antenna_config(&level_args, if_d.as_ref().and_then(|d| d.ant1_level.clone()).unwrap_or("-1.5,-0.5,0.5,1.5".into()), |s| Ok(s.to_string()))?;
    let (levels1, levels2) = (Arc::new(parse_levels(bit1, &lv1_s)?), Arc::new(parse_levels(bit2, &lv2_s)?));
    let ds_s = DEFAULT_SHUFFLE_IN.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(",");
    let shuffle_src = if args.vsrec {
        vec![VSREC_SHUFFLE_EXTERNAL.to_string()]
    } else {
        args.shuffle_in.clone()
    };
    let shuffle_args = normalize_shuffle_args(&shuffle_src)?;
    let (sh1_s, sh2_s) = resolve_per_antenna_config(&shuffle_args, if_d.as_ref().and_then(|d| d.ant1_shuffle.clone()).unwrap_or(ds_s), |s| Ok(s.to_string()))?;
    let (sh1, sh2) = (Arc::new(parse_shuffle(&sh1_s)?), Arc::new(parse_shuffle(&sh2_s)?));
    let sh1_ext: Vec<usize> = sh1_s.split(',').map(|v| v.trim().parse::<usize>().unwrap()).collect();
    let sh2_ext: Vec<usize> = sh2_s.split(',').map(|v| v.trim().parse::<usize>().unwrap()).collect();
    
    let sb_a = if args.sideband.is_empty() { vec!["LSB".into()] } else { args.sideband.clone() };
    let (sb1_s, sb2_s) = resolve_per_antenna_config(&sb_a, if_d.as_ref().and_then(|d| d.ant1_sideband.clone()).unwrap_or("LSB".into()), |s| Ok(s.to_uppercase()))?;
    let (lsb1, lsb2) = (sb1_s == "LSB", sb2_s == "LSB");
    let output_lsb = lsb1 && lsb2;

    let (tsys1, tsys2) = resolve_per_antenna_config(&args.tsys, 1.0, |s| Ok(s.parse()?))?;
    let (dia1, dia2) = resolve_per_antenna_config(&args.diameter, 0.0, |s| Ok(s.parse()?))?;
    let (eta1, eta2) = resolve_per_antenna_config(&args.eta, 0.65, |s| Ok(s.parse()?))?;
    let (gain1, gain2) = resolve_per_antenna_config(&args.gain, 1.0, |s| Ok(s.parse()?))?;
    let (sefd1, sefd2) = resolve_per_antenna_config(&args.sefd, 0.0, |s| Ok(s.parse()?))?;

    let fs = if_d.as_ref().and_then(|d| d.sampling_hz).unwrap_or(args.sampling * 1e6);
    let obs_mhz = if_d.as_ref().and_then(|d| d.obsfreq_mhz).unwrap_or(args.obsfreq);
    let ep_i = args.epoch.clone().or_else(|| if_d.as_ref().and_then(|d| d.epoch.clone())).unwrap_or("2000".into());
    let (a1p, a2p, _, c_tag) = resolve_input_paths(&args, &ep_i, if_d.as_ref())?;
    let (c_unix, _) = epoch_to_yyyydddhhmmss(&ep_i)?;
    let (o_dir, o_path) = resolve_output_layout(&args, &c_tag)?;
    let ant1_ecef = if_d.as_ref().and_then(|d| d.ant1_ecef_m).unwrap_or(geom::YAMAGU32_ECEF);
    let ant2_ecef = if_d.as_ref().and_then(|d| d.ant2_ecef_m).unwrap_or(geom::YAMAGU34_ECEF);

    let mut gdi: Option<GDI> = None; struct GDI { ra: f64, dec: f64, mjd: f64 }
    let (mut gd0, mut gr0, mut ga0) = (0.0, 0.0, 0.0);
    if let (Some(ra_s), Some(dec_s)) = (args.ra.clone().or_else(|| if_d.as_ref().map(|d| d.ra.clone())), args.dec.clone().or_else(|| if_d.as_ref().map(|d| d.dec.clone()))) {
        let (ra_raw, dec_raw, mjd) = (
            geom::parse_ra(&ra_s)?,
            geom::parse_dec(&dec_s)?,
            geom::parse_epoch_to_mjd(&ep_i)?,
        );
        // Input sky coordinates are always interpreted as J2000.
        let (ra, dec) = geom::precess_j2000_to_mean_of_date(ra_raw, dec_raw, mjd);
        let (_, _, gd, gr, ga) = geom::calculate_geometric_delay_and_derivatives(ant1_ecef, ant2_ecef, ra, dec, mjd);
        gdi = Some(GDI { ra, dec, mjd }); (gd0, gr0, ga0) = (gd, gr, ga);
    }
    let clock_delay_rel_legacy_s = if_d.as_ref().and_then(|d| d.clock_delay_s).unwrap_or(0.0);
    let clock_rate_rel_legacy_sps = if_d.as_ref().and_then(|d| d.clock_rate_sps).unwrap_or(0.0);
    let clock1_delay_s = if_d.as_ref().and_then(|d| d.ant1_clock_delay_s).unwrap_or(0.0);
    let clock2_delay_s = if_d
        .as_ref()
        .and_then(|d| d.ant2_clock_delay_s)
        .unwrap_or(clock1_delay_s + clock_delay_rel_legacy_s);
    let clock1_rate_sps = if_d.as_ref().and_then(|d| d.ant1_clock_rate_sps).unwrap_or(0.0);
    let clock2_rate_sps = if_d
        .as_ref()
        .and_then(|d| d.ant2_clock_rate_sps)
        .unwrap_or(clock1_rate_sps + clock_rate_rel_legacy_sps);
    let clock_delay_s = clock2_delay_s - clock1_delay_s;
    let clock_rate_sps = clock2_rate_sps - clock1_rate_sps;
    let coarse_delay_s = args.coarse.unwrap_or(DEFAULT_COARSE_DELAY_S);
    let net_d_rel_no_clock0 = gd0 + args.gico3_correct + coarse_delay_s + args.delay/fs;
    let net_d0 = net_d_rel_no_clock0 + clock_delay_s;
    
    let rot1 = if_d.as_ref().and_then(|d| d.ant1_rotation_hz).unwrap_or(0.0);
    let rot2 = if_d.as_ref().and_then(|d| d.ant2_rotation_hz).unwrap_or(0.0);
    let bw = fs / 2e6;
    let obs_hz = obs_mhz * 1e6;
    let bw_hz = fs / 2.0;
    let a1_data_low = obs_mhz + rot1/1e6; let a2_data_low = obs_mhz + rot2/1e6;
    let lo1_hz = a1_data_low * 1e6;
    let lo2_hz = a2_data_low * 1e6;
    let ba = compute_band_alignment(fft_len, fs, a1_data_low + 0.5*bw, a2_data_low + 0.5*bw, bw, bw)?;
    let out_grid = if rot1.abs() <= rot2.abs() { OutputGrid::Ant1 } else { OutputGrid::Ant2 };
    
    let bpf1 = (fft_len * bit1 + 7) / 8; let bpf2 = (fft_len * bit2 + 7) / 8;
    let bpf_o = (fft_len * bit_out + 7) / 8;
    let f1_m = std::fs::metadata(&a1p)?; let f2_m = std::fs::metadata(&a2p)?;
    
    let total_sec = args.length.unwrap_or_else(|| {
        let s1 = f1_m.len() as f64 * 8.0 / bit1 as f64 / fs;
        let s2 = f2_m.len() as f64 * 8.0 / bit2 as f64 / fs;
        s1.min(s2)
    });
    // Do not exceed the requested/available duration.
    // A full FFT frame is required, so we keep only complete frames.
    let total_f = (total_sec * fs / fft_len as f64).floor() as usize;

    let (w1, _, _, _) = resolve_weight(tsys1, gain1, if sefd1 > 0.0 { Some(sefd1) } else { None }, if dia1 > 0.0 { Some(dia1) } else { None }, eta1, "A1")?;
    let (w2, _, _, _) = resolve_weight(tsys2, gain2, if sefd2 > 0.0 { Some(sefd2) } else { None }, if dia2 > 0.0 { Some(dia2) } else { None }, eta2, "A2")?;
    let a1_name = if_d.as_ref().and_then(|d| d.ant1_station_name.as_deref()).unwrap_or("YAMAGU32");
    let a2_name = if_d.as_ref().and_then(|d| d.ant2_station_name.as_deref()).unwrap_or("YAMAGU34");
    let a1_key_opt = if_d.as_ref().and_then(|d| d.ant1_station_key.as_deref());
    let a2_key_opt = if_d.as_ref().and_then(|d| d.ant2_station_key.as_deref());
    let a1_key = a1_key_opt.unwrap_or("-");
    let a2_key = a2_key_opt.unwrap_or("-");
    let a1_code = a1_key_opt.and_then(|k| k.as_bytes().first().copied()).or_else(|| a1_name.as_bytes().first().copied()).unwrap_or(b'1');
    let a2_code = a2_key_opt.and_then(|k| k.as_bytes().first().copied()).or_else(|| a2_name.as_bytes().first().copied()).unwrap_or(b'2');

    println!("Starting phased array processing with the following arguments:");
    println!("--------------------------------------------------");
    println!("  {}:  {} (Size: {} bytes, Estimated Obs Time: {:.2}s)", a1_name, a1p.display(), f1_m.len(), (f1_m.len()*8) as f64 / bit1 as f64 / fs);
    println!("  {}:  {} (Size: {} bytes, Estimated Obs Time: {:.2}s)", a2_name, a2p.display(), f2_m.len(), (f2_m.len()*8) as f64 / bit2 as f64 / fs);
    if let Some(info) = &gdi {
        println!("  ra/dec:     {:.6} / {:.6}", info.ra.to_degrees(), info.dec.to_degrees());
        println!("  source-frame: J2000 (fixed)");
        println!("  epoch:      {}", ep_i);
        println!("  geom-delay: {:.6e} s ({} - {})", gd0, a2_name, a1_name);
        println!("  geom-rate:  {:.6e} s/s ({} - {}) => {:.6} Hz @ obsfreq", gr0, a2_name, a1_name, gr0 * obs_mhz * 1e6);
        println!("  geom-accel: {:.6e} s/s^2 ({} - {})", ga0, a2_name, a1_name);
    }
    println!("  coarse-delay fixed: {:.6e} s (relative pre-align) => applied {:.3} samples", coarse_delay_s, coarse_delay_s*fs);
    println!("  clock-delay {}: {:.6e} s => applied {:.3} samples", a1_name, clock1_delay_s, clock1_delay_s * fs);
    println!("  clock-delay {}: {:.6e} s => applied {:.3} samples", a2_name, clock2_delay_s, clock2_delay_s * fs);
    println!("  clock-rate  {}: {:.6e} s/s", a1_name, clock1_rate_sps);
    println!("  clock-rate  {}: {:.6e} s/s", a2_name, clock2_rate_sps);
    println!("  res-delay input: {} samples (relative pre-align)", args.delay);
    println!("  read-align delay: {:.3} samples ({:.3e} s)", net_d0 * fs, net_d0);
    println!("  delay-model: per-frame midpoint delay + integer/fractional correction");
    let correction_sign = -1.0f64; // Correct ant1 toward ant2 while geometric model is (ant2 - ant1).
    let geometric_rate_hz = correction_sign * gr0 * obs_mhz * 1e6;
    let clock_rate_hz = clock_rate_sps * obs_mhz * 1e6;
    let df_hz = fs / fft_len as f64;
    // rotation is handled as fringe-stop:
    //   integer-bin component -> band-align bin shift
    //   residual (sub-bin) component -> phase-rate correction term
    let rotation_delta_hz = rot1 - rot2; // ant1 - ant2
    let rotation_shift_hz = ba.shift_bins as f64 * df_hz;
    let rotation_residual_hz = rotation_delta_hz - rotation_shift_hz;
    let rotation_fringe_hz = -rotation_residual_hz; // correction applied on ant1 toward ant2
    let total_rate_hz = geometric_rate_hz + clock_rate_hz + rotation_fringe_hz + args.rate;
    println!(
        "  delay-rate: {:.6} Hz (geom {:.6} + clock {:.6} + rot-res {:.6} + user {:.6})",
        total_rate_hz, geometric_rate_hz, clock_rate_hz, rotation_fringe_hz, args.rate
    );
    println!("  obsfreq:    {:.6} MHz (Reference)", obs_mhz);
    println!("  rotation:   {:.3} Hz ({}), {:.3} Hz ({})", rot1, a1_name, rot2, a2_name);
    println!(
        "  phase-freq: {:.6} MHz ({}) / {:.6} MHz ({}) (delay/rate correction carriers)",
        a1_data_low, a1_name, a2_data_low, a2_name
    );
    println!("  bw:         {:.3} MHz", bw);
    println!(
        "  {}-param: key={} sideband={} obsfreq_hz={:.3} rotation_hz={:.3} ref_band_low_hz={:.3} data_band_low_hz={:.3} data_band_center_hz={:.3} bw_hz={:.3}",
        a1_name, a1_key, sb1_s, obs_hz, rot1, obs_hz, a1_data_low * 1e6, (a1_data_low + 0.5 * bw) * 1e6, bw_hz
    );
    println!(
        "  {}-param: key={} sideband={} obsfreq_hz={:.3} rotation_hz={:.3} ref_band_low_hz={:.3} data_band_low_hz={:.3} data_band_center_hz={:.3} bw_hz={:.3}",
        a2_name, a2_key, sb2_s, obs_hz, rot2, obs_hz, a2_data_low * 1e6, (a2_data_low + 0.5 * bw) * 1e6, bw_hz
    );
    println!("  sampling:   {:.0} Hz ({:.6} MHz)", fs, fs/1e6);
    println!("  samples/s:  {:.0}", fs);
    println!("  frames/s:   {:.6}", fs/fft_len as f64);
    println!("  fft:        {}", fft_len);
    println!("  debug:      {}", args.debug);
    println!("  bit:        {}={} {}={} -> out={}", a1_name, bit1, a2_name, bit2, bit_out);
    println!("  vsrec:      {}", args.vsrec);
    println!("  bit-code:   {}=({}) {}=({})", a1_name, format_bit_codes(bit1), a2_name, format_bit_codes(bit2));
    println!("  level:      {}={:?}", a1_name, levels1);
    println!("  level:      {}={:?}", a2_name, levels2);
    println!("  level-map:  {}={}", a1_name, format_level_map(bit1, &levels1));
    println!("  level-map:  {}={}", a2_name, format_level_map(bit2, &levels2));
    println!("  shuffle-in: {}={:?}", a1_name, sh1_ext);
    println!("  shuffle-in: {}={:?}", a2_name, sh2_ext);
    println!("  sideband:   {}={} {}={}", a1_name, sb1_s, a2_name, sb2_s);
    println!("  sideband-normalize: {}={} {}={} (internal USB domain)", a1_name, if lsb1 { "LSB->USB" } else { "USB" }, a2_name, if lsb2 { "LSB->USB" } else { "USB" });
    println!("  sideband-output: {} (YAMAGU66 raw)", if output_lsb { "LSB" } else { "USB" });
    println!("  obs-band:   {:.6} .. {:.6} MHz", obs_mhz, obs_mhz + bw);
    println!("  band-align: {}->{} shift {} bins, overlap {}[{}..{})", a2_name, a1_name, ba.shift_bins, a1_name, ba.a1s, ba.a1e);
    println!("  band-overlap: {} bins ({:.6} MHz, {:.6} MHz/bin)", ba.a1e - ba.a1s, (ba.a1e - ba.a1s) as f64 * fs/fft_len as f64 / 1e6, fs/fft_len as f64 / 1e6);
    println!("  rotation-shift: target {}, grid {}-ref, shift {} bins", a2_name, a1_name, ba.shift_bins);
    println!("  rotation-fringestop: delta_hz={:.6} shift_hz={:.6} residual_hz={:.6}", rotation_delta_hz, rotation_shift_hz, rotation_residual_hz);
    println!("  output-grid: {}", match out_grid { OutputGrid::Ant1 => a1_name, OutputGrid::Ant2 => a2_name });
    println!("  ant-fixed:  {}={:?}, {}={:?}", a1_name, ant1_ecef, a2_name, ant2_ecef);
    println!("  cpu:        {} (compute threads: {})", args.cpu, rayon::current_num_threads());
    println!("  skip:       0");
    let processed_sec = total_f as f64 * fft_len as f64 / fs;
    println!("  length:     {:.6}s requested, {:.6}s processed", total_sec, processed_sec);
    println!("  tsys:       {} ({}), {} ({})", tsys1, a1_name, tsys2, a2_name);
    println!("  eta:        {} ({}), {} ({})", eta1, a1_name, eta2, a2_name);
    println!("  gain:       {} ({}), {} ({})", gain1, a1_name, gain2, a2_name);
    println!("  weight:     {:.6} ({}), {:.6} ({})", w1, a1_name, w2, a2_name);
    println!("  results:    {}", o_dir.display());
    println!("  output:     {}", o_path.display());
    println!("--------------------------------------------------");

    let (mut r1, mut r2) = (BufReader::new(File::open(&a1p)?), BufReader::new(File::open(&a2p)?));
    let (mut s1_s, mut s2_s) = (0, 0); apply_integer_delay_with_forward_seek((net_d0 * fs).round() as i64, bit1, bit2, "init", &mut r1, &mut r2, f1_m.len(), f2_m.len(), &mut s1_s, &mut s2_s)?;
    let d_seek = (s2_s - s1_s) as f64 / fs; let helper = Arc::new(FftHelper::new(fft_len));
    let frame_dt = fft_len as f64 / fs;
    let net_d1_base = (net_d_rel_no_clock0 - clock1_delay_s) - d_seek;
    let net_d2_base = -clock2_delay_s;
    let rate_rel_no_clock_base = (geometric_rate_hz + rotation_fringe_hz + args.rate) / (obs_mhz * 1e6);
    let total_rate1_base = rate_rel_no_clock_base - clock1_rate_sps;
    let total_rate2_base = -clock2_rate_sps;
    let total_accel_base = correction_sign * ga0;
    let total_accel1_base = total_accel_base;
    let total_accel2_base = 0.0;
    let exact_geom_params = gdi.as_ref().map(|v| (v.ra, v.dec, v.mjd));
    let extra_delay_rate_sps = (rotation_fringe_hz + args.rate) / (obs_mhz * 1e6);
    if exact_geom_params.is_some() {
        println!("[info] Delay model refinement: per-frame geometric delay re-evaluation enabled");
    }
    let source_name = if_d.as_ref().and_then(|d| d.source.clone()).unwrap_or_else(|| "UNKNOWN".into());

    if args.fringe {
        let mut ac = vec![Complex::new(0.0, 0.0); fft_len/2+1];
        let (tx, rx) = mpsc::sync_channel::<Vec<Vec<u8>>>(4);
        let (r1_p, r2_p, s1_c, s2_c) = (a1p.clone(), a2p.clone(), s1_s, s2_s);
        thread::spawn(move || {
            let (mut rd1, mut rd2) = (BufReader::new(File::open(r1_p).unwrap()), BufReader::new(File::open(r2_p).unwrap()));
            rd1.seek(SeekFrom::Start(s1_c as u64 * bit1 as u64 / 8)).unwrap(); rd2.seek(SeekFrom::Start(s2_c as u64 * bit2 as u64 / 8)).unwrap();
            let mut nf = 0; while nf < total_f {
                let n = (total_f - nf).min(4096); let (mut b1, mut b2) = (vec![0u8; n * bpf1], vec![0u8; n * bpf2]);
                read_with_padding(&mut rd1, &mut b1).unwrap(); read_with_padding(&mut rd2, &mut b2).unwrap();
                tx.send(vec![b1, b2]).unwrap(); nf += n;
            }
        });
        let (dp1, dp2) = (Arc::new(build_decode_plan(bit1, sh1.as_ref())?), Arc::new(build_decode_plan(bit2, sh2.as_ref())?));
        let mut processed = 0; let inv_fft2 = 1.0 / (fft_len as f64).powi(2);
        for bufs in rx {
            let (raw1, raw2) = (&bufs[0], &bufs[1]); let nf = raw1.len() / bpf1;
            let bc = raw1.par_chunks(bpf1).zip(raw2.par_chunks(bpf2)).enumerate().fold(|| vec![Complex::new(0.0, 0.0); fft_len/2+1], |mut c, (i, (r1, r2))| {
                let t0 = (processed + i) as f64 * frame_dt;
                let t_mid = t0 + 0.5 * frame_dt;
                let (mut f1, mut f2, mut s1, mut s2) = (vec![0.0; fft_len], vec![0.0; fft_len], vec![Complex::new(0.0, 0.0); fft_len/2+1], vec![Complex::new(0.0, 0.0); fft_len/2+1]);
                decode_block_into_with_plan(r1, &levels1, fft_len, &dp1, &mut Vec::new(), &mut f1, lsb1).unwrap();
                decode_block_into_with_plan(r2, &levels2, fft_len, &dp2, &mut Vec::new(), &mut f2, lsb2).unwrap();
                let (tau1, tau2) = if let Some((ra, dec, mjd0)) = exact_geom_params {
                    let mjd_t = mjd0 + t_mid / 86400.0;
                    let (_, _, gd_t, _, _) = geom::calculate_geometric_delay_and_derivatives(ant1_ecef, ant2_ecef, ra, dec, mjd_t);
                    let net_d_rel_no_clock_t =
                        gd_t + args.gico3_correct + coarse_delay_s + args.delay / fs + extra_delay_rate_sps * t_mid;
                    let clock1_t = clock1_delay_s + clock1_rate_sps * t_mid;
                    let clock2_t = clock2_delay_s + clock2_rate_sps * t_mid;
                    ((net_d_rel_no_clock_t - clock1_t) - d_seek, -clock2_t)
                } else {
                    let net_d1 = net_d1_base;
                    let net_d2 = net_d2_base;
                    let total_rate1 = total_rate1_base;
                    let total_rate2 = total_rate2_base;
                    (
                        delay_seconds_at_time(net_d1, total_rate1, total_accel1_base, t_mid),
                        delay_seconds_at_time(net_d2, total_rate2, total_accel2_base, t_mid),
                    )
                };
                let (int_shift1, frac_delay1) = split_delay_to_integer_and_fractional(tau1, fs);
                let (int_shift2, frac_delay2) = split_delay_to_integer_and_fractional(tau2, fs);
                apply_integer_sample_shift_zerofill(&mut f1, int_shift1);
                apply_integer_sample_shift_zerofill(&mut f2, int_shift2);
                helper.forward_r2c_process(&mut f1, &mut s1).unwrap();
                helper.forward_r2c_process(&mut f2, &mut s2).unwrap();
                
                // Use observing reference frequency in .cor header (not upper band edge)
                // so residual-rate definition matches frinZ delay/rate search.
                let fr_lo1 = carrier_phase_from_delay(lo1_hz, tau1);
                let fr_lo2 = carrier_phase_from_delay(lo2_hz, tau2);
                apply_delay_and_rate_regular_bins(&mut s1, fft_len, fs/fft_len as f64, frac_delay1, 0.0, 0.0, 0.0, false);
                apply_delay_and_rate_regular_bins(&mut s2, fft_len, fs/fft_len as f64, frac_delay2, 0.0, 0.0, 0.0, false);
                
                for k in 0..(ba.a1e - ba.a1s) {
                    let i1 = ba.a1s + k;
                    let i2 = ba.a2s + k;
                    match out_grid {
                        OutputGrid::Ant1 => c[i1] += (s1[i1] * fr_lo1) * (s2[i2] * fr_lo2).conj(),
                        OutputGrid::Ant2 => c[i2] += (s1[i1] * fr_lo1) * (s2[i2] * fr_lo2).conj(),
                    }
                }
                c
            }).reduce(|| vec![Complex::new(0.0, 0.0); fft_len/2+1], |mut c1, c2| { for k in 0..c1.len() { c1[k] += c2[k]; } c1 });
            for k in 0..ac.len() { ac[k] += bc[k] * inv_fft2; }
            processed += nf; print!("\rCorrelating ({}/{})", processed, total_f); std::io::stdout().flush().unwrap();
        }
        println!();
        let xcf_res = finalize_cross_spectrum(
            &mut ac,
            helper.as_ref(),
            fft_len,
            fft_len / 2 + 1,
            &(-(fft_len as i32 / 2)..fft_len as i32 / 2).collect::<Vec<_>>(),
            fs,
            obs_mhz,
            true,
            &o_dir,
            false,
        )?;
        if let Some(delay_s) = xcf_res.delay_seconds_from_phase {
            println!("[info] Delay estimate from phase slope: {:.9e} s", delay_s);
        }
    }

    if total_f > 0 {
        println!("[info] Synthesising output...");
        let mut wr = BufWriter::new(File::create(&o_path)?);
        let frame_sec = fft_len as f64 / fs;
        let total_duration_sec = total_f as f64 * frame_sec;
        // "second sectors" should track actual covered duration without
        // creating an extra tail sector from tiny floating-point excess.
        let mut sector_count = total_duration_sec.round().max(1.0) as usize;
        sector_count = sector_count.min(total_f.max(1));
        let base = total_f / sector_count;
        let extra = total_f % sector_count;
        let mut sec_counts = Vec::with_capacity(sector_count);
        for si in 0..sector_count {
            let nf = base + if si < extra { 1 } else { 0 };
            if nf > 0 {
                sec_counts.push(nf);
            }
        }
        let dp1 = Arc::new(build_decode_plan(bit1, sh1.as_ref())?); let dp2 = Arc::new(build_decode_plan(bit2, sh2.as_ref())?);
        let mut emitted = 0;
        let cor_h_freq_hz = obs_mhz * 1e6;
        let mut cw_ph = CorWriter::create(&o_dir.join(format!("YAMAGU66_YAMAGU66_{}_phasedarray.cor", c_tag)), &CorHeaderConfig { sampling_speed_hz: fs.round() as i32, observing_frequency_hz: cor_h_freq_hz, fft_point: fft_len as i32, number_of_sector_hint: sec_counts.len() as i32, clock_reference_unix_sec: c_unix, source_name: source_name.clone(), source_ra_rad: gdi.as_ref().map(|v| v.ra).unwrap_or(0.0), source_dec_rad: gdi.as_ref().map(|v| v.dec).unwrap_or(0.0) }, CorStation { name: "YAMAGU66", code: b'M', ecef_m: ant1_ecef }, CorStation { name: "YAMAGU66", code: b'M', ecef_m: ant1_ecef })?;
        let mut cw_11 = CorWriter::create(&o_dir.join(format!("YAMAGU32_YAMAGU32_{}_phasedarray.cor", c_tag)), &CorHeaderConfig { sampling_speed_hz: fs.round() as i32, observing_frequency_hz: cor_h_freq_hz, fft_point: fft_len as i32, number_of_sector_hint: sec_counts.len() as i32, clock_reference_unix_sec: c_unix, source_name: source_name.clone(), source_ra_rad: gdi.as_ref().map(|v| v.ra).unwrap_or(0.0), source_dec_rad: gdi.as_ref().map(|v| v.dec).unwrap_or(0.0) }, CorStation { name: a1_name, code: a1_code, ecef_m: ant1_ecef }, CorStation { name: a1_name, code: a1_code, ecef_m: ant1_ecef })?;
        let mut cw_12 = CorWriter::create(&o_dir.join(format!("YAMAGU32_YAMAGU34_{}_phasedarray.cor", c_tag)), &CorHeaderConfig { sampling_speed_hz: fs.round() as i32, observing_frequency_hz: cor_h_freq_hz, fft_point: fft_len as i32, number_of_sector_hint: sec_counts.len() as i32, clock_reference_unix_sec: c_unix, source_name: source_name.clone(), source_ra_rad: gdi.as_ref().map(|v| v.ra).unwrap_or(0.0), source_dec_rad: gdi.as_ref().map(|v| v.dec).unwrap_or(0.0) }, CorStation { name: a1_name, code: a1_code, ecef_m: ant1_ecef }, CorStation { name: a2_name, code: a2_code, ecef_m: ant2_ecef })?;
        let mut cw_22 = CorWriter::create(&o_dir.join(format!("YAMAGU34_YAMAGU34_{}_phasedarray.cor", c_tag)), &CorHeaderConfig { sampling_speed_hz: fs.round() as i32, observing_frequency_hz: cor_h_freq_hz, fft_point: fft_len as i32, number_of_sector_hint: sec_counts.len() as i32, clock_reference_unix_sec: c_unix, source_name: source_name.clone(), source_ra_rad: gdi.as_ref().map(|v| v.ra).unwrap_or(0.0), source_dec_rad: gdi.as_ref().map(|v| v.dec).unwrap_or(0.0) }, CorStation { name: a2_name, code: a2_code, ecef_m: ant2_ecef }, CorStation { name: a2_name, code: a2_code, ecef_m: ant2_ecef })?;
        
        let mut acc_ph_total = vec![0.0; fft_len/2+1];
        let mut acc_11_total = vec![0.0; fft_len/2+1];
        let mut acc_22_total = vec![0.0; fft_len/2+1];
        if let Some(info) = &gdi {
            let p = o_dir.join(format!("YAMAGU66_{}_geometricdelay.txt", c_tag)); let mut gf = File::create(p)?;
            writeln!(gf, "# Geometric delay per second (ant2 - ant1)\n# sec  delay_seconds  samples")?;
            for sec in 0..sec_counts.len() {
                let (_, _, gd, _, _) = geom::calculate_geometric_delay_and_derivatives(ant1_ecef, ant2_ecef, info.ra, info.dec, info.mjd + sec as f64 / 86400.0);
                writeln!(gf, "{:>4}  {:.9e}  {:.3}", sec, gd + args.gico3_correct, (gd + args.gico3_correct) * fs)?;
            }
        }

        let (mut pr1, mut pr2) = (BufReader::new(File::open(&a1p)?), BufReader::new(File::open(&a2p)?));
        pr1.seek(SeekFrom::Start(s1_s as u64 * bit1 as u64 / 8))?; pr2.seek(SeekFrom::Start(s2_s as u64 * bit2 as u64 / 8))?;

        for (si, &nf) in sec_counts.iter().enumerate() {
            let (mut b1, mut b2) = (vec![0u8; nf * bpf1], vec![0u8; nf * bpf2]); 
            read_with_padding(&mut pr1, &mut b1)?; read_with_padding(&mut pr2, &mut b2)?;
            let mut enc = vec![0u8; nf * bpf_o];
            let (batch_ph, batch_11, batch_12, batch_22) = enc.par_chunks_mut(bpf_o).enumerate().map(|(i, out_f)| {
                let t0 = (emitted + i) as f64 * frame_dt;
                let t_mid = t0 + 0.5 * frame_dt;
                let (mut f1, mut f2, mut s1, mut s2) = (vec![0.0; fft_len], vec![0.0; fft_len], vec![Complex::new(0.0, 0.0); fft_len/2+1], vec![Complex::new(0.0, 0.0); fft_len/2+1]);
                decode_block_into_with_plan(&b1[i*bpf1..(i+1)*bpf1], &levels1, fft_len, &dp1, &mut Vec::new(), &mut f1, lsb1).unwrap();
                decode_block_into_with_plan(&b2[i*bpf2..(i+1)*bpf2], &levels2, fft_len, &dp2, &mut Vec::new(), &mut f2, lsb2).unwrap();
                let (tau1, tau2) = if let Some((ra, dec, mjd0)) = exact_geom_params {
                    let mjd_t = mjd0 + t_mid / 86400.0;
                    let (_, _, gd_t, _, _) = geom::calculate_geometric_delay_and_derivatives(ant1_ecef, ant2_ecef, ra, dec, mjd_t);
                    let net_d_rel_no_clock_t =
                        gd_t + args.gico3_correct + coarse_delay_s + args.delay / fs + extra_delay_rate_sps * t_mid;
                    let clock1_t = clock1_delay_s + clock1_rate_sps * t_mid;
                    let clock2_t = clock2_delay_s + clock2_rate_sps * t_mid;
                    ((net_d_rel_no_clock_t - clock1_t) - d_seek, -clock2_t)
                } else {
                    let net_d1 = net_d1_base;
                    let net_d2 = net_d2_base;
                    let total_rate1 = total_rate1_base;
                    let total_rate2 = total_rate2_base;
                    (
                        delay_seconds_at_time(net_d1, total_rate1, total_accel1_base, t_mid),
                        delay_seconds_at_time(net_d2, total_rate2, total_accel2_base, t_mid),
                    )
                };
                let (int_shift1, frac_delay1) = split_delay_to_integer_and_fractional(tau1, fs);
                let (int_shift2, frac_delay2) = split_delay_to_integer_and_fractional(tau2, fs);
                apply_integer_sample_shift_zerofill(&mut f1, int_shift1);
                apply_integer_sample_shift_zerofill(&mut f2, int_shift2);
                helper.forward_r2c_process(&mut f1, &mut s1).unwrap();
                helper.forward_r2c_process(&mut f2, &mut s2).unwrap();
                // Keep the same reference as cross-correlation path.
                let fr_lo1 = carrier_phase_from_delay(lo1_hz, tau1);
                let fr_lo2 = carrier_phase_from_delay(lo2_hz, tau2);
                apply_delay_and_rate_regular_bins(&mut s1, fft_len, fs/fft_len as f64, frac_delay1, 0.0, 0.0, 0.0, false);
                apply_delay_and_rate_regular_bins(&mut s2, fft_len, fs/fft_len as f64, frac_delay2, 0.0, 0.0, 0.0, false);
                
                let mut cb = vec![Complex::new(0.0, 0.0); fft_len/2+1];
                let mut s1_aligned = vec![Complex::new(0.0, 0.0); fft_len/2+1];
                let mut s2_aligned = vec![Complex::new(0.0, 0.0); fft_len/2+1];
                match out_grid {
                    OutputGrid::Ant1 => {
                        for k in 0..(ba.a1e - ba.a1s) {
                            let i1 = ba.a1s + k;
                            let i2 = ba.a2s + k;
                            s2_aligned[i1] = s2[i2] * fr_lo2;
                        }
                        for k in 0..cb.len() { cb[k] = (s1[k] * fr_lo1) * w1 + s2_aligned[k] * w2; }
                    }
                    OutputGrid::Ant2 => {
                        for k in 0..(ba.a1e - ba.a1s) {
                            let i1 = ba.a1s + k;
                            let i2 = ba.a2s + k;
                            s1_aligned[i2] = s1[i1] * fr_lo1;
                        }
                        for k in 0..cb.len() { cb[k] = s1_aligned[k] * w1 + (s2[k] * fr_lo2) * w2; }
                    }
                }
                let phased_pow = cb.iter().map(|c| c.norm_sqr()).collect::<Vec<_>>();
                cb[0].im = 0.0; if fft_len % 2 == 0 { cb[fft_len/2].im = 0.0; }
                let mut out_t = vec![0.0; fft_len]; helper.inverse_c2r_process(&mut cb, &mut out_t).unwrap();
                if output_lsb {
                    // Convert internal USB-domain waveform back to LSB for output raw format.
                    for odd in out_t.iter_mut().skip(1).step_by(2) { *odd = -*odd; }
                }
                let mut tmp_enc = Vec::new(); quantise_frame(&out_t, bit_out, &levels1, sh1.as_ref(), &mut tmp_enc).unwrap();
                out_f.copy_from_slice(&tmp_enc);
                // Cross-spectrum convention for .cor/frinZ compatibility:
                // X12 = X1(corrected/aligned) * conj(X2(corrected/aligned)).
                match out_grid {
                    OutputGrid::Ant1 => {
                        let s1c: Vec<Complex<f64>> = s1.iter().map(|z| *z * fr_lo1).collect();
                        (
                            phased_pow,
                            s1c.iter().map(|c| c.norm_sqr()).collect::<Vec<_>>(),
                            s1c.iter().zip(s2_aligned.iter()).map(|(z1, z2)| *z1 * z2.conj()).collect::<Vec<_>>(),
                            s2_aligned.iter().map(|c| c.norm_sqr()).collect::<Vec<_>>(),
                        )
                    }
                    OutputGrid::Ant2 => {
                        let s2c: Vec<Complex<f64>> = s2.iter().map(|z| *z * fr_lo2).collect();
                        (
                            phased_pow,
                            s1_aligned.iter().map(|c| c.norm_sqr()).collect::<Vec<_>>(),
                            s1_aligned.iter().zip(s2c.iter()).map(|(z1, z2)| *z1 * z2.conj()).collect::<Vec<_>>(),
                            s2c.iter().map(|c| c.norm_sqr()).collect::<Vec<_>>(),
                        )
                    }
                }
            }).fold(|| (vec![0.0; fft_len/2+1], vec![0.0; fft_len/2+1], vec![Complex::new(0.0, 0.0); fft_len/2+1], vec![0.0; fft_len/2+1]), |mut acc, (p_ph, p_11, p_12, p_22)| {
                for k in 0..acc.0.len() { acc.0[k]+=p_ph[k]; acc.1[k]+=p_11[k]; acc.2[k]+=p_12[k]; acc.3[k]+=p_22[k]; } acc
            }).reduce(|| (vec![0.0; fft_len/2+1], vec![0.0; fft_len/2+1], vec![Complex::new(0.0, 0.0); fft_len/2+1], vec![0.0; fft_len/2+1]), |mut acc1, acc2| {
                for k in 0..acc1.0.len() { acc1.0[k]+=acc2.0[k]; acc1.1[k]+=acc2.1[k]; acc1.2[k]+=acc2.2[k]; acc1.3[k]+=acc2.3[k]; } acc1
            });
            wr.write_all(&enc)?; emitted += nf; println!("[info] Synthesised second {}/{} ({:.2}%)", si+1, sec_counts.len(), (emitted as f64 / total_f as f64) * 100.0);
            for k in 0..acc_ph_total.len() { acc_ph_total[k] += batch_ph[k]; acc_11_total[k] += batch_11[k]; acc_22_total[k] += batch_22[k]; }
            let inv = 1.0 / (nf as f64 * (fft_len as f64).powi(2));
            let s_ph: Vec<Complex<f32>> = batch_ph.iter().take(fft_len/2).map(|&v| Complex::new((v*inv) as f32, 0.0)).collect();
            let s_11: Vec<Complex<f32>> = batch_11.iter().take(fft_len/2).map(|&v| Complex::new((v*inv) as f32, 0.0)).collect();
            let s_22: Vec<Complex<f32>> = batch_22.iter().take(fft_len/2).map(|&v| Complex::new((v*inv) as f32, 0.0)).collect();
            let s_12: Vec<Complex<f32>> = batch_12.iter().take(fft_len/2).map(|v| Complex::new((v.re*inv) as f32, (v.im*inv) as f32)).collect();
            cw_ph.write_sector(c_unix + si as i64, (nf as f64 * fft_len as f64 / fs) as f32, &s_ph)?;
            cw_11.write_sector(c_unix + si as i64, (nf as f64 * fft_len as f64 / fs) as f32, &s_11)?;
            cw_12.write_sector(c_unix + si as i64, (nf as f64 * fft_len as f64 / fs) as f32, &s_12)?;
            cw_22.write_sector(c_unix + si as i64, (nf as f64 * fft_len as f64 / fs) as f32, &s_22)?;
        }
        wr.flush()?; cw_ph.finalize()?; cw_11.finalize()?; cw_12.finalize()?; cw_22.finalize()?;

        println!("[info] Generating phased-array plots...");
        let power_norm = (emitted as f64 * fft_len as f64).max(1.0);
        let phased_auto_mag = finalize_auto_spectrum(&mut acc_ph_total, power_norm)?;
        let a11_auto_mag = finalize_auto_spectrum(&mut acc_11_total, power_norm)?;
        let a22_auto_mag = finalize_auto_spectrum(&mut acc_22_total, power_norm)?;
        // Plotting convention: use [0 .. fft/2-1] for even fft (exclude Nyquist bin).
        let spec_bins = if fft_len % 2 == 0 {
            fft_len / 2
        } else {
            phased_auto_mag.len()
        };
        let phased_plot = &phased_auto_mag[..spec_bins];
        let a11_plot = &a11_auto_mag[..spec_bins];
        let a22_plot = &a22_auto_mag[..spec_bins];

        let df = fs / fft_len as f64 / 1e6;
        let freqs_obs_mhz: Vec<f64> = (0..spec_bins).map(|i| obs_mhz + (i as f64 * df)).collect();

        plot_series_f64_x(&freqs_obs_mhz, phased_plot, "Phased Auto-Spectrum (ObsRef)", &o_dir.join(format!("YAMAGU66_{}_phased_auto_spectrum.png", c_tag)).to_string_lossy(), "Frequency (MHz)", "Power", None, "Auto-Spectrum")?;
        
        let amp_ph: Vec<f64> = phased_plot.iter().map(|&v| v.sqrt()).collect();
        let amp_11: Vec<f64> = a11_plot.iter().map(|&v| v.sqrt()).collect();
        let amp_22: Vec<f64> = a22_plot.iter().map(|&v| v.sqrt()).collect();
        let l1 = format!("{a1_name} (Ref)");
        let l2 = format!("{a2_name} (Target)");
        plot_multi_series_f64_x(&freqs_obs_mhz, &[(&amp_ph, &BLUE, "YAMAGU66 (Phased)"), (&amp_11, &GREEN, &l1), (&amp_22, &RED, &l2)], "Phased Spectrum Amplitude (ObsRef)", &o_dir.join(format!("YAMAGU66_{}_phased_spectrum_amplitude.png", c_tag)).to_string_lossy(), "Frequency (MHz)", "Amplitude", None)?;

        let mut full_spec = vec![Complex::new(0.0, 0.0); fft_len];
        for (i, &v) in phased_plot.iter().enumerate() {
            full_spec[i] = Complex::new(v, 0.0);
            if i > 0 && i < fft_len/2 { full_spec[fft_len - i] = Complex::new(v, 0.0); }
        }
        helper.inverse_c2c(&mut full_spec)?;
        let acf_mag: Vec<f64> = full_spec.iter().map(|c| c.re).collect();
        let acf_shifted: Vec<f64> = acf_mag.iter().cycle().skip(fft_len/2).take(fft_len).copied().collect();
        let lags: Vec<i32> = (-(fft_len as i32 / 2)..(fft_len as i32 / 2)).collect();
        plot_series_with_x(&lags, &[(&acf_shifted, &BLUE)], "Phased Autocorrelation", &o_dir.join(format!("YAMAGU66_{}_phased_autocorrelation.png", c_tag)).to_string_lossy(), "Lag (samples)", "ACF", None, None)?;
    }
    Ok(())
}
