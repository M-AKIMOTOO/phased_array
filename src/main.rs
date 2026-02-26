mod acf;
mod args;
mod cor;
mod geom;
mod ifile;
mod plot;
mod utils;
mod xcf;

use std::collections::{HashSet, VecDeque};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, Write};
#[cfg(target_family = "unix")]
use std::os::fd::AsRawFd;
use std::path::{Path, PathBuf};
use std::sync::{mpsc, Arc};
use std::thread;

use clap::{CommandFactory, Parser};
use libc;
use num_complex::Complex;
use rayon::prelude::*;

use acf::finalize_auto_spectrum;
use args::{
    parse_levels, parse_shuffle, resolve_weight, DEFAULT_LEVELS_2BIT_CSV, DEFAULT_SHUFFLE_IN,
};
use cor::{
    epoch_to_yyyydddhhmmss, CorClockModel, CorHeaderConfig, CorSectorModel, CorStation, CorWriter,
};
use plot::{plot_multi_series_f64_x, plot_series_f64_x, BLUE, GREEN, RED};
use utils::{
    accumulate_power_add, apply_delay_and_rate_regular_bins, build_decode_plan,
    decode_block_into_with_plan, quantise_frame, rotate_regular_bins_with_step, DynError,
    FftHelper,
};
use xcf::finalize_cross_spectrum;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum DelayReference {
    Ant1,
    Ant2,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RotationShiftTarget {
    Ant1,
    Ant2,
}

fn carrier_phase_rotation(
    observing_frequency_hz: f64,
    delay_seconds: f64,
    delay_rate_sps: f64,
    delay_accel_sps2: f64,
    frame_time_s: f64,
) -> Complex<f64> {
    if observing_frequency_hz == 0.0 {
        return Complex::new(1.0, 0.0);
    }
    let tau_t = delay_seconds + delay_rate_sps * frame_time_s + 0.5 * delay_accel_sps2 * frame_time_s.powi(2);
    let phase = -2.0 * std::f64::consts::PI * observing_frequency_hz * tau_t;
    Complex::from_polar(1.0, phase)
}

struct PhaseRotationRecurrence {
    current: Complex<f64>,
    step: Complex<f64>,
    step2: Complex<f64>,
}

impl PhaseRotationRecurrence {
    fn new(
        observing_frequency_hz: f64,
        delay_seconds: f64,
        delay_rate_sps: f64,
        delay_accel_sps2: f64,
        start_time_s: f64,
        delta_time_s: f64,
    ) -> Self {
        if observing_frequency_hz == 0.0 {
            return Self {
                current: Complex::new(1.0, 0.0),
                step: Complex::new(1.0, 0.0),
                step2: Complex::new(1.0, 0.0),
            };
        }

        let tau0 = delay_seconds
            + delay_rate_sps * start_time_s
            + 0.5 * delay_accel_sps2 * start_time_s * start_time_s;
        let phase0 = -2.0 * std::f64::consts::PI * observing_frequency_hz * tau0;

        let delta_tau0 = delay_rate_sps * delta_time_s
            + delay_accel_sps2 * start_time_s * delta_time_s
            + 0.5 * delay_accel_sps2 * delta_time_s * delta_time_s;
        let delta_phase0 = -2.0 * std::f64::consts::PI * observing_frequency_hz * delta_tau0;

        let delta2_tau = delay_accel_sps2 * delta_time_s * delta_time_s;
        let delta2_phase = -2.0 * std::f64::consts::PI * observing_frequency_hz * delta2_tau;

        Self {
            current: Complex::from_polar(1.0, phase0),
            step: Complex::from_polar(1.0, delta_phase0),
            step2: Complex::from_polar(1.0, delta2_phase),
        }
    }

    #[inline(always)]
    fn current(&self) -> Complex<f64> {
        self.current
    }

    #[inline(always)]
    fn advance(&mut self) {
        self.current *= self.step;
        self.step *= self.step2;
    }
}

#[derive(Clone, Copy, Debug)]
struct BandAlignment {
    shift_bins: isize,
    ant1_start: usize,
    ant1_end: usize,
    ant2_start: usize,
    ant2_end: usize,
    residual_shift_hz: f64,
}

fn compute_band_alignment(
    fft_len: usize,
    sampling_hz: f64,
    ant1_center_mhz: f64,
    ant2_center_mhz: f64,
    ant1_bw_mhz: f64,
    ant2_bw_mhz: f64,
) -> Result<BandAlignment, DynError> {
    if fft_len == 0 {
        return Err("FFT length must be positive for band alignment".into());
    }
    if sampling_hz <= 0.0 {
        return Err("Sampling rate must be positive for band alignment".into());
    }
    if ant1_bw_mhz <= 0.0 || ant2_bw_mhz <= 0.0 {
        return Err("Bandwidth must be positive for band alignment".into());
    }

    let half_spec_len = fft_len / 2 + 1;
    let usable_half_bins = if fft_len % 2 == 0 {
        // Keep Nyquist out of band-mapped processing; it has no conjugate pair.
        half_spec_len.saturating_sub(1)
    } else {
        half_spec_len
    };
    let df_hz = sampling_hz / fft_len as f64;
    let valid1 = (((ant1_bw_mhz * 1_000_000.0) / df_hz).floor() as usize)
        .saturating_add(1)
        .min(usable_half_bins);
    let valid2 = (((ant2_bw_mhz * 1_000_000.0) / df_hz).floor() as usize)
        .saturating_add(1)
        .min(usable_half_bins);
    if valid1 == 0 || valid2 == 0 {
        return Err("No valid bins remain after bandwidth limits".into());
    }

    // Align ant2 onto ant1 frequency grid.
    // Shift is defined only by center-frequency difference so different bandwidths
    // do not spuriously introduce extra bin shifts.
    let center1_hz = ant1_center_mhz * 1_000_000.0;
    let center2_hz = ant2_center_mhz * 1_000_000.0;
    let shift_real = (center1_hz - center2_hz) / df_hz;
    let shift_bins = shift_real.round() as isize;
    let residual_shift_hz = (shift_real - shift_bins as f64) * df_hz;

    let ant1_start_i = 0isize.max(-shift_bins);
    let ant1_end_i = (valid1 as isize).min(valid2 as isize - shift_bins);
    if ant1_end_i <= ant1_start_i {
        return Err("No overlapping frequency bins between ant1 and ant2".into());
    }
    let ant2_start_i = ant1_start_i + shift_bins;
    let ant2_end_i = ant1_end_i + shift_bins;
    if ant2_start_i < 0 || ant2_end_i > valid2 as isize {
        return Err("Band alignment produced out-of-range ant2 bins".into());
    }

    Ok(BandAlignment {
        shift_bins,
        ant1_start: ant1_start_i as usize,
        ant1_end: ant1_end_i as usize,
        ant2_start: ant2_start_i as usize,
        ant2_end: ant2_end_i as usize,
        residual_shift_hz,
    })
}

fn format_bit_codes(bit_depth: usize) -> String {
    if bit_depth == 0 {
        return "n/a".to_string();
    }
    let Some(code_count) = 1usize.checked_shl(bit_depth as u32) else {
        return "n/a".to_string();
    };
    if code_count == 0 {
        return "n/a".to_string();
    }
    if code_count > 1024 {
        return format!("too-many-codes({code_count})");
    }
    (0..code_count)
        .map(|code| format!("{code:0bit_depth$b}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn format_level_map(bit_depth: usize, levels: &[f64]) -> String {
    if bit_depth == 0 {
        return "n/a".to_string();
    }
    let Some(code_count) = 1usize.checked_shl(bit_depth as u32) else {
        return "n/a".to_string();
    };
    let limit = code_count.min(levels.len());
    if limit == 0 {
        return "n/a".to_string();
    }
    if limit > 256 {
        return format!("too-many-entries({limit})");
    }
    (0..limit)
        .map(|code| format!("{code:0bit_depth$b}->{:.6}", levels[code]))
        .collect::<Vec<_>>()
        .join(", ")
}

fn seek_forward_samples(
    reader: &mut BufReader<File>,
    samples_to_seek: u64,
    bit_depth: usize,
    file_len_bytes: u64,
    total_samples_sought: &mut isize,
    antenna_label: &str,
    reason: &str,
) -> Result<(), DynError> {
    if samples_to_seek == 0 {
        return Ok(());
    }

    let bytes_to_seek = (samples_to_seek * bit_depth as u64 + 7) / 8;
    if bytes_to_seek >= file_len_bytes {
        return Err(
            format!(
                "Requested seek of {samples_to_seek} samples ({bytes_to_seek} bytes) for {reason} exceeds file length for {antenna_label}",
            )
            .into(),
        );
    }

    let actual_samples_seeked = ((bytes_to_seek * 8) / bit_depth as u64) as isize;
    reader.seek(std::io::SeekFrom::Current(bytes_to_seek as i64))?;
    *total_samples_sought += actual_samples_seeked;
    println!(
        "[info] Seeking {antenna_label} forward by {} bytes (approx. {} samples) for {}",
        bytes_to_seek, actual_samples_seeked, reason
    );
    Ok(())
}

fn apply_integer_delay_with_forward_seek(
    desired_delay_samples: i64,
    bit_depth: usize,
    reason: &str,
    delay_reference: DelayReference,
    reader_ant1: &mut BufReader<File>,
    reader_ant2: &mut BufReader<File>,
    ant1_len_bytes: u64,
    ant2_len_bytes: u64,
    total_samples_sought_ant1: &mut isize,
    total_samples_sought_ant2: &mut isize,
) -> Result<(), DynError> {
    if desired_delay_samples == 0 {
        return Ok(());
    }

    match delay_reference {
        DelayReference::Ant2 => {
            if desired_delay_samples > 0 {
                seek_forward_samples(
                    reader_ant2,
                    desired_delay_samples as u64,
                    bit_depth,
                    ant2_len_bytes,
                    total_samples_sought_ant2,
                    "ant2",
                    reason,
                )?;
            } else {
                let samples_to_seek = (-desired_delay_samples) as u64;
                seek_forward_samples(
                    reader_ant1,
                    samples_to_seek,
                    bit_depth,
                    ant1_len_bytes,
                    total_samples_sought_ant1,
                    "ant1",
                    &format!(
                        "{} (advancing ant1 to compensate negative delay on ant2)",
                        reason
                    ),
                )?;
            }
        }
        DelayReference::Ant1 => {
            if desired_delay_samples > 0 {
                seek_forward_samples(
                    reader_ant1,
                    desired_delay_samples as u64,
                    bit_depth,
                    ant1_len_bytes,
                    total_samples_sought_ant1,
                    "ant1",
                    reason,
                )?;
            } else {
                let samples_to_seek = (-desired_delay_samples) as u64;
                seek_forward_samples(
                    reader_ant2,
                    samples_to_seek,
                    bit_depth,
                    ant2_len_bytes,
                    total_samples_sought_ant2,
                    "ant2",
                    &format!(
                        "{} (advancing ant2 to compensate negative delay on ant1)",
                        reason
                    ),
                )?;
            }
        }
    }

    Ok(())
}

#[cfg(target_family = "unix")]
fn advise_file_sequential(file: &File) {
    let fd = file.as_raw_fd();
    unsafe {
        let _ = libc::posix_fadvise(fd, 0, 0, libc::POSIX_FADV_SEQUENTIAL);
    }
}

#[cfg(not(target_family = "unix"))]
fn advise_file_sequential(_file: &File) {}

#[cfg(target_family = "unix")]
fn drop_file_cache_range(file: &File, offset: u64, len: u64) {
    let Ok(offset_i64) = i64::try_from(offset) else {
        return;
    };
    let Ok(len_i64) = i64::try_from(len) else {
        return;
    };
    let fd = file.as_raw_fd();
    unsafe {
        let _ = libc::posix_fadvise(
            fd,
            offset_i64 as libc::off_t,
            len_i64 as libc::off_t,
            libc::POSIX_FADV_DONTNEED,
        );
    }
}

#[cfg(not(target_family = "unix"))]
fn drop_file_cache_range(_file: &File, _offset: u64, _len: u64) {}

#[cfg(target_family = "unix")]
fn drop_file_cache_all(file: &File) {
    let fd = file.as_raw_fd();
    unsafe {
        let _ = libc::posix_fadvise(fd, 0, 0, libc::POSIX_FADV_DONTNEED);
    }
}

#[cfg(not(target_family = "unix"))]
fn drop_file_cache_all(_file: &File) {}

fn read_block_partial(
    reader: &mut BufReader<File>,
    buffer: &mut Vec<u8>,
) -> Result<usize, DynError> {
    use std::io::ErrorKind;

    let mut total_read = 0usize;
    while total_read < buffer.len() {
        match reader.read(&mut buffer[total_read..]) {
            Ok(0) => break,
            Ok(n) => total_read += n,
            Err(ref e) if e.kind() == ErrorKind::Interrupted => continue,
            Err(e) => return Err(e.into()),
        }
    }
    buffer.truncate(total_read);
    Ok(total_read)
}

struct SecondRawChunk {
    second_index: usize,
    start_frame: usize,
    frames: usize,
    raw1: Vec<u8>,
    raw2: Vec<u8>,
}

fn read_second_chunk(
    reader1: &mut BufReader<File>,
    reader2: &mut BufReader<File>,
    frames_requested: usize,
    bytes_per_frame: usize,
    second_index: usize,
    start_frame: usize,
) -> Result<Option<SecondRawChunk>, DynError> {
    if frames_requested == 0 {
        return Ok(None);
    }

    let mut raw1 = vec![0u8; frames_requested * bytes_per_frame];
    let mut raw2 = vec![0u8; frames_requested * bytes_per_frame];
    let bytes_read1 = read_block_partial(reader1, &mut raw1)?;
    let bytes_read2 = read_block_partial(reader2, &mut raw2)?;
    let frames_read = (bytes_read1 / bytes_per_frame).min(bytes_read2 / bytes_per_frame);
    if frames_read == 0 {
        return Ok(None);
    }

    raw1.truncate(frames_read * bytes_per_frame);
    raw2.truncate(frames_read * bytes_per_frame);
    Ok(Some(SecondRawChunk {
        second_index,
        start_frame,
        frames: frames_read,
        raw1,
        raw2,
    }))
}

#[allow(dead_code)]
const K_BOLTZMANN: f64 = 1.380_649e-23; // J/K
#[allow(dead_code)]
const SI_TO_JY: f64 = 1.0e26; // multiply W/m^2/Hz to express in Jy
const DEFAULT_SAMPLING_HZ: f64 = 1_024_000_000.0;
const DEFAULT_SAMPLING_MHZ: f64 = DEFAULT_SAMPLING_HZ / 1_000_000.0;
const DEFAULT_BIT_DEPTH: usize = 2;
const YAMAGUCHI_FIXED_COARSE_DELAY_S: f64 = 1.7e-6;
const YAMAGUCHI_FIXED_RATE_SPS: f64 = 0.0;

fn resolve_input_paths(
    args: &args::Args,
    preferred_epoch_input: &str,
    ifile_data: Option<&ifile::IFileData>,
) -> Result<(PathBuf, PathBuf, String, String), DynError> {
    match (&args.ant1, &args.ant2) {
        (Some(ant1), Some(ant2)) => {
            let (_, epoch_tag) = epoch_to_yyyydddhhmmss(preferred_epoch_input)?;
            Ok((
                ant1.clone(),
                ant2.clone(),
                preferred_epoch_input.to_string(),
                epoch_tag,
            ))
        }
        (None, None) => {
            if args.ifile.is_none() {
                return Err(
                    "When --ant1/--ant2 are omitted, --ifile must be provided for input auto-resolution"
                        .into(),
                );
            }
            let mut epoch_candidates = vec![preferred_epoch_input.to_string()];
            if args.epoch.is_none() {
                if let Some(meta) = ifile_data {
                    for epoch in &meta.process_epochs {
                        if !epoch_candidates.iter().any(|e| e == epoch) {
                            epoch_candidates.push(epoch.clone());
                        }
                    }
                }
            }

            let data_dir = args
                .data
                .clone()
                .or_else(|| {
                    args.ifile
                        .as_ref()
                        .and_then(|p| p.parent().map(|pp| pp.to_path_buf()))
                })
                .unwrap_or_else(|| PathBuf::from("."));
            let mut tried_messages: Vec<String> = Vec::new();
            for epoch_candidate in epoch_candidates {
                let (_, epoch_tag) = epoch_to_yyyydddhhmmss(&epoch_candidate)?;
                let mut candidates: Vec<(PathBuf, PathBuf, String)> = Vec::new();
                let mut seen_pairs: HashSet<(PathBuf, PathBuf)> = HashSet::new();

                let mut push_candidate = |ant1_prefix: &str, ant2_prefix: &str, label: &str| {
                    let ant1_prefix = ant1_prefix.trim();
                    let ant2_prefix = ant2_prefix.trim();
                    if ant1_prefix.is_empty() || ant2_prefix.is_empty() {
                        return;
                    }
                let ant1 = data_dir.join(format!("{ant1_prefix}_{epoch_tag}.raw"));
                let ant2 = data_dir.join(format!("{ant2_prefix}_{epoch_tag}.raw"));
                if seen_pairs.insert((ant1.clone(), ant2.clone())) {
                    candidates.push((ant1, ant2, label.to_string()));
                }
                };

                if let Some(meta) = ifile_data {
                    if let (Some(ant1_name), Some(ant2_name)) = (
                        meta.ant1_station_name.as_deref(),
                        meta.ant2_station_name.as_deref(),
                    ) {
                        push_candidate(ant1_name, ant2_name, "xml station names");
                    }
                    if let (Some(ant1_key), Some(ant2_key)) =
                        (meta.ant1_station_key.as_deref(), meta.ant2_station_key.as_deref())
                    {
                        let ant1_key = ant1_key.trim();
                        let ant2_key = ant2_key.trim();
                        let mut key1_variants = vec![ant1_key.to_string()];
                        let mut key2_variants = vec![ant2_key.to_string()];
                        let ant1_key_upper = ant1_key.to_ascii_uppercase();
                        let ant2_key_upper = ant2_key.to_ascii_uppercase();
                        let ant1_key_lower = ant1_key.to_ascii_lowercase();
                        let ant2_key_lower = ant2_key.to_ascii_lowercase();
                        if !key1_variants.contains(&ant1_key_upper) {
                            key1_variants.push(ant1_key_upper);
                        }
                        if !key1_variants.contains(&ant1_key_lower) {
                            key1_variants.push(ant1_key_lower);
                        }
                        if !key2_variants.contains(&ant2_key_upper) {
                            key2_variants.push(ant2_key_upper);
                        }
                        if !key2_variants.contains(&ant2_key_lower) {
                            key2_variants.push(ant2_key_lower);
                        }
                        for key1 in &key1_variants {
                            for key2 in &key2_variants {
                                push_candidate(key1, key2, "xml station keys");
                            }
                        }
                        if let Some(ant1_name) = meta.ant1_station_name.as_deref() {
                            for key2 in &key2_variants {
                                push_candidate(ant1_name, key2, "xml ant1 name + ant2 key");
                            }
                        }
                        if let Some(ant2_name) = meta.ant2_station_name.as_deref() {
                            for key1 in &key1_variants {
                                push_candidate(key1, ant2_name, "xml ant1 key + ant2 name");
                            }
                        }
                    }
                }

                // Fixed fallback for legacy Yamaguchi naming.
                push_candidate("YAMAGU32", "YAMAGU34", "legacy YAMAGU fallback");

                for (ant1, ant2, label) in &candidates {
                    if ant1.exists() && ant2.exists() {
                        println!(
                            "[info] Auto-resolved inputs ({}): {} / {}",
                            label,
                            ant1.display(),
                            ant2.display()
                        );
                        if epoch_candidate != preferred_epoch_input {
                            println!(
                                "[info] Selected process epoch {} (skip missing earlier process files)",
                                epoch_candidate
                            );
                        }
                        return Ok((ant1.clone(), ant2.clone(), epoch_candidate, epoch_tag));
                    }
                }

                let tried_this_epoch = candidates
                    .iter()
                    .map(|(ant1, ant2, label)| {
                        format!(
                            "[epoch {} / {}] {} / {}",
                            epoch_candidate,
                            label,
                            ant1.display(),
                            ant2.display()
                        )
                    })
                    .collect::<Vec<_>>()
                    .join("; ");
                tried_messages.push(tried_this_epoch);
            }
            Err(format!(
                "Auto-resolved input files not found for process epochs. Tried: {}. Use --data or pass --ant1/--ant2 explicitly.",
                tried_messages.join("; ")
            )
            .into())
        }
        _ => Err("Specify both --ant1 and --ant2, or omit both and use --ifile (optionally with --data)".into()),
    }
}

fn resolve_output_layout(args: &args::Args, epoch_tag: &str) -> Result<(PathBuf, PathBuf), DynError> {
    let run_stem = format!("YAMAGU66_{epoch_tag}");
    if let Some(path) = args.output.clone() {
        let output_dir = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| Path::new("."))
            .to_path_buf();
        std::fs::create_dir_all(&output_dir)?;
        return Ok((output_dir, path));
    }
    let output_dir = PathBuf::from("phased_array").join(&run_stem);
    std::fs::create_dir_all(&output_dir)?;
    let output_path = output_dir.join(format!("{run_stem}.raw"));
    Ok((output_dir, output_path))
}

fn main() -> Result<(), DynError> {
    if std::env::args_os().len() == 1 {
        args::Args::command().print_help()?;
        println!();
        return Ok(());
    }

    let mut args = args::Args::parse();
    let ifile_data = if let Some(ifile_path) = &args.ifile {
        Some(ifile::parse_ifile(ifile_path)?)
    } else {
        None
    };
    let cli_fft_is_default = args.fft == args::DEFAULT_FFT;
    if cli_fft_is_default {
        if let Some(ifile_fft) = ifile_data.as_ref().and_then(|d| d.fft) {
            args.fft = ifile_fft;
        }
    }

    // Fixed Yamaguchi mode: YAMAGU32 is reference, delay is applied to YAMAGU34 (ant2).
    let delay_reference = DelayReference::Ant2;
    if !args.delay_reference.eq_ignore_ascii_case("ant2") {
        eprintln!(
            "[warn] --delay-reference is ignored in fixed Yamaguchi mode; forcing ant2 (YAMAGU34)"
        );
    }
    let delay_reference_label = "ant2(YAMAGU34, ref=YAMAGU32)";
    let coarse_delay_s = args.coarse.unwrap_or(YAMAGUCHI_FIXED_COARSE_DELAY_S);

    let available_cores = unsafe { libc::sysconf(libc::_SC_NPROCESSORS_ONLN) } as usize;
    if args.cpu > available_cores {
        return Err(format!(
            "--cpu value ({}) exceeds the number of available cores ({})",
            args.cpu, available_cores
        )
        .into());
    }

    if args.cpu == 0 {
        return Err("--cpu must be at least 1".into());
    }
    if args.fft == 0 {
        return Err("--fft must be at least 1".into());
    }
    if !args.fft.is_power_of_two() {
        return Err("--fft must be a power of two (2^n)".into());
    }

    // --- Geometric Delay Calculation ---
    let mut geometric_delay_s_initial: f64 = 0.0;
    let mut geometric_rate_s_initial: f64 = 0.0;
    let mut geometric_accel_s_initial: f64 = 0.0;
    let mut geometric_delay_geom_only: Option<f64> = None;

    #[derive(Clone, Copy)]
    struct GeomDebugInfo {
        ra_rad: f64,
        dec_rad: f64,
        initial_mjd: f64,
    }
    let mut geom_debug_info: Option<GeomDebugInfo> = None;

    let cli_sideband_is_default = args.sideband.eq_ignore_ascii_case("usb");
    let resolved_sideband = if cli_sideband_is_default {
        ifile_data
            .as_ref()
            .and_then(|d| d.sideband.clone())
            .unwrap_or_else(|| args.sideband.clone())
    } else {
        args.sideband.clone()
    };
    if !resolved_sideband.eq_ignore_ascii_case("usb")
        && !resolved_sideband.eq_ignore_ascii_case("lsb")
    {
        return Err("Resolved sideband must be usb or lsb".into());
    }
    // For USB-sampled data, no alternating-sign conversion is required.
    // Apply conversion only when sampler sideband is LSB.
    let lsb_to_usb = resolved_sideband.eq_ignore_ascii_case("lsb");

    let cli_sampling_hz = args.sampling * 1_000_000.0;
    let cli_sampling_is_default = (args.sampling - DEFAULT_SAMPLING_MHZ).abs() < 1e-12;
    let ifile_sampling_hz = ifile_data.as_ref().and_then(|d| d.sampling_hz);
    let sampling_rate_hz = match (cli_sampling_is_default, ifile_sampling_hz) {
        (true, Some(v)) => v,
        _ => cli_sampling_hz,
    };
    if sampling_rate_hz <= 0.0 {
        return Err("Resolved sampling frequency must be positive".into());
    }

    let ifile_bit = ifile_data.as_ref().and_then(|d| d.bit);
    let bit_depth = match (args.bit == DEFAULT_BIT_DEPTH, ifile_bit) {
        (true, Some(v)) => v,
        _ => args.bit,
    };
    if bit_depth == 0 {
        return Err("Resolved bit depth must be at least 1".into());
    }

    let cli_obsfreq_mhz = args.sky_freq;
    let cli_obsfreq_is_default = cli_obsfreq_mhz.abs() < 1e-12;
    let ifile_obsfreq_mhz = ifile_data.as_ref().and_then(|d| d.obsfreq_mhz);
    let observing_frequency_mhz = match (cli_obsfreq_is_default, ifile_obsfreq_mhz) {
        (true, Some(v)) => v,
        _ => cli_obsfreq_mhz,
    };
    if observing_frequency_mhz < 0.0 {
        return Err("Resolved observing frequency must be non-negative".into());
    }

    let ra_input = args
        .ra
        .clone()
        .or_else(|| ifile_data.as_ref().map(|d| d.ra.clone()));
    let dec_input = args
        .dec
        .clone()
        .or_else(|| ifile_data.as_ref().map(|d| d.dec.clone()));
    let epoch_input = args
        .epoch
        .clone()
        .or_else(|| ifile_data.as_ref().and_then(|d| d.epoch.clone()))
        .unwrap_or_else(|| "2000".to_string());
    let (ant1_path, ant2_path, epoch_input, cor_epoch_tag) =
        resolve_input_paths(&args, &epoch_input, ifile_data.as_ref())?;
    let (cor_epoch_unix_s, _) = epoch_to_yyyydddhhmmss(&epoch_input)?;
    let (output_dir, output_path) = resolve_output_layout(&args, &cor_epoch_tag)?;

    if ra_input.is_some() ^ dec_input.is_some() {
        return Err("Specify both --ra and --dec together".into());
    }

    if let (Some(ra_text), Some(dec_text)) = (ra_input.as_ref(), dec_input.as_ref()) {
        let ra_rad = geom::parse_ra(ra_text)?;
        let dec_rad = geom::parse_dec(dec_text)?;
        let initial_mjd = geom::parse_epoch_to_mjd(&epoch_input)?;

        let (_delay1, _delay2, geom_delay, geom_rate, geom_accel) =
            geom::calculate_geometric_delay_and_derivatives(
                geom::YAMAGU32_ECEF,
                geom::YAMAGU34_ECEF,
                ra_rad,
                dec_rad,
                initial_mjd,
            );

        geometric_delay_geom_only = Some(geom_delay);
        geom_debug_info = Some(GeomDebugInfo {
            ra_rad,
            dec_rad,
            initial_mjd,
        });
        geometric_delay_s_initial = geom_delay;
        geometric_rate_s_initial = geom_rate;
        geometric_accel_s_initial = geom_accel;
    }

    geometric_delay_s_initial += args.gico3_correct;

    let clock_delay_s = ifile_data
        .as_ref()
        .and_then(|d| d.clock_delay_s)
        .unwrap_or(0.0);
    let clock_rate_sps = ifile_data
        .as_ref()
        .and_then(|d| d.clock_rate_sps)
        .unwrap_or(YAMAGUCHI_FIXED_RATE_SPS);

    // Target inputs are treated as values relative to clock reference.
    let residual_delay_s = args.delay / sampling_rate_hz;
    let net_delay_s_diff =
        geometric_delay_s_initial + coarse_delay_s + clock_delay_s + residual_delay_s;

    let geometric_delay_samples_float = geometric_delay_s_initial * sampling_rate_hz;
    let geometric_integer_delay_samples = geometric_delay_samples_float.round() as isize;

    let coarse_delay_samples_float = coarse_delay_s * sampling_rate_hz;
    let clock_delay_samples_float = clock_delay_s * sampling_rate_hz;
    let residual_delay_samples_float = args.delay;
    let coarse_delay_applied_samples = if delay_reference == DelayReference::Ant2 {
        coarse_delay_samples_float
    } else {
        -coarse_delay_samples_float
    };
    let clock_delay_applied_samples = if delay_reference == DelayReference::Ant2 {
        clock_delay_samples_float
    } else {
        -clock_delay_samples_float
    };
    let residual_delay_applied_samples = if delay_reference == DelayReference::Ant2 {
        residual_delay_samples_float
    } else {
        -residual_delay_samples_float
    };

    let coarse_res_delay_samples_float =
        coarse_delay_samples_float + clock_delay_samples_float + residual_delay_samples_float;
    let coarse_res_integer_delay_samples = coarse_res_delay_samples_float.round() as isize;

    let mut total_samples_sought_ant1: isize = 0;
    let mut total_samples_sought_ant2: isize = 0;

    let cli_level_is_default = args
        .level
        .as_ref()
        .map(|v| v.replace(' ', "") == DEFAULT_LEVELS_2BIT_CSV)
        .unwrap_or(true);
    let resolved_level = if cli_level_is_default {
        ifile_data
            .as_ref()
            .and_then(|d| d.level.clone())
            .or_else(|| args.level.clone())
    } else {
        args.level.clone()
    };
    let levels = Arc::new(parse_levels(bit_depth, resolved_level)?);
    let resolved_shuffle = args
        .shuffle_in
        .clone()
        .or_else(|| ifile_data.as_ref().and_then(|d| d.shuffle.clone()));
    let shuffle_in = Arc::new(parse_shuffle(resolved_shuffle, &DEFAULT_SHUFFLE_IN)?);
    let shuffle_display_external: Vec<usize> = (0..32).map(|i| shuffle_in[31 - i]).collect();

    let compute_threads = if args.debug_corr && args.cpu > 1 {
        args.cpu - 1
    } else {
        args.cpu
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(compute_threads)
        .build_global()
        .map_err(|_| "Failed to initialise rayon thread pool")?;

    let seconds_per_frame = args.fft as f64 / sampling_rate_hz;
    // 1 秒あたりのフレーム数を推算し、I/O と FFT を秒単位でまとめて処理する
    let frames_per_second_f64 = if seconds_per_frame > 0.0 {
        1.0 / seconds_per_frame
    } else {
        sampling_rate_hz / args.fft as f64
    };
    let mut frames_per_second = frames_per_second_f64.round() as usize;
    if frames_per_second == 0 {
        frames_per_second = 1;
    }
    if (frames_per_second_f64 - frames_per_second as f64).abs() > 1e-6 {
        println!(
            "[warn] Non-integer frames-per-second {:.6}; using nearest integer {}",
            frames_per_second_f64, frames_per_second
        );
    }
    let frames_per_second = frames_per_second;
    let sampling_mhz = sampling_rate_hz / 1_000_000.0;
    let bandwidth_mhz = sampling_mhz / 2.0;
    let ant1_rotation_hz = ifile_data
        .as_ref()
        .and_then(|d| d.ant1_rotation_hz)
        .unwrap_or(0.0);
    let ant2_rotation_hz = ifile_data
        .as_ref()
        .and_then(|d| d.ant2_rotation_hz)
        .unwrap_or(0.0);
    let ant1_bw_mhz = ifile_data
        .as_ref()
        .and_then(|d| d.ant1_bw_mhz)
        .unwrap_or(bandwidth_mhz);
    let ant2_bw_mhz = ifile_data
        .as_ref()
        .and_then(|d| d.ant2_bw_mhz)
        .unwrap_or(bandwidth_mhz);
    // sky-freq / XML <frequency> is treated as obsfreq center [MHz].
    let ant1_obsfreq_mhz = observing_frequency_mhz + ant1_rotation_hz / 1_000_000.0;
    let ant2_obsfreq_mhz = observing_frequency_mhz + ant2_rotation_hz / 1_000_000.0;
    let ant1_center_mhz = ifile_data
        .as_ref()
        .and_then(|d| d.ant1_center_mhz)
        .unwrap_or(ant1_obsfreq_mhz);
    let ant2_center_mhz = ifile_data
        .as_ref()
        .and_then(|d| d.ant2_center_mhz)
        .unwrap_or(ant2_obsfreq_mhz);
    let ant1_low_mhz = ant1_center_mhz - 0.5 * ant1_bw_mhz;
    let ant2_low_mhz = ant2_center_mhz - 0.5 * ant2_bw_mhz;
    let ant1_high_mhz = ant1_center_mhz + 0.5 * ant1_bw_mhz;
    let ant2_high_mhz = ant2_center_mhz + 0.5 * ant2_bw_mhz;
    let ant1_observing_frequency_hz = ant1_low_mhz * 1_000_000.0;
    let ant2_observing_frequency_hz = ant2_low_mhz * 1_000_000.0;
    if ant1_observing_frequency_hz <= 0.0 || ant2_observing_frequency_hz <= 0.0 {
        return Err("Resolved per-antenna observing frequency must be positive".into());
    }
    // Re-resolve target correction frequency now that low-edge is determined.
    let phase_correction_frequency_hz = match delay_reference {
        DelayReference::Ant2 => ant1_observing_frequency_hz,
        DelayReference::Ant1 => ant2_observing_frequency_hz,
    };
    // Recalculate target_rate_hz in s/s using the correct reference frequency.
    let user_rate_sps = if args.rate.abs() > 0.0 {
        args.rate / phase_correction_frequency_hz
    } else {
        0.0
    };
    let fixed_rate_delay_s_per_s = clock_rate_sps + user_rate_sps;

    let band_alignment = compute_band_alignment(
        args.fft,
        sampling_rate_hz,
        ant1_center_mhz,
        ant2_center_mhz,
        ant1_bw_mhz,
        ant2_bw_mhz,
    )?;
    const ROTATION_EPS_HZ: f64 = 1e-6;
    let ant1_has_rotation = ant1_rotation_hz.abs() > ROTATION_EPS_HZ;
    let ant2_has_rotation = ant2_rotation_hz.abs() > ROTATION_EPS_HZ;
    let rotation_shift_target = if ant1_has_rotation && !ant2_has_rotation {
        RotationShiftTarget::Ant1
    } else {
        RotationShiftTarget::Ant2
    };
    let align_ref_is_ant1 = matches!(rotation_shift_target, RotationShiftTarget::Ant2);
    let align_shift_bins = if align_ref_is_ant1 {
        band_alignment.shift_bins
    } else {
        -band_alignment.shift_bins
    };
    let align_ref_label = if align_ref_is_ant1 { "ant1" } else { "ant2" };
    let align_shifted_label = if align_ref_is_ant1 { "ant2" } else { "ant1" };
    let fft_bin_width_mhz = sampling_mhz / args.fft as f64;
    let overlap_bins = band_alignment.ant1_end.saturating_sub(band_alignment.ant1_start);
    let overlap_bw_mhz = overlap_bins as f64 * fft_bin_width_mhz;
    let overlap_low_mhz = ant1_low_mhz + band_alignment.ant1_start as f64 * fft_bin_width_mhz;
    let overlap_center_mhz = overlap_low_mhz + 0.5 * overlap_bw_mhz;
    // .cor header observing frequency must remain the base/sky frequency
    // from XML <frequency> or --sky-freq/--obsfreq.
    let cor_observing_frequency_hz = observing_frequency_mhz * 1_000_000.0;

    let (weight1, _gain1_effective, sefd1_used, aeff1) = resolve_weight(
        args.tsys1,
        args.gain1,
        args.sefd1,
        args.diameter1,
        args.eta1,
        "Antenna 1",
    )?;
    let (weight2, _gain2_effective, sefd2_used, aeff2) = resolve_weight(
        args.tsys2,
        args.gain2,
        args.sefd2,
        args.diameter2,
        args.eta2,
        "Antenna 2",
    )?;

    let length_display = args
        .length
        .map(|v| format!("{v}"))
        .unwrap_or_else(|| "n/a".to_string());

    let bits_per_frame = args
        .fft
        .checked_mul(bit_depth)
        .ok_or("Overflow computing bits per frame")?;
    let bytes_per_frame = (bits_per_frame + 7) / 8;
    if bytes_per_frame == 0 {
        return Err("Byte count per frame evaluated to zero".into());
    }
    if bits_per_frame % 32 != 0 {
        return Err("Product of --fft and --bit must be divisible by 32".into());
    }

    let file1_meta = std::fs::metadata(&ant1_path)?;
    let file2_meta = std::fs::metadata(&ant2_path)?;
    let frames_file1 = (file1_meta.len() / bytes_per_frame as u64) as usize;
    let frames_file2 = (file2_meta.len() / bytes_per_frame as u64) as usize;
    let mut total_frames_to_process = frames_file1.min(frames_file2);
    if total_frames_to_process == 0 {
        return Err("Insufficient data for a single FFT frame in one of the inputs".into());
    }

    let estimate_seconds_file1 = (frames_file1 * args.fft) as f64 / sampling_rate_hz;
    let estimate_seconds_file2 = (frames_file2 * args.fft) as f64 / sampling_rate_hz;

    println!("Starting phased array processing with the following arguments:");
    println!("--------------------------------------------------");
    println!(
        "  ant1:       {} (Size: {} bytes, Estimated Obs Time: {:.2}s)",
        ant1_path.display(),
        file1_meta.len(),
        estimate_seconds_file1
    );
    println!(
        "  ant2:       {} (Size: {} bytes, Estimated Obs Time: {:.2}s)",
        ant2_path.display(),
        file2_meta.len(),
        estimate_seconds_file2
    );
    if let (Some(ra_text), Some(dec_text)) = (ra_input.as_ref(), dec_input.as_ref()) {
        println!("  ra/dec:     {} / {}", ra_text, dec_text);
        println!("  epoch:      {}", epoch_input);
        if let Some(base_geom_delay) = geometric_delay_geom_only {
            println!("  geom-delay: {:.6e} s (ant2 - ant1)", base_geom_delay);
            println!(
                "  geom-rate:  {:.6e} s/s (ant2 - ant1) => {:.6} Hz @ sky",
                geometric_rate_s_initial,
                geometric_rate_s_initial * phase_correction_frequency_hz
            );
            println!(
                "  geom-accel: {:.6e} s/s^2 (ant2 - ant1)",
                geometric_accel_s_initial
            );
            if args.gico3_correct != 0.0 {
                println!(
                    "  geom-delay gico3-correct: +{:.6e} s => total {:.6e} s",
                    args.gico3_correct, geometric_delay_s_initial
                );
            }
        }
    } else {
        println!("  [warn] ra/dec not set; GEOMETRIC DELAY CORRECTION IS DISABLED.");
    }
    println!("  delay-ref:  {}", delay_reference_label);
    let coarse_label = if args.coarse.is_some() {
        "coarse-delay user"
    } else {
        "coarse-delay fixed"
    };
    println!(
        "  {}: {:.6e} s (ant2 - ant1) => applied {:.3} samples to {}",
        coarse_label, coarse_delay_s, coarse_delay_applied_samples, delay_reference_label
    );
    println!(
        "  clock-delay ref: {:.6e} s (ant2 - ant1) => applied {:.3} samples to {}",
        clock_delay_s, clock_delay_applied_samples, delay_reference_label
    );
    println!(
        "  res-delay input (target rel): {} samples (ant2 - ant1) => applied {:.3} samples to {}",
        args.delay, residual_delay_applied_samples, delay_reference_label
    );
    println!(
        "  total-delay (ant2 - ant1): {:.3} samples ({:.3e} s)",
        net_delay_s_diff * sampling_rate_hz,
        net_delay_s_diff
    );
    println!(
        "  delay-rate: {:.6} Hz (clock-rel {:.6} + target-rel {:.6})",
        fixed_rate_delay_s_per_s * phase_correction_frequency_hz,
        clock_rate_sps * phase_correction_frequency_hz,
        args.rate
    );
    println!("  sky-freq:   {:.6} MHz (obsfreq base)", observing_frequency_mhz);
    println!(
        "  rotation:   {:.3} Hz (ant1), {:.3} Hz (ant2)",
        ant1_rotation_hz, ant2_rotation_hz
    );
    println!("  obsfreq:    {:.6} MHz (ant1), {:.6} MHz (ant2)", ant1_center_mhz, ant2_center_mhz);
    println!(
        "  phase-freq: {:.6} MHz (delay/rate correction target)",
        phase_correction_frequency_hz / 1_000_000.0
    );
    println!("  bw:         {:.3} MHz", bandwidth_mhz);
    println!(
        "  sampling:   {:.0} Hz ({:.6} MHz)",
        sampling_rate_hz, sampling_mhz
    );
    println!("  samples/s:  {:.0}", sampling_rate_hz);
    println!(
        "  frames/s:   {:.6} (fft {})",
        frames_per_second_f64, args.fft
    );
    println!(
        "  debug:      {}",
        if args.debug { "true" } else { "false" }
    );
    if args.debug {
        println!("  debug-frames:{}", args.debug_frames);
    }
    println!("  bit:        {} ({})", bit_depth, format_bit_codes(bit_depth));
    println!("  level:      {:?}", levels);
    println!("  level-map:  {}", format_level_map(bit_depth, &levels));
    println!("  shuffle-in: {:?}", shuffle_display_external);
    println!("  sideband:   {}", resolved_sideband);
    println!("  obs-band1:  {:.6} .. {:.6} MHz", ant1_low_mhz, ant1_high_mhz);
    println!("  obs-band2:  {:.6} .. {:.6} MHz", ant2_low_mhz, ant2_high_mhz);
    println!(
        "  band-ant1:  low {:.6} MHz, high {:.6} MHz, center {:.6} MHz, bw {:.6} MHz",
        ant1_low_mhz, ant1_high_mhz, ant1_center_mhz, ant1_bw_mhz
    );
    println!(
        "  band-ant2:  low {:.6} MHz, high {:.6} MHz, center {:.6} MHz, bw {:.6} MHz",
        ant2_low_mhz, ant2_high_mhz, ant2_center_mhz, ant2_bw_mhz
    );
    println!(
        "  band-align: ant2->ant1 shift {} bins, overlap ant1[{}..{}), ant2[{}..{})",
        band_alignment.shift_bins,
        band_alignment.ant1_start,
        band_alignment.ant1_end,
        band_alignment.ant2_start,
        band_alignment.ant2_end
    );
    println!(
        "  band-overlap: {} bins ({:.6} MHz, {:.6} MHz/bin), center {:.6} MHz",
        overlap_bins, overlap_bw_mhz, fft_bin_width_mhz, overlap_center_mhz
    );
    println!(
        "  rotation-shift: target {} (XML rotation), grid {}-ref, shift {} bins",
        align_shifted_label, align_ref_label, align_shift_bins
    );
    if ant1_has_rotation == ant2_has_rotation && ant1_rotation_hz != ant2_rotation_hz {
        println!(
            "  [warn] Both/no station rotation keys are active; using {} as shift target by default",
            align_shifted_label
        );
    }
    if band_alignment.residual_shift_hz.abs() > (sampling_rate_hz / args.fft as f64) * 0.05 {
        println!(
            "  [warn] residual fractional band shift {:.3} Hz (integer-bin alignment applied)",
            band_alignment.residual_shift_hz
        );
    }
    println!(
        "  ant-fixed:  YAMAGU32={:?}, YAMAGU34={:?}",
        geom::YAMAGU32_ECEF,
        geom::YAMAGU34_ECEF
    );
    println!(
        "  cpu:        {} (compute threads: {})",
        args.cpu, compute_threads
    );
    println!("  skip:       0");
    println!("  length:     {}", length_display);
    println!("  tsys:       {} (ant1), {} (ant2)", args.tsys1, args.tsys2);
    if args.diameter1.is_some() || args.diameter2.is_some() {
        println!(
            "  diameter:   {} (ant1), {} (ant2)",
            args.diameter1
                .map_or_else(|| "n/a".to_string(), |d| format!("{d}")),
            args.diameter2
                .map_or_else(|| "n/a".to_string(), |d| format!("{d}"))
        );
    }
    println!("  eta:        {} (ant1), {} (ant2)", args.eta1, args.eta2);
    if aeff1.is_some() || aeff2.is_some() {
        println!(
            "  aeff:       {} (ant1), {} (ant2)",
            aeff1.map_or_else(|| "n/a".to_string(), |a| format!("{a}")),
            aeff2.map_or_else(|| "n/a".to_string(), |a| format!("{a}"))
        );
    }
    if sefd1_used.is_some() || sefd2_used.is_some() {
        println!(
            "  sefd:       {} (ant1), {} (ant2)",
            sefd1_used.map_or_else(|| "n/a".to_string(), |s| format!("{s}")),
            sefd2_used.map_or_else(|| "n/a".to_string(), |s| format!("{s}"))
        );
    }
    println!("  gain:       {} (ant1), {} (ant2)", args.gain1, args.gain2);
    println!("  weight:     {} (ant1), {} (ant2)", weight1, weight2);
    println!("  results:    {}", output_dir.display());
    println!("  output:     {}", output_path.display());
    println!("--------------------------------------------------");

    if let Some(length_seconds) = args.length {
        if length_seconds <= 0.0 {
            return Err("--length must be positive".into());
        }
        let total_samples_by_length = (sampling_rate_hz * length_seconds).floor();
        if !total_samples_by_length.is_finite() {
            return Err("Invalid --length value".into());
        }
        let frames_by_length = (total_samples_by_length as usize) / args.fft;
        if frames_by_length == 0 {
            return Err("Requested length shorter than one FFT frame".into());
        }
        if frames_by_length > total_frames_to_process {
            let max_seconds = (total_frames_to_process * args.fft) as f64 / sampling_rate_hz;
            println!(
                "[warn] --length {:.3}s exceeds available overlap {:.6}s; truncating to available data.",
                length_seconds, max_seconds
            );
        }
        total_frames_to_process = total_frames_to_process.min(frames_by_length);
        // Apply length limit
    }

    let helper = Arc::new(FftHelper::new(args.fft));
    let mut reader1 = BufReader::new(File::open(&ant1_path)?);
    let mut reader2 = BufReader::new(File::open(&ant2_path)?);
    advise_file_sequential(reader1.get_ref());
    advise_file_sequential(reader2.get_ref());

    if geometric_integer_delay_samples != 0 {
        let desired_delay_for_reference = if delay_reference == DelayReference::Ant2 {
            geometric_integer_delay_samples as i64
        } else {
            (-geometric_integer_delay_samples) as i64
        };
        println!(
            "[info] Applying geometric delay (reference {}): {} samples (raw diff {})",
            delay_reference_label, desired_delay_for_reference, geometric_integer_delay_samples
        );
        apply_integer_delay_with_forward_seek(
            desired_delay_for_reference,
            bit_depth,
            "geometric delay",
            delay_reference,
            &mut reader1,
            &mut reader2,
            file1_meta.len(),
            file2_meta.len(),
            &mut total_samples_sought_ant1,
            &mut total_samples_sought_ant2,
        )?;
    }

    if coarse_res_integer_delay_samples != 0 {
        let desired_delay_for_reference = if delay_reference == DelayReference::Ant2 {
            coarse_res_integer_delay_samples as i64
        } else {
            (-coarse_res_integer_delay_samples) as i64
        };
        println!(
            "[info] Applying coarse + res-delay (reference {}): {} samples (raw diff {})",
            delay_reference_label, desired_delay_for_reference, coarse_res_integer_delay_samples
        );
        apply_integer_delay_with_forward_seek(
            desired_delay_for_reference,
            bit_depth,
            "coarse+res delay",
            delay_reference,
            &mut reader1,
            &mut reader2,
            file1_meta.len(),
            file2_meta.len(),
            &mut total_samples_sought_ant1,
            &mut total_samples_sought_ant2,
        )?;
    }

    // Calculate the total delay applied by seeking in seconds
    let delay_applied_by_seek_s =
        (total_samples_sought_ant2 - total_samples_sought_ant1) as f64 / sampling_rate_hz;

    let ant1_pos_bytes = reader1.stream_position()?;
    let ant2_pos_bytes = reader2.stream_position()?;
    let remaining_bytes_ant1 = file1_meta.len().saturating_sub(ant1_pos_bytes);
    let remaining_bytes_ant2 = file2_meta.len().saturating_sub(ant2_pos_bytes);
    let remaining_frames_ant1 = (remaining_bytes_ant1 / bytes_per_frame as u64) as usize;
    let remaining_frames_ant2 = (remaining_bytes_ant2 / bytes_per_frame as u64) as usize;

    let total_frames = remaining_frames_ant1
        .min(remaining_frames_ant2)
        .min(total_frames_to_process);

    if total_frames == 0 {
        return Err("Insufficient overlapping data after coarse delay adjustment".into());
    }
    println!(
        "[info] Delay corrected by seek: {:.3e} s (equivalent to {:.3} samples)",
        delay_applied_by_seek_s,
        delay_applied_by_seek_s * sampling_rate_hz
    );
    let remaining_delay_diff_s = net_delay_s_diff - delay_applied_by_seek_s;
    let remaining_delay_for_reference_s = if delay_reference == DelayReference::Ant2 {
        remaining_delay_diff_s
    } else {
        -remaining_delay_diff_s
    };
    println!(
        "[info] Remaining delay for phase correction: {:.3e} s (equivalent to {:.3} samples, reference {})",
        remaining_delay_for_reference_s,
        remaining_delay_for_reference_s * sampling_rate_hz,
        delay_reference_label
    );

    println!("[info] Processing {} overlapping frames.", total_frames);

    if args.debug_corr {
        if args.fft != 1024 {
            println!(
                "[info] Delay search is typically sufficient with --fft 1024; running requested --fft {}",
                args.fft
            );
        }
        let fft_len = args.fft;
        let half = fft_len / 2;
        let lags: Vec<i32> = (-(half as isize)..half as isize)
            .map(|i| i as i32)
            .collect();
        let half_spec_len = args.fft / 2 + 1;

        let freqs: Vec<f64> = (0..half_spec_len)
            .map(|i| (i as f64 * sampling_rate_hz) / args.fft as f64)
            .collect();
        let freqs_mhz: Vec<f64> = freqs.iter().map(|f| f / 1_000_000.0).collect();
        let fft_bin_step_hz_for_corr = sampling_rate_hz / args.fft as f64;

        // 相関解析の蓄積バッファ（クロスは複素数、自己は実数パワー）
        let mut accumulated_cross_spec = vec![Complex::new(0.0, 0.0); half_spec_len];
        let mut accumulated_auto_spec1 = vec![0.0f64; half_spec_len];
        let mut accumulated_auto_spec2 = vec![0.0f64; half_spec_len];
        let mut processed_frames: usize = 0;

        println!("Correlating...");

        const CHANNEL_BUFFER_SIZE: usize = 4; // Buffer up to 4 I/O chunks
        const MAX_FRAMES_PER_IO_CHUNK: usize = 4096;
        let frames_per_io_chunk = frames_per_second.min(MAX_FRAMES_PER_IO_CHUNK).max(1);
        if frames_per_io_chunk < frames_per_second {
            println!(
                "[info] Capping correlation I/O chunk size to {} frames (from {} frames ~= 1 s).",
                frames_per_io_chunk, frames_per_second
            );
        }
        let (tx, rx) = mpsc::sync_channel::<Result<Vec<Vec<u8>>, String>>(CHANNEL_BUFFER_SIZE);

        // --- I/O Producer Thread ---
        let producer_thread = {
            let total_frames = total_frames;
            let bytes_per_frame = bytes_per_frame;
            let cpu_arg = args.cpu;
            let frames_per_io_chunk = frames_per_io_chunk;
            let reader1_start_pos = reader1.stream_position()?;
            let reader2_start_pos = reader2.stream_position()?;
            thread::spawn(move || {
                if cpu_arg >= 2 {
                    unsafe {
                        let mut cpu_set = std::mem::zeroed();
                        libc::CPU_ZERO(&mut cpu_set);
                        libc::CPU_SET(0, &mut cpu_set); // Pin to Core 0
                        if libc::sched_setaffinity(
                            0,
                            std::mem::size_of::<libc::cpu_set_t>(),
                            &cpu_set,
                        ) != 0
                        {
                            eprintln!("[warn] Failed to set affinity for I/O thread");
                        }
                    }
                }

                let mut frames_read = 0;
                // Read in large chunks to keep the disk busy and feed the consumer
                while frames_read < total_frames {
                    let frames_to_read_this_chunk =
                        (total_frames - frames_read).min(frames_per_io_chunk);
                    if frames_to_read_this_chunk == 0 {
                        break;
                    }
                    let bytes_to_read = frames_to_read_this_chunk * bytes_per_frame;

                    let mut read_buffer1 = vec![0u8; bytes_to_read];
                    if let Err(e) = reader1.read_exact(&mut read_buffer1) {
                        let _ = tx.send(Err(format!(
                            "\n[warn] Failed to read from ant1 file: {}. Stopping correlation.",
                            e
                        )));
                        break;
                    }

                    let mut read_buffer2 = vec![0u8; bytes_to_read];
                    if let Err(e) = reader2.read_exact(&mut read_buffer2) {
                        let _ = tx.send(Err(format!(
                            "\n[warn] Failed to read from ant2 file: {}. Stopping correlation.",
                            e
                        )));
                        break;
                    }

                    if tx.send(Ok(vec![read_buffer1, read_buffer2])).is_err() {
                        // Receiver has hung up, no point in reading more.
                        break;
                    }
                    let chunk_offset = frames_read as u64 * bytes_per_frame as u64;
                    drop_file_cache_range(
                        reader1.get_ref(),
                        reader1_start_pos.saturating_add(chunk_offset),
                        bytes_to_read as u64,
                    );
                    drop_file_cache_range(
                        reader2.get_ref(),
                        reader2_start_pos.saturating_add(chunk_offset),
                        bytes_to_read as u64,
                    );
                    frames_read += frames_to_read_this_chunk;
                }
            })
        };

        // --- CPU Consumer (Main Thread) ---
        if args.cpu >= 2 {
            unsafe {
                let mut cpu_set = std::mem::zeroed();
                libc::CPU_ZERO(&mut cpu_set);
                for i in 1..args.cpu {
                    libc::CPU_SET(i, &mut cpu_set);
                }
                if libc::sched_setaffinity(0, std::mem::size_of::<libc::cpu_set_t>(), &cpu_set) != 0
                {
                    eprintln!("[warn] Failed to set affinity for main/compute thread");
                }
            }
        }

        let levels_for_corr = Arc::new(levels.clone());
        let shuffle_in_for_corr = Arc::new(shuffle_in.clone());
        let helper_for_corr = helper.clone();
        let args_for_corr = Arc::new(args.clone()); // Capture args
        let bit_depth_for_corr = bit_depth;
        let half_spec_len_for_corr = half_spec_len; // Capture half_spec_len
        let observing_frequency_hz_for_corr = phase_correction_frequency_hz;
        let seconds_per_frame_for_corr = seconds_per_frame; // Capture seconds_per_frame
        let delay_applied_by_seek_s_for_corr = delay_applied_by_seek_s; // Capture delay_applied_by_seek_s
        let residual_delay_s_for_corr = residual_delay_s; // Capture residual_delay_s
        let coarse_delay_s_for_corr = coarse_delay_s; // Capture fixed coarse delay
        let clock_delay_s_for_corr = clock_delay_s; // Capture fixed clock delay
        let fixed_rate_delay_sps_for_corr = fixed_rate_delay_s_per_s; // Capture fixed delay-rate term
        let delay_reference_for_corr = delay_reference; // Capture delay reference selection
        let frames_per_batch_for_corr = args.cpu.max(1);
        let bytes_per_batch_for_corr = frames_per_batch_for_corr * bytes_per_frame;
        let inv_fft2_corr = 1.0 / ((args.fft as f64) * (args.fft as f64));
        let decode_plan_for_corr = Arc::new(build_decode_plan(
            bit_depth_for_corr,
            shuffle_in_for_corr.as_ref(),
        )?);
        let overlap_len_for_corr = overlap_bins;
        let ant1_overlap_start_for_corr = band_alignment.ant1_start;
        let ant2_overlap_start_for_corr = band_alignment.ant2_start;
        let mut last_progress_second: i64 = -1;

        for received_result in rx {
            let buffers = match received_result {
                Ok(b) => b,
                Err(e) => {
                    eprintln!("{}", e);
                    break;
                }
            };
            let read_buffer1 = &buffers[0];
            let read_buffer2 = &buffers[1];
            let frames_in_block = read_buffer1.len() / bytes_per_frame;
            let block_start_frame_idx = processed_frames;

            // Calculate these values for the current block
            let block_current_geometric_delay_s = geometric_delay_s_initial;
            let block_current_geometric_rate_s = geometric_rate_s_initial;
            let block_current_geometric_accel_s = geometric_accel_s_initial;

            let (block_cross_spec_sum, block_auto_spec1_sum, block_auto_spec2_sum): (
                Vec<Complex<f64>>,
                Vec<f64>,
                Vec<f64>,
            ) = read_buffer1
                .par_chunks(bytes_per_batch_for_corr)
                .zip(read_buffer2.par_chunks(bytes_per_batch_for_corr))
                .enumerate()
                .fold(
                    || {
                        (
                            vec![Complex::new(0.0, 0.0); half_spec_len_for_corr],
                            vec![0.0f64; half_spec_len_for_corr],
                            vec![0.0f64; half_spec_len_for_corr],
                        )
                    },
                    |(mut acc_cross, mut acc_auto1, mut acc_auto2),
                     (batch_idx, (raw1_batch, raw2_batch))| {
                        let args_for_corr = Arc::clone(&args_for_corr);
                        let levels_for_corr = Arc::clone(&levels_for_corr);
                        let helper_for_corr = Arc::clone(&helper_for_corr);
                        let decode_plan_for_corr = Arc::clone(&decode_plan_for_corr);

                        let base_frame_idx = batch_idx * frames_per_batch_for_corr;

                        let mut frame1_f64 = vec![0.0f64; args_for_corr.fft];
                        let mut frame2_f64 = vec![0.0f64; args_for_corr.fft];
                        let mut s1_complex_half =
                            vec![Complex::new(0.0, 0.0); half_spec_len_for_corr];
                        let mut s2_complex_half =
                            vec![Complex::new(0.0, 0.0); half_spec_len_for_corr];
                        let mut bit_buffer = Vec::with_capacity(64);

                        for (local_idx, (raw1, raw2)) in raw1_batch
                            .chunks(bytes_per_frame)
                            .zip(raw2_batch.chunks(bytes_per_frame))
                            .enumerate()
                        {
                            decode_block_into_with_plan(
                                raw1,
                                levels_for_corr.as_ref(),
                                args_for_corr.fft,
                                decode_plan_for_corr.as_ref(),
                                &mut bit_buffer,
                                &mut frame1_f64,
                                lsb_to_usb,
                            )
                            .expect("Failed to decode block 1");
                            decode_block_into_with_plan(
                                raw2,
                                levels_for_corr.as_ref(),
                                args_for_corr.fft,
                                decode_plan_for_corr.as_ref(),
                                &mut bit_buffer,
                                &mut frame2_f64,
                                lsb_to_usb,
                            )
                            .expect("Failed to decode block 2");

                            // Keep raw DC component in .cor generation path.

                            // Perform R2C FFT
                            helper_for_corr
                                .forward_r2c_process(&mut frame1_f64, &mut s1_complex_half)
                                .expect("R2C FFT failed");
                            helper_for_corr
                                .forward_r2c_process(&mut frame2_f64, &mut s2_complex_half)
                                .expect("R2C FFT failed");

                            // Full per-antenna auto-spectrum for plotting (relative/baseband axis).
                            accumulate_power_add(&mut acc_auto1, &s1_complex_half);
                            accumulate_power_add(&mut acc_auto2, &s2_complex_half);

                            let frame_idx_global =
                                block_start_frame_idx + base_frame_idx + local_idx;
                            let frame_time_since_start =
                                (frame_idx_global as f64 + 0.5) * seconds_per_frame_for_corr;

                            let total_delay_diff = (block_current_geometric_delay_s
                                + coarse_delay_s_for_corr
                                + clock_delay_s_for_corr
                                + residual_delay_s_for_corr)
                                - delay_applied_by_seek_s_for_corr;
                            let total_rate_diff =
                                block_current_geometric_rate_s + fixed_rate_delay_sps_for_corr;
                            let total_accel_diff = block_current_geometric_accel_s;

                            let (ref_delay_s, ref_rate_sps, ref_accel_sps2) =
                                match delay_reference_for_corr {
                                    DelayReference::Ant2 => {
                                        (total_delay_diff, total_rate_diff, total_accel_diff)
                                    }
                                    DelayReference::Ant1 => {
                                        (-total_delay_diff, -total_rate_diff, -total_accel_diff)
                                    }
                                };
                            let fringe_rot = carrier_phase_rotation(
                                observing_frequency_hz_for_corr,
                                ref_delay_s,
                                ref_rate_sps,
                                ref_accel_sps2,
                                frame_time_since_start,
                            );

                            match delay_reference_for_corr {
                                DelayReference::Ant2 => apply_delay_and_rate_regular_bins(
                                    &mut s1_complex_half,
                                    args_for_corr.fft,
                                    fft_bin_step_hz_for_corr,
                                    total_delay_diff,
                                    total_rate_diff,
                                    total_accel_diff,
                                    frame_time_since_start,
                                ),
                                DelayReference::Ant1 => apply_delay_and_rate_regular_bins(
                                    &mut s2_complex_half,
                                    args_for_corr.fft,
                                    fft_bin_step_hz_for_corr,
                                    -total_delay_diff,
                                    -total_rate_diff,
                                    -total_accel_diff,
                                    frame_time_since_start,
                                ),
                            }

                            for k in 0..overlap_len_for_corr {
                                let i1 = ant1_overlap_start_for_corr + k;
                                let i2 = ant2_overlap_start_for_corr + k;
                                let s1_aligned = s1_complex_half[i1];
                                let s2_aligned = s2_complex_half[i2];
                                match delay_reference_for_corr {
                                    DelayReference::Ant2 => {
                                        acc_cross[i1] +=
                                            (s1_aligned * fringe_rot) * s2_aligned.conj();
                                    }
                                    DelayReference::Ant1 => {
                                        acc_cross[i1] +=
                                            (s2_aligned * fringe_rot) * s1_aligned.conj();
                                    }
                                }
                            }
                        }
                        (acc_cross, acc_auto1, acc_auto2)
                    },
                )
                .reduce(
                    || {
                        (
                            vec![Complex::new(0.0, 0.0); half_spec_len_for_corr],
                            vec![0.0f64; half_spec_len_for_corr],
                            vec![0.0f64; half_spec_len_for_corr],
                        )
                    },
                    |(mut acc1_cross, mut acc1_auto1, mut acc1_auto2),
                     (acc2_cross, acc2_auto1, acc2_auto2)| {
                        for i in 0..half_spec_len_for_corr {
                            acc1_cross[i] += acc2_cross[i];
                            acc1_auto1[i] += acc2_auto1[i];
                            acc1_auto2[i] += acc2_auto2[i];
                        }
                        (acc1_cross, acc1_auto1, acc1_auto2)
                    },
                );

            for i in 0..half_spec_len {
                accumulated_cross_spec[i] += block_cross_spec_sum[i] * inv_fft2_corr;
                accumulated_auto_spec1[i] += block_auto_spec1_sum[i] * inv_fft2_corr;
                accumulated_auto_spec2[i] += block_auto_spec2_sum[i] * inv_fft2_corr;
            }

            processed_frames += frames_in_block;

            print!("\rCorrelating ({}/{})", processed_frames, total_frames);
            std::io::stdout().flush()?;

            let processed_seconds = processed_frames as f64 * seconds_per_frame;
            let progress_second = processed_seconds.floor() as i64;
            if progress_second > last_progress_second {
                last_progress_second = progress_second;
            }
        }

        if producer_thread.join().is_err() {
            return Err("I/O producer thread panicked".into());
        }
        println!();
        println!(
            "[info] Final processed_frames after consumer loop: {}",
            processed_frames
        ); // Debug print

        if processed_frames > 0 {
            accumulated_auto_spec1[0] = 0.0;
            accumulated_auto_spec2[0] = 0.0;

            let power_norm_corr = (processed_frames as f64 * args.fft as f64).max(1.0);

            let invert_lag_axis = delay_reference == DelayReference::Ant1;
            let cross_result = finalize_cross_spectrum(
                &mut accumulated_cross_spec,
                helper.as_ref(),
                args.fft,
                half_spec_len,
                &lags,
                &freqs_mhz,
                invert_lag_axis,
            )?;

            let _auto_mag1 = finalize_auto_spectrum(&mut accumulated_auto_spec1, power_norm_corr)?;
            let _auto_mag2 = finalize_auto_spectrum(&mut accumulated_auto_spec2, power_norm_corr)?;

            if let Some(delay_from_phase) = cross_result.delay_seconds_from_phase {
                let delay_for_reference = if delay_reference == DelayReference::Ant2 {
                    delay_from_phase
                } else {
                    -delay_from_phase
                };
                println!(
                    "[info] Residual delay from unwrapped phase slope: {:.6e} s ({:.6} samples)",
                    delay_for_reference,
                    delay_for_reference * sampling_rate_hz
                );
            } else {
                println!(
                    "[warn] Unable to derive residual delay from phase slope (insufficient span)."
                );
            }

            println!("--------------------------------------------------");
            println!("Correlation analysis complete.");
            println!("Inspect plots for delay peak and phase slope.");
            println!("--------------------------------------------------");
        } else {
            println!("No frames were processed for correlation.");
        }
    }

    let samples_skipped_ant1 = total_samples_sought_ant1.max(0) as usize;
    let samples_skipped_ant2 = total_samples_sought_ant2.max(0) as usize;

    if total_frames == 0 {
        println!("[warn] No overlapping frames available for phased array synthesis.");
    } else {
        println!("[info] Synthesising phased array output...");

        let mut phase_reader1 = BufReader::new(File::open(&ant1_path)?);
        let mut phase_reader2 = BufReader::new(File::open(&ant2_path)?);
        let mut writer = BufWriter::new(File::create(&output_path)?);
        advise_file_sequential(phase_reader1.get_ref());
        advise_file_sequential(phase_reader2.get_ref());
        advise_file_sequential(writer.get_ref());

        if samples_skipped_ant1 > 0 {
            let bits_to_skip = samples_skipped_ant1 as u64 * bit_depth as u64;
            if bits_to_skip % 8 != 0 {
                return Err("Attempted to seek a non byte-aligned offset for ant1".into());
            }
            let bytes_to_skip = bits_to_skip / 8;
            phase_reader1.seek(std::io::SeekFrom::Current(bytes_to_skip as i64))?;
        }
        if samples_skipped_ant2 > 0 {
            let bits_to_skip = samples_skipped_ant2 as u64 * bit_depth as u64;
            if bits_to_skip % 8 != 0 {
                return Err("Attempted to seek a non byte-aligned offset for ant2".into());
            }
            let bytes_to_skip = bits_to_skip / 8;
            phase_reader2.seek(std::io::SeekFrom::Current(bytes_to_skip as i64))?;
        }
        let phase_reader1_start_pos = phase_reader1.stream_position()?;
        let phase_reader2_start_pos = phase_reader2.stream_position()?;

        let weight_sum = weight1 + weight2;
        if weight_sum == 0.0 {
            return Err("Combined antenna weights sum to zero; cannot scale phased output".into());
        }

        let mut frames_emitted = 0usize;
        // フェーズ合成後のスペクトル統計はフレーム全体で蓄積する
        let half_spec_len = args.fft / 2 + 1;
        let cor_bins = args.fft / 2;
        if cor_bins == 0 {
            return Err("FFT must be at least 2 for phased-array correlation bins".into());
        }
        let overlap_len = overlap_bins.min(cor_bins);
        let ant1_overlap_start = band_alignment.ant1_start.min(cor_bins);
        let ant2_overlap_start = band_alignment.ant2_start.min(cor_bins);
        let fft_bin_step_hz_for_synth = sampling_rate_hz / args.fft as f64;
        let mut total_combined_auto_spec = vec![0.0f64; half_spec_len];
        let mut total_ant1_auto_spec = vec![0.0f64; half_spec_len];
        let mut total_ant2_auto_spec = vec![0.0f64; half_spec_len];
        let geom_delay_fixed = geometric_delay_s_initial;
        let geom_rate_fixed = geometric_rate_s_initial;
        let geom_accel_fixed = geometric_accel_s_initial;

        let helper_arc = Arc::clone(&helper);
        const RAM_SECONDS_BUFFER: usize = 4;
        const TARGET_L3_TILE_BYTES_PAIR: usize = 4 * 1024 * 1024;
        let frames_per_batch_synth = (TARGET_L3_TILE_BYTES_PAIR / (bytes_per_frame * 2)).max(1);
        println!(
            "[info] L3 tile target per rayon task: {} bytes (~{} frames)",
            TARGET_L3_TILE_BYTES_PAIR, frames_per_batch_synth
        );

        let total_duration_s = total_frames as f64 * seconds_per_frame;
        let total_seconds_target =
            ((total_frames + frames_per_second - 1) / frames_per_second).max(1);
        let bytes_per_second_pair = frames_per_second
            .saturating_mul(bytes_per_frame)
            .saturating_mul(2);
        let target_ram_bytes = bytes_per_second_pair.saturating_mul(RAM_SECONDS_BUFFER);
        println!(
            "[info] RAM ring target: {} s, {} bytes (~{:.3} GiB)",
            RAM_SECONDS_BUFFER,
            target_ram_bytes,
            target_ram_bytes as f64 / (1024.0 * 1024.0 * 1024.0)
        );

        let mut second_frame_counts = Vec::with_capacity(total_seconds_target);
        let mut frames_left_for_seconds = total_frames;
        while frames_left_for_seconds > 0 {
            let frames_this_second = frames_left_for_seconds.min(frames_per_second);
            second_frame_counts.push(frames_this_second);
            frames_left_for_seconds -= frames_this_second;
        }

        let cor_dir = &output_dir;
        let cor_path_ant11 = cor_dir.join(format!(
            "YAMAGU32_YAMAGU32_{}_phasedarray.cor",
            cor_epoch_tag
        ));
        let cor_path_ant12 = cor_dir.join(format!(
            "YAMAGU32_YAMAGU34_{}_phasedarray.cor",
            cor_epoch_tag
        ));
        let cor_path_ant22 = cor_dir.join(format!(
            "YAMAGU34_YAMAGU34_{}_phasedarray.cor",
            cor_epoch_tag
        ));
        let cor_path_phased = cor_dir.join(format!(
            "YAMAGU66_YAMAGU66_{}_phasedarray.cor",
            cor_epoch_tag
        ));

        let sampling_speed_hz_i32 = i32::try_from(sampling_rate_hz.round() as i64)
            .map_err(|_| "sampling rate out of i32 range for .cor header")?;
        let fft_point_i32 =
            i32::try_from(args.fft).map_err(|_| "fft point out of i32 range for .cor header")?;
        let number_of_sector_hint_i32 = i32::try_from(total_seconds_target)
            .map_err(|_| "number_of_sector out of i32 range for .cor header")?;
        let source_ra_rad = geom_debug_info.map(|v| v.ra_rad).unwrap_or(0.0);
        let source_dec_rad = geom_debug_info.map(|v| v.dec_rad).unwrap_or(0.0);
        let source_name = match ifile_data.as_ref().and_then(|d| d.source.as_deref()) {
            Some(name) if !name.trim().is_empty() => name.trim().to_string(),
            _ => {
                println!(
                    "[warn] ifile source key is missing; using fallback source name PHASEDARRAY."
                );
                "PHASEDARRAY".to_string()
            }
        };

        let cor_header_cfg = CorHeaderConfig {
            sampling_speed_hz: sampling_speed_hz_i32,
            observing_frequency_hz: cor_observing_frequency_hz,
            fft_point: fft_point_i32,
            number_of_sector_hint: number_of_sector_hint_i32,
            clock_reference_unix_sec: cor_epoch_unix_s,
            source_name,
            source_ra_rad,
            source_dec_rad,
        };
        let ant1_station = CorStation {
            name: "YAMAGU32",
            code: b'K',
            ecef_m: geom::YAMAGU32_ECEF,
        };
        let ant2_station = CorStation {
            name: "YAMAGU34",
            code: b'L',
            ecef_m: geom::YAMAGU34_ECEF,
        };
        let phased_station = CorStation {
            name: "YAMAGU66",
            code: b'M',
            ecef_m: [
                0.5 * (geom::YAMAGU32_ECEF[0] + geom::YAMAGU34_ECEF[0]),
                0.5 * (geom::YAMAGU32_ECEF[1] + geom::YAMAGU34_ECEF[1]),
                0.5 * (geom::YAMAGU32_ECEF[2] + geom::YAMAGU34_ECEF[2]),
            ],
        };
        let mut cor_writer_ant11 =
            CorWriter::create(&cor_path_ant11, &cor_header_cfg, ant1_station, ant1_station)?;
        let mut cor_writer_ant12 =
            CorWriter::create(&cor_path_ant12, &cor_header_cfg, ant1_station, ant2_station)?;
        let mut cor_writer_ant22 =
            CorWriter::create(&cor_path_ant22, &cor_header_cfg, ant2_station, ant2_station)?;
        let mut cor_writer_phased = CorWriter::create(
            &cor_path_phased,
            &cor_header_cfg,
            phased_station,
            phased_station,
        )?;

        if let Some(info) = geom_debug_info {
            let total_seconds = total_duration_s.ceil() as usize;
            if total_seconds > 0 {
                let total_seconds = total_seconds.max(1);
                println!("[info] Geometric delay per second (ant2 - ant1):");
                for sec in 0..total_seconds {
                    let mjd = info.initial_mjd + sec as f64 / 86400.0;
                    let (_, _, geom_delay, _, _) = geom::calculate_geometric_delay_and_derivatives(
                        geom::YAMAGU32_ECEF,
                        geom::YAMAGU34_ECEF,
                        info.ra_rad,
                        info.dec_rad,
                        mjd,
                    );
                    let base_samples = geom_delay * sampling_rate_hz;
                    let adjusted_delay = geom_delay + args.gico3_correct;
                    let adjusted_samples = adjusted_delay * sampling_rate_hz;
                    println!(
                        "  t+{:>3} s: base {:.9e} s ({:.3} samples) | +gico3 {:.9e} s ({:.3} samples)",
                        sec,
                        geom_delay,
                        base_samples,
                        adjusted_delay,
                        adjusted_samples
                    );
                }
            }
        }

        struct SynthAccum {
            combined_auto_sum: Vec<f64>,
            ant1_auto_sum: Vec<f64>,
            ant2_auto_sum: Vec<f64>,
            corr_cross_sum: Vec<Complex<f64>>,
            corr_auto1_sum: Vec<f64>,
            corr_auto2_sum: Vec<f64>,
            frame_count: usize,
        }
        impl SynthAccum {
            fn new(half_spec_len: usize) -> Self {
                Self {
                    combined_auto_sum: vec![0.0; half_spec_len],
                    ant1_auto_sum: vec![0.0; half_spec_len],
                    ant2_auto_sum: vec![0.0; half_spec_len],
                    corr_cross_sum: vec![Complex::new(0.0, 0.0); half_spec_len],
                    corr_auto1_sum: vec![0.0; half_spec_len],
                    corr_auto2_sum: vec![0.0; half_spec_len],
                    frame_count: 0,
                }
            }

            fn merge_from(&mut self, other: Self) {
                for (acc, v) in self
                    .combined_auto_sum
                    .iter_mut()
                    .zip(other.combined_auto_sum.iter())
                {
                    *acc += *v;
                }
                for (acc, v) in self.ant1_auto_sum.iter_mut().zip(other.ant1_auto_sum.iter()) {
                    *acc += *v;
                }
                for (acc, v) in self.ant2_auto_sum.iter_mut().zip(other.ant2_auto_sum.iter()) {
                    *acc += *v;
                }
                for (acc, v) in self.corr_cross_sum.iter_mut().zip(other.corr_cross_sum.iter()) {
                    *acc += *v;
                }
                for (acc, v) in self.corr_auto1_sum.iter_mut().zip(other.corr_auto1_sum.iter()) {
                    *acc += *v;
                }
                for (acc, v) in self.corr_auto2_sum.iter_mut().zip(other.corr_auto2_sum.iter()) {
                    *acc += *v;
                }
                self.frame_count += other.frame_count;
            }
        }

        #[derive(Default)]
        struct SynthWorkBuffers {
            frame1_f64: Vec<f64>,
            frame2_f64: Vec<f64>,
            spectrum1_half: Vec<Complex<f64>>,
            spectrum2_half: Vec<Complex<f64>>,
            shifted_half: Vec<Complex<f64>>,
            combined_buffer_half: Vec<Complex<f64>>,
            time_domain_output: Vec<f64>,
            bit_buffer: Vec<u8>,
            encoded_frame: Vec<u8>,
        }

        let mut chunk_queue: VecDeque<SecondRawChunk> = VecDeque::with_capacity(RAM_SECONDS_BUFFER);
        let decode_plan_for_synth = Arc::new(build_decode_plan(bit_depth, shuffle_in.as_ref())?);
        let mut next_second_to_load = 0usize;
        let mut next_start_frame_to_load = 0usize;
        while chunk_queue.len() < RAM_SECONDS_BUFFER
            && next_second_to_load < second_frame_counts.len()
        {
            let frames_requested = second_frame_counts[next_second_to_load];
            match read_second_chunk(
                &mut phase_reader1,
                &mut phase_reader2,
                frames_requested,
                bytes_per_frame,
                next_second_to_load,
                next_start_frame_to_load,
            )? {
                Some(chunk) => {
                    if chunk.frames < frames_requested {
                        println!(
                            "[warn] Second {} short read: requested {} frames, got {}",
                            next_second_to_load + 1,
                            frames_requested,
                            chunk.frames
                        );
                    }
                    next_start_frame_to_load += chunk.frames;
                    next_second_to_load += 1;
                    chunk_queue.push_back(chunk);
                }
                None => break,
            }
        }
        println!(
            "[info] RAM ring primed: {} / {} seconds",
            chunk_queue.len(),
            RAM_SECONDS_BUFFER
        );

        while let Some(chunk) = chunk_queue.pop_front() {
            let second_index = chunk.second_index;
            let second_start_frame = chunk.start_frame;
            let frames_this_second = chunk.frames;
            let chunk_offset_bytes = second_start_frame as u64 * bytes_per_frame as u64;
            let chunk_len_bytes = frames_this_second as u64 * bytes_per_frame as u64;
            drop_file_cache_range(
                phase_reader1.get_ref(),
                phase_reader1_start_pos.saturating_add(chunk_offset_bytes),
                chunk_len_bytes,
            );
            drop_file_cache_range(
                phase_reader2.get_ref(),
                phase_reader2_start_pos.saturating_add(chunk_offset_bytes),
                chunk_len_bytes,
            );
            let raw_block1 = chunk.raw1;
            let raw_block2 = chunk.raw2;
            let mut second_encoded = vec![0u8; frames_this_second * bytes_per_frame];
            let base_total_delay_diff =
                (geom_delay_fixed + coarse_delay_s + clock_delay_s + residual_delay_s)
                    - delay_applied_by_seek_s;
            let base_total_rate_diff = geom_rate_fixed + fixed_rate_delay_s_per_s;
            let base_total_accel_diff = geom_accel_fixed;

            let make_synth_work_buffers = || SynthWorkBuffers {
                frame1_f64: vec![0.0f64; args.fft],
                frame2_f64: vec![0.0f64; args.fft],
                spectrum1_half: vec![Complex::new(0.0, 0.0); args.fft / 2 + 1],
                spectrum2_half: vec![Complex::new(0.0, 0.0); args.fft / 2 + 1],
                shifted_half: vec![Complex::new(0.0, 0.0); args.fft / 2 + 1],
                combined_buffer_half: vec![Complex::new(0.0, 0.0); args.fft / 2 + 1],
                time_domain_output: vec![0.0f64; args.fft],
                bit_buffer: Vec::with_capacity(64),
                encoded_frame: Vec::with_capacity(bytes_per_frame),
            };

            let synth_acc = raw_block1
                .par_chunks(frames_per_batch_synth * bytes_per_frame)
                .zip(raw_block2.par_chunks(frames_per_batch_synth * bytes_per_frame))
                .zip(second_encoded.par_chunks_mut(frames_per_batch_synth * bytes_per_frame))
                .enumerate()
                .try_fold(
                    || {
                        (
                            make_synth_work_buffers(),
                            SynthAccum::new(half_spec_len),
                        )
                    },
                    |(mut buffers, mut acc),
                     (batch_idx, ((raw1_batch, raw2_batch), encoded_batch))| {
                        let decode_plan_for_synth = Arc::clone(&decode_plan_for_synth);
                        let base_frame_idx_in_second = batch_idx * frames_per_batch_synth;
                        let frame_time_since_start =
                            (second_start_frame + base_frame_idx_in_second) as f64 * seconds_per_frame
                                + 0.5 * seconds_per_frame;
                        let total_delay_diff = base_total_delay_diff;
                        let total_rate_diff = base_total_rate_diff;
                        let total_accel_diff = base_total_accel_diff;
                        let (fringe_delay_s, fringe_rate_sps, fringe_accel_sps2) = match delay_reference
                        {
                            DelayReference::Ant2 => {
                                (total_delay_diff, total_rate_diff, total_accel_diff)
                            }
                            DelayReference::Ant1 => {
                                (-total_delay_diff, -total_rate_diff, -total_accel_diff)
                            }
                        };
                        let (apply_delay_s, apply_rate_sps, apply_accel_sps2) = match delay_reference
                        {
                            DelayReference::Ant2 => {
                                (total_delay_diff, total_rate_diff, total_accel_diff)
                            }
                            DelayReference::Ant1 => {
                                (-total_delay_diff, -total_rate_diff, -total_accel_diff)
                            }
                        };
                        let mut fringe_rot_rec = PhaseRotationRecurrence::new(
                            phase_correction_frequency_hz,
                            fringe_delay_s,
                            fringe_rate_sps,
                            fringe_accel_sps2,
                            frame_time_since_start,
                            seconds_per_frame,
                        );
                        let mut bin_step_rec = PhaseRotationRecurrence::new(
                            fft_bin_step_hz_for_synth,
                            apply_delay_s,
                            apply_rate_sps,
                            apply_accel_sps2,
                            frame_time_since_start,
                            seconds_per_frame,
                        );
                        let mut frame_time_s = frame_time_since_start;
                        const PHASE_REBASE_INTERVAL: usize = 4096;

                        for (local_idx, (raw1, raw2)) in raw1_batch
                            .chunks(bytes_per_frame)
                            .zip(raw2_batch.chunks(bytes_per_frame))
                            .enumerate()
                        {
                            let frame_idx_in_second = base_frame_idx_in_second + local_idx;
                            if frame_idx_in_second >= frames_this_second {
                                break;
                            }
                            if local_idx > 0 && (local_idx & (PHASE_REBASE_INTERVAL - 1)) == 0 {
                                fringe_rot_rec = PhaseRotationRecurrence::new(
                                    phase_correction_frequency_hz,
                                    fringe_delay_s,
                                    fringe_rate_sps,
                                    fringe_accel_sps2,
                                    frame_time_s,
                                    seconds_per_frame,
                                );
                                bin_step_rec = PhaseRotationRecurrence::new(
                                    fft_bin_step_hz_for_synth,
                                    apply_delay_s,
                                    apply_rate_sps,
                                    apply_accel_sps2,
                                    frame_time_s,
                                    seconds_per_frame,
                                );
                            }

                            let frame1_f64_buf = &mut buffers.frame1_f64;
                            let frame2_f64_buf = &mut buffers.frame2_f64;
                            let spectrum1_half_buf = &mut buffers.spectrum1_half;
                            let spectrum2_half_buf = &mut buffers.spectrum2_half;
                            let shifted_half_buf = &mut buffers.shifted_half;
                            let combined_buffer_half = &mut buffers.combined_buffer_half;
                            let time_domain_output = &mut buffers.time_domain_output;
                            let bit_buffer = &mut buffers.bit_buffer;

                            decode_block_into_with_plan(
                                raw1,
                                &levels,
                                args.fft,
                                decode_plan_for_synth.as_ref(),
                                bit_buffer,
                                frame1_f64_buf.as_mut_slice(),
                                lsb_to_usb,
                            )?;
                            decode_block_into_with_plan(
                                raw2,
                                &levels,
                                args.fft,
                                decode_plan_for_synth.as_ref(),
                                bit_buffer,
                                frame2_f64_buf.as_mut_slice(),
                                lsb_to_usb,
                            )?;

                            helper_arc
                                .forward_r2c_process(frame1_f64_buf, spectrum1_half_buf)
                                .expect("R2C FFT failed");
                            helper_arc
                                .forward_r2c_process(frame2_f64_buf, spectrum2_half_buf)
                                .expect("R2C FFT failed");

                            let fringe_rot = fringe_rot_rec.current();
                            let bin_rot_step = bin_step_rec.current();
                            match delay_reference {
                                DelayReference::Ant2 => {
                                    rotate_regular_bins_with_step(spectrum1_half_buf, args.fft, bin_rot_step)
                                }
                                DelayReference::Ant1 => {
                                    rotate_regular_bins_with_step(spectrum2_half_buf, args.fft, bin_rot_step)
                                }
                            }

                            // Full per-antenna auto-spectrum for plotting (relative/baseband axis).
                            accumulate_power_add(&mut acc.ant1_auto_sum, spectrum1_half_buf);
                            accumulate_power_add(&mut acc.ant2_auto_sum, spectrum2_half_buf);

                            shifted_half_buf.fill(Complex::new(0.0, 0.0));
                            if align_ref_is_ant1 {
                                for k in 0..overlap_len {
                                    let dst = ant1_overlap_start + k;
                                    let src = ant2_overlap_start + k;
                                    if dst < cor_bins && src < cor_bins {
                                        shifted_half_buf[dst] = spectrum2_half_buf[src];
                                    }
                                }
                                if cor_bins < shifted_half_buf.len() && cor_bins < spectrum2_half_buf.len()
                                {
                                    shifted_half_buf[cor_bins] = spectrum2_half_buf[cor_bins];
                                }
                                for k in 0..cor_bins {
                                    let s1 = spectrum1_half_buf[k];
                                    let s2 = shifted_half_buf[k];
                                    acc.corr_auto1_sum[k] += s1.norm_sqr();
                                    acc.corr_auto2_sum[k] += s2.norm_sqr();
                                    acc.corr_cross_sum[k] += match delay_reference {
                                        DelayReference::Ant2 => (s1 * fringe_rot) * s2.conj(),
                                        DelayReference::Ant1 => (s2 * fringe_rot) * s1.conj(),
                                    };
                                }
                                for k in 0..half_spec_len {
                                    let mut s1_for_phase = spectrum1_half_buf[k];
                                    let mut s2_for_phase = shifted_half_buf[k];
                                    match delay_reference {
                                        DelayReference::Ant2 => {
                                            s1_for_phase *= fringe_rot;
                                        }
                                        DelayReference::Ant1 => {
                                            s2_for_phase *= fringe_rot;
                                        }
                                    }
                                    combined_buffer_half[k] =
                                        s1_for_phase * weight1 + s2_for_phase * weight2;
                                }
                            } else {
                                for k in 0..overlap_len {
                                    let dst = ant2_overlap_start + k;
                                    let src = ant1_overlap_start + k;
                                    if dst < cor_bins && src < cor_bins {
                                        shifted_half_buf[dst] = spectrum1_half_buf[src];
                                    }
                                }
                                if cor_bins < shifted_half_buf.len() && cor_bins < spectrum1_half_buf.len()
                                {
                                    shifted_half_buf[cor_bins] = spectrum1_half_buf[cor_bins];
                                }
                                for k in 0..cor_bins {
                                    let s1 = shifted_half_buf[k];
                                    let s2 = spectrum2_half_buf[k];
                                    acc.corr_auto1_sum[k] += s1.norm_sqr();
                                    acc.corr_auto2_sum[k] += s2.norm_sqr();
                                    acc.corr_cross_sum[k] += match delay_reference {
                                        DelayReference::Ant2 => (s1 * fringe_rot) * s2.conj(),
                                        DelayReference::Ant1 => (s2 * fringe_rot) * s1.conj(),
                                    };
                                }
                                for k in 0..half_spec_len {
                                    let mut s1_for_phase = shifted_half_buf[k];
                                    let mut s2_for_phase = spectrum2_half_buf[k];
                                    match delay_reference {
                                        DelayReference::Ant2 => {
                                            s1_for_phase *= fringe_rot;
                                        }
                                        DelayReference::Ant1 => {
                                            s2_for_phase *= fringe_rot;
                                        }
                                    }
                                    combined_buffer_half[k] =
                                        s1_for_phase * weight1 + s2_for_phase * weight2;
                                }
                            }
                            accumulate_power_add(&mut acc.combined_auto_sum, combined_buffer_half);
                            if !combined_buffer_half.is_empty() {
                                combined_buffer_half[0].im = 0.0;
                                if args.fft % 2 == 0 {
                                    let nyquist_idx = args.fft / 2;
                                    if nyquist_idx < combined_buffer_half.len() {
                                        combined_buffer_half[nyquist_idx].im = 0.0;
                                    }
                                }
                            }

                            helper_arc
                                .inverse_c2r_process(combined_buffer_half, time_domain_output)?;
                            let encoded_frame = &mut buffers.encoded_frame;
                            quantise_frame(
                                time_domain_output,
                                bit_depth,
                                &levels,
                                shuffle_in.as_ref(),
                                encoded_frame,
                            )?;
                            let encoded_offset = local_idx * bytes_per_frame;
                            if encoded_offset + bytes_per_frame <= encoded_batch.len() {
                                encoded_batch[encoded_offset..encoded_offset + bytes_per_frame]
                                    .copy_from_slice(encoded_frame.as_slice());
                            }
                            acc.frame_count += 1;
                            frame_time_s += seconds_per_frame;
                            fringe_rot_rec.advance();
                            bin_step_rec.advance();
                        }
                        Ok((buffers, acc))
                    },
                )
                .map(|state: Result<(SynthWorkBuffers, SynthAccum), DynError>| {
                    state.map(|(_, acc)| acc)
                })
                .try_reduce(
                    || SynthAccum::new(half_spec_len),
                    |mut acc_a, acc_b| {
                        acc_a.merge_from(acc_b);
                        Ok(acc_a)
                    },
                )?;

            let SynthAccum {
                combined_auto_sum,
                ant1_auto_sum,
                ant2_auto_sum,
                corr_cross_sum: second_corr_cross_sum,
                corr_auto1_sum: second_corr_auto1_sum,
                corr_auto2_sum: second_corr_auto2_sum,
                frame_count: second_corr_frame_count,
            } = synth_acc;

            for (acc_bin, value) in total_combined_auto_spec
                .iter_mut()
                .zip(combined_auto_sum.iter())
            {
                *acc_bin += *value;
            }
            for (acc_bin, value) in total_ant1_auto_spec.iter_mut().zip(ant1_auto_sum.iter()) {
                *acc_bin += *value;
            }
            for (acc_bin, value) in total_ant2_auto_spec.iter_mut().zip(ant2_auto_sum.iter()) {
                *acc_bin += *value;
            }
            writer.write_all(&second_encoded)?;
            frames_emitted += second_corr_frame_count;

            if second_corr_frame_count > 0 {
                let fft_norm = (args.fft as f64) * (args.fft as f64);
                let inv = 1.0 / (second_corr_frame_count as f64 * fft_norm);
                let second_cross_spec: Vec<Complex<f32>> = second_corr_cross_sum
                    .iter()
                    .take(cor_bins)
                    .map(|v| Complex::new((v.re * inv) as f32, (v.im * inv) as f32))
                    .collect();
                let second_auto1_spec: Vec<Complex<f32>> = second_corr_auto1_sum
                    .iter()
                    .take(cor_bins)
                    .map(|v| Complex::new((v * inv) as f32, 0.0))
                    .collect();
                let second_auto2_spec: Vec<Complex<f32>> = second_corr_auto2_sum
                    .iter()
                    .take(cor_bins)
                    .map(|v| Complex::new((v * inv) as f32, 0.0))
                    .collect();
                let second_phased_auto_spec: Vec<Complex<f32>> = combined_auto_sum
                    .iter()
                    .take(cor_bins)
                    .map(|v| Complex::new((v * inv) as f32, 0.0))
                    .collect();
                let sector_timestamp = cor_epoch_unix_s + second_index as i64;
                let effective_integ_time_s =
                    (second_corr_frame_count as f64 * seconds_per_frame) as f32;
                let sector_sec_u32 = u32::try_from(sector_timestamp)
                    .map_err(|_| "sector timestamp out of u32 range for .cor model header")?;
                let mut ant1_model = CorClockModel {
                    sec: sector_sec_u32,
                    nsec: 0,
                    ..CorClockModel::default()
                };
                let mut ant2_model = CorClockModel {
                    sec: sector_sec_u32,
                    nsec: 0,
                    ..CorClockModel::default()
                };
                if let Some(info) = geom_debug_info {
                    let sector_start_time_s = second_start_frame as f64 * seconds_per_frame;
                    let sector_mjd = info.initial_mjd + sector_start_time_s / 86400.0;
                    let (delay1, rate1, acel1) = geom::calculate_antenna_delay_and_derivatives(
                        geom::YAMAGU32_ECEF,
                        info.ra_rad,
                        info.dec_rad,
                        sector_mjd,
                    );
                    let (delay2, rate2, acel2) = geom::calculate_antenna_delay_and_derivatives(
                        geom::YAMAGU34_ECEF,
                        info.ra_rad,
                        info.dec_rad,
                        sector_mjd,
                    );
                    ant1_model.delay = delay1;
                    ant1_model.rate = rate1;
                    ant1_model.acel = acel1;
                    ant2_model.delay = delay2;
                    ant2_model.rate = rate2;
                    ant2_model.acel = acel2;
                }
                let sector_model_11 = CorSectorModel {
                    station1: ant1_model,
                    station2: ant1_model,
                    amp: [0.0, 0.0],
                    phs: [0, 0],
                };
                let sector_model_12 = CorSectorModel {
                    station1: ant1_model,
                    station2: ant2_model,
                    amp: [0.0, 0.0],
                    phs: [0, 0],
                };
                let sector_model_22 = CorSectorModel {
                    station1: ant2_model,
                    station2: ant2_model,
                    amp: [0.0, 0.0],
                    phs: [0, 0],
                };
                let phased_model = CorClockModel {
                    sec: sector_sec_u32,
                    nsec: 0,
                    delay: 0.5 * (ant1_model.delay + ant2_model.delay),
                    rate: 0.5 * (ant1_model.rate + ant2_model.rate),
                    acel: 0.5 * (ant1_model.acel + ant2_model.acel),
                    jerk: 0.0,
                    snap: 0.0,
                };
                let sector_model_66 = CorSectorModel {
                    station1: phased_model,
                    station2: phased_model,
                    amp: [0.0, 0.0],
                    phs: [0, 0],
                };
                cor_writer_ant11.write_sector_with_model(
                    sector_timestamp,
                    effective_integ_time_s,
                    &second_auto1_spec,
                    Some(sector_model_11),
                )?;
                cor_writer_ant12.write_sector_with_model(
                    sector_timestamp,
                    effective_integ_time_s,
                    &second_cross_spec,
                    Some(sector_model_12),
                )?;
                cor_writer_ant22.write_sector_with_model(
                    sector_timestamp,
                    effective_integ_time_s,
                    &second_auto2_spec,
                    Some(sector_model_22),
                )?;
                cor_writer_phased.write_sector_with_model(
                    sector_timestamp,
                    effective_integ_time_s,
                    &second_phased_auto_spec,
                    Some(sector_model_66),
                )?;
            }

            let processed_seconds =
                (frames_emitted as f64 * seconds_per_frame).min(total_duration_s);
            let percentage = if total_duration_s > 0.0 {
                (processed_seconds / total_duration_s) * 100.0
            } else {
                100.0
            };
            println!(
                "[info] Synthesised second {}/{} ({:.2}%: {:.6} / {:.6} s)",
                second_index + 1,
                total_seconds_target,
                percentage,
                processed_seconds,
                total_duration_s
            );
            while chunk_queue.len() < RAM_SECONDS_BUFFER
                && next_second_to_load < second_frame_counts.len()
            {
                let frames_requested = second_frame_counts[next_second_to_load];
                match read_second_chunk(
                    &mut phase_reader1,
                    &mut phase_reader2,
                    frames_requested,
                    bytes_per_frame,
                    next_second_to_load,
                    next_start_frame_to_load,
                )? {
                    Some(next_chunk) => {
                        if next_chunk.frames < frames_requested {
                            println!(
                                "[warn] Second {} short read: requested {} frames, got {}",
                                next_second_to_load + 1,
                                frames_requested,
                                next_chunk.frames
                            );
                        }
                        next_start_frame_to_load += next_chunk.frames;
                        next_second_to_load += 1;
                        chunk_queue.push_back(next_chunk);
                    }
                    None => break,
                }
            }
        }

        if frames_emitted < total_frames {
            let pad_frames = total_frames - frames_emitted;
            println!(
                "[info] Input exhausted early; padding {} frame(s) with zeros",
                pad_frames
            );
            let zero_time = vec![0.0f64; args.fft];
            let mut encoded_pad = Vec::with_capacity(bytes_per_frame);
            quantise_frame(
                &zero_time,
                bit_depth,
                levels.as_ref(),
                shuffle_in.as_ref(),
                &mut encoded_pad,
            )?;
            if encoded_pad.len() != bytes_per_frame {
                encoded_pad.resize(bytes_per_frame, 0);
            }
            for _ in 0..pad_frames {
                writer.write_all(&encoded_pad)?;
            }
            frames_emitted += pad_frames;
            println!(
                "[info] Zero padding complete; synthesised total {:.2}% ({:.2} / {:.2} s)",
                100.0, total_duration_s, total_duration_s
            );
        }

        let cor_written_ant11 = cor_writer_ant11.finalize()?;
        let cor_written_ant12 = cor_writer_ant12.finalize()?;
        let cor_written_ant22 = cor_writer_ant22.finalize()?;
        let cor_written_phased = cor_writer_phased.finalize()?;
        println!(
            "[info] Wrote phased-array .cor files: {}, {}, {}, {}",
            cor_written_ant11.display(),
            cor_written_ant12.display(),
            cor_written_ant22.display(),
            cor_written_phased.display()
        );

        writer.flush()?;
        drop_file_cache_all(writer.get_ref());
        println!(
            "[info] Phased array output complete: {} frames ({:.3} s) written to {}",
            frames_emitted,
            frames_emitted as f64 * seconds_per_frame,
            output_path.display()
        );

        if frames_emitted > 0 {
            let half_spec = args.fft / 2;
            let plot_bin_end_exclusive = if args.fft % 2 == 0 {
                // Exclude Nyquist bin from spectrum plots.
                half_spec
            } else {
                half_spec + 1
            };
            let mut freqs_plot_mhz = Vec::with_capacity(plot_bin_end_exclusive);
            let sample_rate_mhz = sampling_rate_hz / 1_000_000.0;
            let plot_axis_low_mhz = if align_ref_is_ant1 {
                ant1_low_mhz
            } else {
                ant2_low_mhz
            };

            let power_norm = (frames_emitted as f64 * args.fft as f64).max(1.0);

            for i in 0..plot_bin_end_exclusive {
                let freq_base_mhz = (i as f64 * sample_rate_mhz) / args.fft as f64;
                freqs_plot_mhz.push(plot_axis_low_mhz + freq_base_mhz);
            }

            let amplitude_plot: Vec<f64> = total_combined_auto_spec
                .iter()
                .take(plot_bin_end_exclusive)
                .map(|&power_sum| (power_sum / power_norm).sqrt())
                .collect();

            if !freqs_plot_mhz.is_empty() {
                let stem = output_path
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("phased_array");
                let amp_path = output_dir.join(format!("{}_phased_spectrum_amplitude.png", stem));
                let auto_corr_path = output_dir.join(format!("{}_phased_autocorrelation.png", stem));
                let auto_spec_path = output_dir.join(format!("{}_phased_auto_spectrum.png", stem));
                let amp_path_s = amp_path.to_string_lossy().into_owned();
                let auto_corr_path_s = auto_corr_path.to_string_lossy().into_owned();
                let auto_spec_path_s = auto_spec_path.to_string_lossy().into_owned();

                let max_amp = amplitude_plot.iter().copied().fold(0.0_f64, f64::max);
                let amp_upper = if max_amp > 0.0 { max_amp } else { 1.0 };
                plot_series_f64_x(
                    &freqs_plot_mhz,
                    &amplitude_plot,
                    "Phased Array Spectrum Amplitude",
                    &amp_path_s,
                    "Frequency (MHz)",
                    "Amplitude",
                    Some((0.0, amp_upper)),
                    "Phased Array Amplitude",
                )?;
                let auto_amp_combined = amplitude_plot.clone();

                let auto_amp_ant1_full: Vec<f64> = total_ant1_auto_spec[0..plot_bin_end_exclusive]
                    .iter()
                    .map(|&power_sum| (power_sum / power_norm).sqrt())
                    .collect();
                let auto_amp_ant2_full: Vec<f64> = total_ant2_auto_spec[0..plot_bin_end_exclusive]
                    .iter()
                    .map(|&power_sum| (power_sum / power_norm).sqrt())
                    .collect();
                let mut auto_amp_ant1 = vec![0.0f64; plot_bin_end_exclusive];
                let mut auto_amp_ant2 = vec![0.0f64; plot_bin_end_exclusive];
                if matches!(rotation_shift_target, RotationShiftTarget::Ant1) {
                    auto_amp_ant2.copy_from_slice(&auto_amp_ant2_full);
                    for k in 0..overlap_len {
                        let dst = ant2_overlap_start + k;
                        let src = ant1_overlap_start + k;
                        if dst < cor_bins && src < auto_amp_ant1_full.len() && dst < auto_amp_ant1.len() {
                            auto_amp_ant1[dst] = auto_amp_ant1_full[src];
                        }
                    }
                } else {
                    auto_amp_ant1.copy_from_slice(&auto_amp_ant1_full);
                    for k in 0..overlap_len {
                        let dst = ant1_overlap_start + k;
                        let src = ant2_overlap_start + k;
                        if dst < cor_bins && src < auto_amp_ant2_full.len() && dst < auto_amp_ant2.len() {
                            auto_amp_ant2[dst] = auto_amp_ant2_full[src];
                        }
                    }
                }

                let mut max_val = 0.0f64;
                for value in auto_amp_combined
                    .iter()
                    .chain(auto_amp_ant1.iter())
                    .chain(auto_amp_ant2.iter())
                {
                    if *value > max_val {
                        max_val = *value;
                    }
                }
                let y_upper = if max_val > 0.0 { max_val * 1.05 } else { 1.0 };

                let ant1_label = if matches!(rotation_shift_target, RotationShiftTarget::Ant1)
                    && band_alignment.shift_bins != 0
                {
                    "ant1(shifted)"
                } else {
                    "ant1"
                };
                let ant2_label = if matches!(rotation_shift_target, RotationShiftTarget::Ant2)
                    && band_alignment.shift_bins != 0
                {
                    "ant2(shifted)"
                } else {
                    "ant2"
                };
                let auto_spec_series = [
                    (auto_amp_combined.as_slice(), &BLUE, "phased"),
                    (auto_amp_ant1.as_slice(), &RED, ant1_label),
                    (auto_amp_ant2.as_slice(), &GREEN, ant2_label),
                ];

                plot_multi_series_f64_x(
                    &freqs_plot_mhz,
                    &auto_spec_series,
                    "Auto-Spectrum Amplitude",
                    &auto_spec_path_s,
                    "Frequency (MHz)",
                    "Amplitude",
                    Some((0.0, y_upper)),
                )?;

                let auto_corr_spectrum_half: Vec<Complex<f64>> = total_combined_auto_spec
                    .iter()
                    .map(|&power| Complex::new(power, 0.0))
                    .collect();

                let mut full_auto_corr_spec = vec![Complex::new(0.0, 0.0); args.fft];
                full_auto_corr_spec[..half_spec_len].copy_from_slice(&auto_corr_spectrum_half);
                if args.fft > 1 {
                    let mirror_limit = if args.fft % 2 == 0 {
                        half_spec_len - 1
                    } else {
                        half_spec_len
                    };
                    for mirrored in 1..mirror_limit {
                        full_auto_corr_spec[args.fft - mirrored] =
                            auto_corr_spectrum_half[mirrored].conj();
                    }
                }
                helper.as_ref().inverse_c2c(&mut full_auto_corr_spec)?;
                let auto_corr_mag: Vec<f64> =
                    full_auto_corr_spec.iter().map(|c| c.norm()).collect();
                let shift = args.fft / 2;
                let auto_corr_shifted: Vec<f64> = auto_corr_mag
                    .iter()
                    .cycle()
                    .skip(shift)
                    .take(args.fft)
                    .copied()
                    .collect();
                let half = args.fft / 2;
                let lags: Vec<i32> = (-(half as isize)..half as isize)
                    .map(|i| {
                        if delay_reference == DelayReference::Ant1 {
                            -i as i32
                        } else {
                            i as i32
                        }
                    })
                    .collect();

                let norm_corr_factor = (frames_emitted as f64 * args.fft as f64).max(1.0);
                let auto_corr_phased: Vec<f64> = auto_corr_shifted
                    .iter()
                    .map(|v| v / norm_corr_factor)
                    .collect();

                let auto_corr_ant1_half: Vec<Complex<f64>> = total_ant1_auto_spec
                    .iter()
                    .map(|&power| Complex::new(power, 0.0))
                    .collect();
                let mut full_auto_corr_ant1_spec = vec![Complex::new(0.0, 0.0); args.fft];
                full_auto_corr_ant1_spec[..half_spec_len].copy_from_slice(&auto_corr_ant1_half);
                if args.fft > 1 {
                    let mirror_limit = if args.fft % 2 == 0 {
                        half_spec_len - 1
                    } else {
                        half_spec_len
                    };
                    for mirrored in 1..mirror_limit {
                        full_auto_corr_ant1_spec[args.fft - mirrored] =
                            auto_corr_ant1_half[mirrored].conj();
                    }
                }
                helper.as_ref().inverse_c2c(&mut full_auto_corr_ant1_spec)?; // Use inverse_c2c for full spectrum
                let auto_corr_ant1_mag: Vec<f64> =
                    full_auto_corr_ant1_spec.iter().map(|c| c.norm()).collect();
                let auto_corr_ant1_shifted: Vec<f64> = auto_corr_ant1_mag
                    .iter()
                    .cycle()
                    .skip(shift)
                    .take(args.fft)
                    .copied()
                    .collect();

                let auto_corr_ant1_norm: Vec<f64> = auto_corr_ant1_shifted
                    .iter()
                    .map(|v| v / norm_corr_factor)
                    .collect();

                let auto_corr_ant2_half: Vec<Complex<f64>> = total_ant2_auto_spec
                    .iter()
                    .map(|&power| Complex::new(power, 0.0))
                    .collect();
                let mut full_auto_corr_ant2_spec = vec![Complex::new(0.0, 0.0); args.fft];
                full_auto_corr_ant2_spec[..half_spec_len].copy_from_slice(&auto_corr_ant2_half);
                if args.fft > 1 {
                    let mirror_limit = if args.fft % 2 == 0 {
                        half_spec_len - 1
                    } else {
                        half_spec_len
                    };
                    for mirrored in 1..mirror_limit {
                        full_auto_corr_ant2_spec[args.fft - mirrored] =
                            auto_corr_ant2_half[mirrored].conj();
                    }
                }
                helper.as_ref().inverse_c2c(&mut full_auto_corr_ant2_spec)?; // Use inverse_c2c for full spectrum
                let auto_corr_ant2_mag: Vec<f64> =
                    full_auto_corr_ant2_spec.iter().map(|c| c.norm()).collect();
                let auto_corr_ant2_shifted: Vec<f64> = auto_corr_ant2_mag
                    .iter()
                    .cycle()
                    .skip(shift)
                    .take(args.fft)
                    .copied()
                    .collect();

                let auto_corr_ant2_norm: Vec<f64> = auto_corr_ant2_shifted
                    .iter()
                    .map(|v| v / norm_corr_factor)
                    .collect();

                // Find and print ACF peak positions
                let (phased_peak_idx, &phased_peak_val) = auto_corr_phased
                    .iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    .unwrap_or((0, &0.0));
                println!(
                    "[info] Phased ACF Peak: lag = {}, value = {:.6}",
                    lags[phased_peak_idx], phased_peak_val
                );

                let (ant1_peak_idx, &ant1_peak_val) = auto_corr_ant1_norm
                    .iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    .unwrap_or((0, &0.0));
                println!(
                    "[info] Ant1 ACF Peak: lag = {}, value = {:.6}",
                    lags[ant1_peak_idx], ant1_peak_val
                );

                let (ant2_peak_idx, &ant2_peak_val) = auto_corr_ant2_norm
                    .iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    .unwrap_or((0, &0.0));
                println!(
                    "[info] Ant2 ACF Peak: lag = {}, value = {:.6}",
                    lags[ant2_peak_idx], ant2_peak_val
                );

                fn shift_series(mut data: Vec<f64>, offset: isize) -> Vec<f64> {
                    if data.is_empty() {
                        return data;
                    }
                    let len = data.len() as isize;
                    let offset_mod = ((offset % len) + len) % len;
                    if offset_mod > 0 {
                        data.rotate_right(offset_mod as usize);
                    } else if offset_mod < 0 {
                        data.rotate_left((-offset_mod) as usize);
                    }
                    data
                }

                let auto_corr_ant1_offset = shift_series(auto_corr_ant1_norm, -100);
                let auto_corr_phased_offset = auto_corr_phased.clone();
                let auto_corr_ant2_offset = shift_series(auto_corr_ant2_norm, 100);

                let auto_corr_series = [
                    (auto_corr_phased_offset.as_slice(), &BLUE, "phased"),
                    (
                        auto_corr_ant1_offset.as_slice(),
                        &RED,
                        "ant1 (lag shift -100)",
                    ),
                    (
                        auto_corr_ant2_offset.as_slice(),
                        &GREEN,
                        "ant2 (lag shift +100)",
                    ),
                ];
                let auto_corr_upper = auto_corr_series
                    .iter()
                    .flat_map(|(series, _, _)| series.iter().copied())
                    .fold(0.0_f64, |acc, v| if v > acc { v } else { acc });
                let auto_corr_upper = if auto_corr_upper > 0.0 {
                    auto_corr_upper * 1.05
                } else {
                    1.0
                };
                let lags_f64: Vec<f64> = lags.iter().map(|&v| v as f64).collect();
                plot_multi_series_f64_x(
                    &lags_f64,
                    &auto_corr_series,
                    "Auto-Correlation Magnitude",
                    &auto_corr_path_s,
                    "Lag (samples, fftshifted)",
                    "Power / (frames × FFT)",
                    Some((0.0, auto_corr_upper)),
                )?;

                println!(
                    "[info] Wrote phased spectrum plots: amplitude -> {}, auto-spec -> {}, autocorr -> {}",
                    amp_path.display(),
                    auto_spec_path.display(),
                    auto_corr_path.display()
                );
            } else {
                println!("[warn] Insufficient spectrum bins to plot after removing DC component.");
            }
        } else {
            println!("[warn] No frames available for phased spectrum analysis.");
        }
    }

    println!("Processing finished.");

    Ok(())
}
