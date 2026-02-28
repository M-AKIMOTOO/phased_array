use crate::utils::DynError;
use clap::Parser;
use std::f64::consts::PI;
use std::path::PathBuf;

const K_BOLTZMANN: f64 = 1.380_649e-23; // J/K
const SI_TO_JY: f64 = 1.0e26; // multiply W/m^2/Hz to express in Jy
pub const DEFAULT_FFT: usize = 16384;
pub const DEFAULT_SHUFFLE_IN: [usize; 32] = [
    31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8,
    7, 6, 5, 4, 3, 2, 1, 0,
];

#[derive(Parser, Debug, Clone)]
#[command(
    author,
    version,
    about = "Simple two-antenna phased-array aligner",
    long_about = None,
    arg_required_else_help = true,
    after_help = "Examples:\n  phased_array --ifile G9.62.ifile --data raw --fft 1024 --length 5\n  phased_array --ifile 3C345.ifile --sideband ant1:LSB ant2:LSB --obsfreq 21812\n"
)]
pub struct Args {
    #[arg(long)]
    pub ant1: Option<PathBuf>,

    #[arg(long)]
    pub ant2: Option<PathBuf>,

    #[arg(long)]
    pub ifile: Option<PathBuf>,

    #[arg(long)]
    pub data: Option<PathBuf>,

    #[arg(long)]
    pub ra: Option<String>,

    #[arg(long, allow_hyphen_values = true)]
    pub dec: Option<String>,

    #[arg(long)]
    pub epoch: Option<String>,

    #[arg(long, default_value_t = DEFAULT_FFT)]
    pub fft: usize,

    /// Bits per sample (e.g. "2" or "ant1:2 ant2:4")
    #[arg(long, num_args = 1..)]
    pub bit: Vec<String>,

    /// Residual delay to apply to antenna 2 (samples)
    #[arg(long, allow_hyphen_values = true, default_value_t = 0.0)]
    pub delay: f64,

    #[arg(long, allow_hyphen_values = true)]
    pub coarse: Option<f64>,

    #[arg(long, allow_hyphen_values = true, default_value_t = 0.0)]
    pub gico3_correct: f64,

    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    pub rate: f64,

    #[arg(long, default_value_t = 0.0)]
    pub obsfreq: f64,

    #[arg(long, default_value_t = 1024.0)]
    pub sampling: f64,

    #[arg(long)]
    pub length: Option<f64>,

    #[arg(long)]
    pub output: Option<PathBuf>,

    /// Quantisation levels (e.g. "ant1:-1.5,0.5,-0.5,1.5" or "--level -1.5 0.5 -0.5 1.5")
    #[arg(long, allow_negative_numbers = true, num_args = 1..)]
    pub level: Vec<String>,

    /// Bit shuffle map (e.g. "ant1:31,30...0" or 32 space-separated values)
    #[arg(long = "shuffle", num_args = 1..=32)]
    pub shuffle_in: Vec<String>,

    /// Force VSREC decode map (override level/shuffle)
    #[arg(long)]
    pub vsrec: bool,

    /// Sampler sideband (e.g. "ant1:LSB ant2:USB")
    #[arg(long, num_args = 1.., default_value = "LSB")]
    pub sideband: Vec<String>,

    /// System temperature in Kelvin (e.g. "ant1:50 ant2:60")
    #[arg(long, num_args = 1..)]
    pub tsys: Vec<String>,

    /// Antenna diameter in metres (e.g. "ant1:32 ant2:34")
    #[arg(long, num_args = 1..)]
    pub diameter: Vec<String>,

    /// Aperture efficiency (e.g. "ant1:0.6 ant2:0.65")
    #[arg(long, num_args = 1..)]
    pub eta: Vec<String>,

    /// Gain / Scaling (e.g. "ant1:1 ant2:1")
    #[arg(long, num_args = 1..)]
    pub gain: Vec<String>,

    /// SEFD in Jy (e.g. "ant1:10000 ant2:12000")
    #[arg(long, num_args = 1..)]
    pub sefd: Vec<String>,

    #[arg(long, default_value_t = 2)]
    pub cpu: usize,

    /// Enable detailed debug logging to file
    #[arg(long)]
    pub debug: bool,

    #[arg(long, default_value_t = 4)]
    pub debug_frames: usize,

    /// Enable fringe search (cross-correlation analysis and plots)
    #[arg(long)]
    pub fringe: bool,

    #[arg(long, hide = true, default_value = "ant2", value_parser = clap::builder::PossibleValuesParser::new(["ant1", "ant2"]))]
    pub delay_reference: String,
}

pub fn resolve_per_antenna_config<T: Clone>(
    args_list: &[String],
    default_val: T,
    parser: impl Fn(&str) -> Result<T, DynError>,
) -> Result<(T, T), DynError> {
    let mut ant1_val = default_val.clone();
    let mut ant2_val = default_val;
    for entry in args_list {
        if entry.contains("ant1:") || entry.contains("ant2:") {
            let parts: Vec<&str> = entry.split(|c: char| c == ',' || c == ' ').filter(|s| !s.is_empty()).collect();
            for part in parts {
                if let Some(val_str) = part.strip_prefix("ant1:") { ant1_val = parser(val_str.trim())?; }
                else if let Some(val_str) = part.strip_prefix("ant2:") { ant2_val = parser(val_str.trim())?; }
                else { let v = parser(part.trim())?; ant1_val = v.clone(); ant2_val = v; }
            }
        } else {
            let v = parser(entry.trim())?; ant1_val = v.clone(); ant2_val = v;
        }
    }
    Ok((ant1_val, ant2_val))
}

pub fn parse_levels(bit: usize, list: &str) -> Result<Vec<f64>, DynError> {
    if bit == 0 { return Err("Bit depth must be at least 1".into()); }
    let expected = 1usize << bit;
    let levels = list.split(',').map(|v| v.trim().parse::<f64>()).collect::<Result<Vec<_>, _>>()?;
    if levels.len() != expected { return Err(format!("Expected {expected} levels, received {}", levels.len()).into()); }
    Ok(levels)
}

pub fn parse_shuffle(list: &str) -> Result<Vec<usize>, DynError> {
    let values_external: Vec<usize> = list.split(',').map(|v| v.trim().parse::<usize>()).collect::<Result<Vec<_>, _>>()?;
    if values_external.len() != 32 { return Err("Shuffle map must contain exactly 32 entries".into()); }
    let mut sorted = values_external.clone(); sorted.sort_unstable();
    for (expected, &found) in (0usize..32).zip(sorted.iter()) { if expected != found { return Err("Shuffle map must be a permutation of 0..31".into()); } }
    let mut values_internal = vec![0usize; 32];
    for (idx_msb_to_lsb, input_bit) in values_external.into_iter().enumerate() { let out_bit_lsb = 31 - idx_msb_to_lsb; values_internal[out_bit_lsb] = input_bit; }
    Ok(values_internal)
}

pub fn resolve_weight(tsys: f64, gain: f64, sefd: Option<f64>, diameter: Option<f64>, eta: f64, label: &str) -> Result<(f64, f64, Option<f64>, Option<f64>), DynError> {
    if tsys == 0.0 { return Err(format!("{label} Tsys must be non-zero").into()); }
    if let Some(sefd_v) = sefd {
        if sefd_v <= 0.0 { return Err(format!("{label} SEFD must be positive").into()); }
        let w = 1.0 / sefd_v; Ok((w, tsys * w, Some(sefd_v), None))
    } else if let Some(dia_v) = diameter {
        if dia_v <= 0.0 || eta <= 0.0 { return Err(format!("{label} diameter/eta must be positive").into()); }
        let geom_area = PI * (dia_v / 2.0).powi(2); let eff_area = eta * geom_area;
        let sefd_si = (2.0 * K_BOLTZMANN * tsys) / eff_area; let sefd_jy = sefd_si * SI_TO_JY;
        let w = 1.0 / sefd_jy; Ok((w, tsys * w, Some(sefd_jy), Some(eff_area)))
    } else {
        if gain <= 0.0 { return Err(format!("{label} gain must be positive").into()); }
        Ok((gain / tsys, gain, None, None))
    }
}
