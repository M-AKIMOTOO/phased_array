use crate::utils::DynError;
use clap::Parser;
use std::f64::consts::PI;
use std::path::PathBuf;

const K_BOLTZMANN: f64 = 1.380_649e-23; // J/K
const SI_TO_JY: f64 = 1.0e26; // multiply W/m^2/Hz to express in Jy
pub const DEFAULT_FFT: usize = 16384;
pub const DEFAULT_LEVELS_2BIT_CSV: &str = "-1.5,-0.5,0.5,1.5";
pub const DEFAULT_SHUFFLE_IN: [usize; 32] = [
    // External convention: MSB=31 ... LSB=0
    // This default is identity in that convention.
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
    after_help = "Examples:\n  phased_array --ifile G9.62.ifile --data raw --fft 1024 --length 5\n  phased_array --ant1 raw/YAMAGU32.raw --ant2 raw/YAMAGU34.raw --ifile G9.62.ifile --fft 1024 --length 5\n  phased_array --ant1 raw/YAMAGU32.raw --ant2 raw/YAMAGU34.raw --ra 18h06m14.659 --dec -20d31m31.57s --epoch 2025-09-29T08:38:00Z --fft 4096\n  phased_array --ant1 raw/YAMAGU32.raw --ant2 raw/YAMAGU34.raw --ra 271.561 --dec -20.5254 --epoch 2000 --fft 16384\n"
)]
pub struct Args {
    /// Path to antenna 1 quantised data (optional when --ifile + --data is used)
    #[arg(long, visible_alias = "a1")]
    pub ant1: Option<PathBuf>,

    /// Path to antenna 2 quantised data (optional when --ifile + --data is used)
    #[arg(long, visible_alias = "a2")]
    pub ant2: Option<PathBuf>,

    /// Input file for geometric delay parameters
    #[arg(long)]
    pub ifile: Option<PathBuf>,

    /// Directory containing raw files for auto input resolution
    #[arg(long)]
    pub data: Option<PathBuf>,

    /// Source RA (supports hhmmss/hms or degrees)
    #[arg(long)]
    pub ra: Option<String>,

    /// Source Dec (supports ddmmss/dms or degrees)
    #[arg(long, allow_hyphen_values = true)]
    pub dec: Option<String>,

    /// Epoch for geometry (default: 2000). Accepts year, MJD, or ISO datetime.
    #[arg(long = "epoch", visible_alias = "ecpoh")]
    pub epoch: Option<String>,

    /// FFT length (number of samples processed)
    #[arg(long, default_value_t = DEFAULT_FFT)]
    pub fft: usize,

    /// Bits per sample in the quantised stream
    #[arg(long, default_value_t = 2)]
    pub bit: usize,

    /// Residual delay to apply to antenna 2 (samples)
    #[arg(long, allow_hyphen_values = true, default_value_t = 0.0)]
    pub delay: f64,

    /// Coarse delay correction (seconds). Default keeps Yamaguchi fixed value.
    #[arg(long, allow_hyphen_values = true)]
    pub coarse: Option<f64>,

    /// Additional coarse delay correction (seconds) inspired by gico3
    #[arg(
        long = "gico3-correct",
        visible_alias = "gico3",
        allow_hyphen_values = true,
        default_value_t = 0.0
    )]
    pub gico3_correct: f64,

    /// Residual fringe rate relative to clock reference (Hz)
    #[arg(long, default_value_t = 0.0, allow_hyphen_values = true)]
    pub rate: f64,

    /// Sky/base observing frequency in MHz (used for .cor header and rate correction)
    #[arg(long = "sky-freq", visible_alias = "obsfreq", default_value_t = 0.0)]
    pub sky_freq: f64,

    /// Sampling frequency of the raw data in MHz (default 1024 MHz = 1024000000 Hz)
    #[arg(long = "sampling", visible_alias = "fs", default_value_t = 1024.0)]
    pub sampling: f64,

    /// Duration of data to process in seconds
    #[arg(long, visible_alias = "sec")]
    pub length: Option<f64>,

    /// Optional output file for phased-array samples
    #[arg(long)]
    pub output: Option<PathBuf>,

    /// Comma-separated list of quantisation levels (2^bit entries)
    #[arg(long, allow_hyphen_values = true, default_value = DEFAULT_LEVELS_2BIT_CSV)]
    pub level: Option<String>,

    /// Optional comma-separated bit shuffle for decoding (32 entries)
    #[arg(long = "shuffle", visible_alias = "shuf", allow_hyphen_values = true)]
    pub shuffle_in: Option<String>,

    /// Sampler sideband (usb or lsb). If lsb, LSB->USB conversion is applied.
    #[arg(
        long = "sideband",
        visible_alias = "analysis-sideband",
        default_value = "usb",
        value_parser = clap::builder::PossibleValuesParser::new(["usb", "lsb"])
    )]
    pub sideband: String,

    /// System temperature for antenna 1 in kelvin (used for weighting)
    #[arg(long, default_value_t = 1.0)]
    pub tsys1: f64,

    /// System temperature for antenna 2 in kelvin (used for weighting)
    #[arg(long, default_value_t = 1.0)]
    pub tsys2: f64,

    /// Antenna diameter for antenna 1 in metres (used with --eta1 when SEFD not provided)
    #[arg(long)]
    pub diameter1: Option<f64>,

    /// Antenna diameter for antenna 2 in metres (used with --eta2 when SEFD not provided)
    #[arg(long)]
    pub diameter2: Option<f64>,

    /// Aperture efficiency for antenna 1 (default 0.65, used with --diameter1)
    #[arg(long, default_value_t = 0.65)]
    pub eta1: f64,

    /// Aperture efficiency for antenna 2 (default 0.65, used with --diameter2)
    #[arg(long, default_value_t = 0.65)]
    pub eta2: f64,

    /// Gain (or effective TA scaling) for antenna 1
    #[arg(long, default_value_t = 1.0)]
    pub gain1: f64,

    /// Gain (or effective TA scaling) for antenna 2
    #[arg(long, default_value_t = 1.0)]
    pub gain2: f64,

    /// System equivalent flux density for antenna 1 (overrides --gain1 when provided)
    #[arg(long)]
    pub sefd1: Option<f64>,

    /// System equivalent flux density for antenna 2 (overrides --gain2 when provided)
    #[arg(long)]
    pub sefd2: Option<f64>,

    /// Number of parallel worker threads
    #[arg(long, default_value_t = 2)]
    pub cpu: usize,

    /// Enable detailed debug logging for initial frames
    #[arg(long)]
    pub debug: bool,

    /// Number of frames to include in debug output when --debug is set
    #[arg(long, default_value_t = 4)]
    pub debug_frames: usize,

    /// Enable detailed correlation debug logging for each frame
    #[arg(long, visible_alias = "dcorr")]
    pub debug_corr: bool,

    /// Enable fringe search analysis
    #[arg(long, visible_alias = "fringe")]
    pub fringe_search: bool,

    /// Legacy option (ignored in fixed Yamaguchi mode). Delay is always applied to ant2 (YAMAGU34).
    #[arg(long, visible_alias = "ref", hide = true, default_value = "ant2", value_parser = clap::builder::PossibleValuesParser::new(["ant1", "ant2"]))]
    pub delay_reference: String,
}

pub fn parse_levels(bit: usize, list: Option<String>) -> Result<Vec<f64>, DynError> {
    if bit == 0 {
        return Err("Bit depth must be at least 1".into());
    }
    let expected = 1usize << bit;
    let levels = if let Some(values) = list {
        values
            .split(',')
            .map(|v| v.trim().parse::<f64>())
            .collect::<Result<Vec<_>, _>>()?
    } else {
        if expected == 4 {
            // Canonical 2-bit order: 00, 01, 10, 11
            vec![-1.5, -0.5, 0.5, 1.5]
        } else {
            return Err(
                "Provide --level with 2^bit comma-separated values for bit depths other than 2"
                    .into(),
            );
        }
    };
    if levels.len() != expected {
        return Err(format!(
            "Expected {expected} quantisation levels for {bit} bits, received {}",
            levels.len()
        )
        .into());
    }
    Ok(levels)
}

pub fn parse_shuffle(list: Option<String>, default: &[usize; 32]) -> Result<Vec<usize>, DynError> {
    let values_external = if let Some(values) = list {
        values
            .split(',')
            .map(|v| v.trim().parse::<usize>())
            .collect::<Result<Vec<_>, _>>()?
    } else {
        default.to_vec()
    };
    if values_external.len() != 32 {
        return Err("Shuffle map must contain exactly 32 entries".into());
    }
    let mut sorted = values_external.clone();
    sorted.sort_unstable();
    for (expected, &found) in (0usize..32).zip(sorted.iter()) {
        if expected != found {
            return Err("Shuffle map must be a permutation of 0..31".into());
        }
    }
    // External order is gico-style: entries are listed from output bit 31 down to 0.
    // Internal form is indexed by output bit number in LSB order:
    //   internal[out_bit_lsb] = input_bit_lsb
    let mut values_internal = vec![0usize; 32];
    for (idx_msb_to_lsb, input_bit) in values_external.into_iter().enumerate() {
        let out_bit_lsb = 31 - idx_msb_to_lsb;
        values_internal[out_bit_lsb] = input_bit;
    }
    Ok(values_internal)
}

pub fn resolve_weight(
    tsys: f64,
    gain: f64,
    sefd: Option<f64>,
    diameter: Option<f64>,
    eta: f64,
    label: &str,
) -> Result<(f64, f64, Option<f64>, Option<f64>), DynError> {
    if tsys == 0.0 {
        return Err(format!("{label} Tsys must be non-zero").into());
    }
    if let Some(sefd_value) = sefd {
        if sefd_value <= 0.0 {
            return Err(format!("{label} SEFD must be positive").into());
        }
        let weight = 1.0 / sefd_value;
        let effective_gain = tsys * weight;
        Ok((weight, effective_gain, Some(sefd_value), None))
    } else {
        if let Some(diameter_value) = diameter {
            if diameter_value <= 0.0 {
                return Err(format!("{label} diameter must be positive").into());
            }
            if eta <= 0.0 {
                return Err(format!("{label} aperture efficiency must be positive").into());
            }
            let radius = diameter_value / 2.0;
            let geom_area = PI * radius * radius;
            let effective_area = eta * geom_area;
            if effective_area <= 0.0 {
                return Err(format!("{label} effective area must be positive").into());
            }
            let sefd_si = (2.0 * K_BOLTZMANN * tsys) / effective_area;
            let sefd_jy = sefd_si * SI_TO_JY;
            let weight = 1.0 / sefd_jy;
            let effective_gain = tsys * weight;
            Ok((weight, effective_gain, Some(sefd_jy), Some(effective_area)))
        } else {
            if gain <= 0.0 {
                return Err(format!("{label} gain must be positive").into());
            }
            Ok((gain / tsys, gain, None, None))
        }
    }
}
