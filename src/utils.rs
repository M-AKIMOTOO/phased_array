use std::f64::consts::PI;
use std::sync::Arc;

use realfft::{ComplexToReal, RealFftPlanner, RealToComplex};
use rustfft::num_complex::Complex;
use std::error::Error;

pub type DynError = Box<dyn Error + Send + Sync>;

pub struct FftHelper {
    len: usize,
    pub forward_r2c: Arc<dyn RealToComplex<f32>>,
    pub inverse_c2r: Arc<dyn ComplexToReal<f32>>,
}

pub struct FftScratch {
    forward_r2c: Vec<Complex<f32>>,
    inverse_c2r: Vec<Complex<f32>>,
}

impl FftHelper {
    pub fn new(len: usize) -> Self {
        let mut planner_r2c = RealFftPlanner::new();
        let forward_r2c = planner_r2c.plan_fft_forward(len);
        let inverse_c2r = planner_r2c.plan_fft_inverse(len);
        Self { len, forward_r2c, inverse_c2r }
    }
    pub fn make_scratch(&self) -> FftScratch {
        FftScratch {
            forward_r2c: self.forward_r2c.make_scratch_vec(),
            inverse_c2r: self.inverse_c2r.make_scratch_vec(),
        }
    }
    pub fn forward_r2c_process_with_scratch(
        &self,
        input: &mut [f32],
        output: &mut [Complex<f32>],
        scratch: &mut FftScratch,
    ) -> Result<(), DynError> {
        if input.len() != self.len || output.len() != self.len / 2 + 1 { return Err("Length mismatch".into()); }
        self.forward_r2c.process_with_scratch(input, output, &mut scratch.forward_r2c)?;
        Ok(())
    }
    pub fn inverse_c2r_process_with_scratch(
        &self,
        spectrum: &mut [Complex<f32>],
        output: &mut [f32],
        scratch: &mut FftScratch,
    ) -> Result<(), DynError> {
        if spectrum.len() != self.len / 2 + 1 || output.len() != self.len { return Err("Length mismatch".into()); }
        self.inverse_c2r.process_with_scratch(spectrum, output, &mut scratch.inverse_c2r)?;
        let scale = 1.0_f32 / self.len as f32;
        for value in output.iter_mut() { *value *= scale; }
        Ok(())
    }
}

pub fn apply_delay_and_rate_regular_bins(
    spectrum: &mut [Complex<f32>],
    fft_len: usize,
    freq_step_hz: f64,
    delay_seconds: f64,
    delay_rate: f64,
    geometric_accel: f64,
    frame_time: f64,
    is_lsb: bool,
) {
    apply_delay_and_rate_regular_bins_range(
        spectrum,
        fft_len,
        freq_step_hz,
        delay_seconds,
        delay_rate,
        geometric_accel,
        frame_time,
        is_lsb,
        0,
        spectrum.len(),
    );
}

pub fn apply_delay_and_rate_regular_bins_range(
    spectrum: &mut [Complex<f32>],
    fft_len: usize,
    freq_step_hz: f64,
    delay_seconds: f64,
    delay_rate: f64,
    geometric_accel: f64,
    frame_time: f64,
    is_lsb: bool,
    start_bin: usize,
    end_bin: usize,
) {
    if spectrum.is_empty() { return; }
    let end = end_bin.min(spectrum.len());
    let start = start_bin.min(end);
    if start >= end { return; }
    let total_delay_at_time = delay_seconds + delay_rate * frame_time + 0.5 * geometric_accel * frame_time.powi(2);
    // Standard VLBI sign: exp(-j * 2*PI * f * tau)
    // If we have flipped LSB to USB in time-domain, we use USB sign (-).
    // If not flipped, LSB uses positive sign (+).
    let sign = if is_lsb { 1.0 } else { -1.0 };
    let phase_step = sign * 2.0 * PI * freq_step_hz * total_delay_at_time;
    let mut rot = Complex::from_polar(1.0_f32, (phase_step * start as f64) as f32);
    let bin_rot_step = Complex::from_polar(1.0_f32, phase_step as f32);
    for bin in spectrum[start..end].iter_mut() { *bin *= rot; rot *= bin_rot_step; }
    if fft_len % 2 == 0 {
        let nyquist_idx = fft_len / 2;
        if nyquist_idx >= start && nyquist_idx < end && nyquist_idx < spectrum.len() {
            spectrum[nyquist_idx].im = 0.0_f32;
        }
    }
}

pub fn apply_integer_sample_shift_zerofill(samples: &mut [f32], shift_samples: i64) {
    if samples.is_empty() || shift_samples == 0 {
        return;
    }
    let n = samples.len();
    let s_abs = shift_samples.unsigned_abs() as usize;
    if s_abs >= n {
        samples.fill(0.0_f32);
        return;
    }
    if shift_samples > 0 {
        // Positive delay: shift right (later in time), zero-fill head.
        samples.copy_within(..(n - s_abs), s_abs);
        samples[..s_abs].fill(0.0_f32);
    } else {
        // Negative delay: shift left (earlier in time), zero-fill tail.
        samples.copy_within(s_abs.., 0);
        samples[(n - s_abs)..].fill(0.0_f32);
    }
}

pub struct DecodePlan {
    bits: usize,
    shuffle_kind: ShuffleKind,
    input_shifts: [u32; 32],
    fast_kind: DecodeFastKind,
    levels_f32: Vec<f32>,
    packed2_lut16_usb: Option<Box<[[f32; 8]]>>,
    packed2_lut16_lsb_flip: Option<Box<[[f32; 8]]>>,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)] enum ShuffleKind { Identity, PairSwap, Generic }
#[derive(Clone, Copy, Debug, PartialEq, Eq)] enum DecodeFastKind { Generic, Packed1, Packed2, Packed4, Packed8 }

#[inline]
fn apply_shuffle_word(mut word: u32, plan: &DecodePlan) -> u32 {
    if plan.shuffle_kind == ShuffleKind::PairSwap {
        word = ((word & 0xAAAA_AAAA) >> 1) | ((word & 0x5555_5555) << 1);
    } else if plan.shuffle_kind == ShuffleKind::Generic {
        let mut shuffled = 0u32;
        for (out_bit, &in_shift) in plan.input_shifts.iter().enumerate() {
            shuffled |= ((word >> in_shift) & 1) << out_bit;
        }
        word = shuffled;
    }
    word
}

pub fn build_decode_plan(bits: usize, shuffle_in: &[usize], levels: &[f64]) -> Result<DecodePlan, DynError> {
    let mut input_shifts = [0u32; 32];
    for (idx, &mapped) in shuffle_in.iter().enumerate() { input_shifts[idx] = mapped as u32; }
    let is_identity = shuffle_in.iter().enumerate().all(|(i, &v)| i == v);
    let is_pair_swap = shuffle_in.iter().enumerate().all(|(i, &v)| (i ^ 1) == v);
    let kind = if is_identity { ShuffleKind::Identity } else if is_pair_swap { ShuffleKind::PairSwap } else { ShuffleKind::Generic };
    let fast_kind = match bits {
        1 => DecodeFastKind::Packed1,
        2 => DecodeFastKind::Packed2,
        4 => DecodeFastKind::Packed4,
        8 => DecodeFastKind::Packed8,
        _ => DecodeFastKind::Generic,
    };
    let levels_f32 = levels.iter().map(|&v| v as f32).collect::<Vec<_>>();
    let (packed2_lut16_usb, packed2_lut16_lsb_flip) = if bits == 2 && kind == ShuffleKind::Identity {
        let mut lut_usb = vec![[0.0_f32; 8]; 1 << 16].into_boxed_slice();
        let mut lut_lsb_flip = vec![[0.0_f32; 8]; 1 << 16].into_boxed_slice();
        for word in 0u32..(1u32 << 16) {
            let w = word as usize;
            for i in 0..8 {
                let code = (w >> (2 * i)) & 0x3;
                let val = levels_f32[code];
                lut_usb[w][i] = val;
                // LSB->USB normalization flips odd-indexed time samples.
                lut_lsb_flip[w][i] = if (i & 1) == 1 { -val } else { val };
            }
        }
        (Some(lut_usb), Some(lut_lsb_flip))
    } else {
        (None, None)
    };
    Ok(DecodePlan {
        bits,
        shuffle_kind: kind,
        input_shifts,
        fast_kind,
        levels_f32,
        packed2_lut16_usb,
        packed2_lut16_lsb_flip,
    })
}

pub fn decode_block_into_with_plan(
    raw: &[u8],
    samples: usize,
    plan: &DecodePlan,
    output: &mut [f32],
    lsb_to_usb: bool,
) -> Result<(), DynError> {
    if output.len() < samples {
        return Err(format!("output buffer too short: {} < {}", output.len(), samples).into());
    }
    let level_map = plan.levels_f32.as_slice();

    let mut out_idx = 0usize;
    let mut odd = false;

    match plan.fast_kind {
        DecodeFastKind::Packed1 => {
            for chunk in raw.chunks_exact(4) {
                let mut word = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
                word = apply_shuffle_word(word, plan);
                for _ in 0..32 {
                    if out_idx >= samples { break; }
                    let code = (word & 0x1) as usize;
                    let mut val = level_map[code];
                    if lsb_to_usb && odd { val = -val; }
                    output[out_idx] = val;
                    out_idx += 1;
                    odd = !odd;
                    word >>= 1;
                }
            }
        }
        DecodeFastKind::Packed2 => {
            if plan.shuffle_kind == ShuffleKind::Identity {
                let lut = if lsb_to_usb {
                    plan.packed2_lut16_lsb_flip
                        .as_ref()
                        .ok_or("internal error: missing packed2(16b) LSB lookup table")?
                } else {
                    plan.packed2_lut16_usb
                        .as_ref()
                        .ok_or("internal error: missing packed2(16b) lookup table")?
                };
                let mut pairs = raw.chunks_exact(2);
                for pair in &mut pairs {
                    if out_idx + 8 <= samples {
                        let word = u16::from_le_bytes([pair[0], pair[1]]) as usize;
                        output[out_idx..out_idx + 8].copy_from_slice(&lut[word]);
                        out_idx += 8;
                    } else {
                        // Rare tail path when sample count is not aligned by 8 samples.
                        let mut word = u16::from_le_bytes([pair[0], pair[1]]) as usize;
                        for _ in 0..8 {
                            if out_idx >= samples {
                                break;
                            }
                            let code = word & 0x3;
                            let mut val = level_map[code];
                            if lsb_to_usb && odd {
                                val = -val;
                            }
                            output[out_idx] = val;
                            out_idx += 1;
                            odd = !odd;
                            word >>= 2;
                        }
                    }
                }
                for &byte in pairs.remainder() {
                    let mut packed = byte as usize;
                    for _ in 0..4 {
                        if out_idx >= samples {
                            break;
                        }
                        let code = packed & 0x3;
                        let mut val = level_map[code];
                        if lsb_to_usb && odd {
                            val = -val;
                        }
                        output[out_idx] = val;
                        out_idx += 1;
                        odd = !odd;
                        packed >>= 2;
                    }
                }
            } else {
                for chunk in raw.chunks_exact(4) {
                    let mut word = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
                    word = apply_shuffle_word(word, plan);
                    for _ in 0..16 {
                        if out_idx >= samples { break; }
                        let code = (word & 0x3) as usize;
                        let mut val = level_map[code];
                        if lsb_to_usb && odd { val = -val; }
                        output[out_idx] = val;
                        out_idx += 1;
                        odd = !odd;
                        word >>= 2;
                    }
                }
            }
        }
        DecodeFastKind::Packed4 => {
            for chunk in raw.chunks_exact(4) {
                let mut word = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
                word = apply_shuffle_word(word, plan);
                for _ in 0..8 {
                    if out_idx >= samples { break; }
                    let code = (word & 0xF) as usize;
                    let mut val = level_map[code];
                    if lsb_to_usb && odd { val = -val; }
                    output[out_idx] = val;
                    out_idx += 1;
                    odd = !odd;
                    word >>= 4;
                }
            }
        }
        DecodeFastKind::Packed8 => {
            for chunk in raw.chunks_exact(4) {
                let mut word = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
                word = apply_shuffle_word(word, plan);
                for _ in 0..4 {
                    if out_idx >= samples { break; }
                    let code = (word & 0xFF) as usize;
                    let mut val = level_map[code];
                    if lsb_to_usb && odd { val = -val; }
                    output[out_idx] = val;
                    out_idx += 1;
                    odd = !odd;
                    word >>= 8;
                }
            }
        }
        DecodeFastKind::Generic => {
            let bits = plan.bits;
            let code_mask = (1u64 << bits) - 1;
            let mut acc = 0u64;
            let mut acc_bits = 0;
            for chunk in raw.chunks_exact(4) {
                let word = apply_shuffle_word(u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]), plan);
                acc |= (word as u64) << acc_bits;
                acc_bits += 32;
                while acc_bits >= bits && out_idx < samples {
                    let code = (acc & code_mask) as usize;
                    let mut val = level_map[code];
                    if lsb_to_usb && odd { val = -val; }
                    output[out_idx] = val;
                    out_idx += 1;
                    odd = !odd;
                    acc_bits -= bits;
                    acc >>= bits;
                }
            }
        }
    }

    if out_idx != samples {
        return Err(format!("decoded {} samples, expected {}", out_idx, samples).into());
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum QuantizeShuffleKind {
    Identity,
    PairSwap,
    Generic,
}

pub struct QuantizePlan {
    bits: usize,
    code_mask: u64,
    levels_f32: Vec<f32>,
    shuffle_kind: QuantizeShuffleKind,
    shuffle_out: [u8; 32],
    fast_kind: QuantizeFastKind,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum QuantizeFastKind {
    Generic,
    Bits2DefaultLevels,
}

pub fn build_quantize_plan(bits: usize, levels: &[f64], shuffle_out: &[usize]) -> Result<QuantizePlan, DynError> {
    if bits == 0 {
        return Err("quantize bits must be >= 1".into());
    }
    if levels.is_empty() {
        return Err("quantize levels must not be empty".into());
    }
    if shuffle_out.len() != 32 {
        return Err("shuffle_out must contain exactly 32 entries".into());
    }
    let mut so = [0u8; 32];
    for (idx, &v) in shuffle_out.iter().enumerate() {
        if v >= 32 {
            return Err("shuffle_out values must be in 0..31".into());
        }
        so[idx] = v as u8;
    }
    let is_identity = shuffle_out.iter().enumerate().all(|(i, &v)| i == v);
    let is_pair_swap = shuffle_out.iter().enumerate().all(|(i, &v)| (i ^ 1) == v);
    let shuffle_kind = if is_identity {
        QuantizeShuffleKind::Identity
    } else if is_pair_swap {
        QuantizeShuffleKind::PairSwap
    } else {
        QuantizeShuffleKind::Generic
    };
    let levels_f32 = levels.iter().map(|&v| v as f32).collect::<Vec<_>>();
    let fast_kind = if bits == 2
        && levels.len() == 4
        && (levels[0] + 1.5).abs() < 1e-9
        && (levels[1] - 0.5).abs() < 1e-9
        && (levels[2] + 0.5).abs() < 1e-9
        && (levels[3] - 1.5).abs() < 1e-9
    {
        QuantizeFastKind::Bits2DefaultLevels
    } else {
        QuantizeFastKind::Generic
    };
    Ok(QuantizePlan {
        bits,
        code_mask: (1u64 << bits) - 1,
        levels_f32,
        shuffle_kind,
        shuffle_out: so,
        fast_kind,
    })
}

#[inline]
fn shuffle_word_quantize(word: u32, plan: &QuantizePlan) -> u32 {
    match plan.shuffle_kind {
        QuantizeShuffleKind::Identity => word,
        QuantizeShuffleKind::PairSwap => ((word & 0xAAAA_AAAA) >> 1) | ((word & 0x5555_5555) << 1),
        QuantizeShuffleKind::Generic => {
            let mut v = 0u32;
            for (pos, &target) in plan.shuffle_out.iter().enumerate() {
                v |= ((word >> pos) & 1) << target;
            }
            v
        }
    }
}

#[inline]
fn nearest_level_code(val: f32, levels: &[f32], fast_kind: QuantizeFastKind) -> u32 {
    if fast_kind == QuantizeFastKind::Bits2DefaultLevels {
        // Default 2-bit levels: code->value = 0:-1.5, 1:0.5, 2:-0.5, 3:1.5
        return if val < -1.0 {
            0
        } else if val < 0.0 {
            2
        } else if val < 1.0 {
            1
        } else {
            3
        };
    }
    let mut best_idx = 0usize;
    let mut min_err = f32::MAX;
    for (idx, &lv) in levels.iter().enumerate() {
        let err = (val - lv).abs();
        if err < min_err {
            min_err = err;
            best_idx = idx;
        }
    }
    best_idx as u32
}

#[inline]
fn quantize_code_2bit_default(val: f32) -> u32 {
    if val < -1.0 {
        0
    } else if val < 0.0 {
        2
    } else if val < 1.0 {
        1
    } else {
        3
    }
}

pub fn quantise_frame_with_plan(
    samples: &[f32],
    plan: &QuantizePlan,
    output: &mut Vec<u8>,
) -> Result<(), DynError> {
    if plan.fast_kind == QuantizeFastKind::Bits2DefaultLevels {
        let words = (samples.len() + 15) / 16;
        let target_bytes = words * 4;
        output.clear();
        if output.capacity() < target_bytes {
            output.reserve(target_bytes - output.capacity());
        }
        output.resize(target_bytes, 0);
        for wi in 0..words {
            let start = wi * 16;
            let end = (start + 16).min(samples.len());
            let mut word = 0u32;
            for (j, &val) in samples[start..end].iter().enumerate() {
                let code = quantize_code_2bit_default(val);
                word |= code << (2 * j);
            }
            let shuffled = shuffle_word_quantize(word, plan);
            let off = wi * 4;
            output[off..off + 4].copy_from_slice(&shuffled.to_le_bytes());
        }
        return Ok(());
    }

    output.clear();
    let bits = plan.bits;
    let target_bytes = ((samples.len().saturating_mul(bits) + 31) / 32).saturating_mul(4);
    if output.capacity() < target_bytes {
        output.reserve(target_bytes - output.capacity());
    }
    output.resize(target_bytes, 0);

    let code_mask = plan.code_mask;
    let mut bit_acc = 0u64;
    let mut acc_bits = 0;
    let mut out_words = 0usize;
    for &val in samples.iter() {
        let code = nearest_level_code(val, &plan.levels_f32, plan.fast_kind) as u64;
        bit_acc |= (code & code_mask) << acc_bits;
        acc_bits += bits;
        while acc_bits >= 32 {
            let word = (bit_acc & 0xFFFF_FFFF) as u32;
            let shuffled = shuffle_word_quantize(word, plan);
            let off = out_words * 4;
            output[off..off + 4].copy_from_slice(&shuffled.to_le_bytes());
            out_words += 1;
            bit_acc >>= 32; acc_bits -= 32;
        }
    }
    if acc_bits > 0 {
        // Zero-pad the tail to the next 32-bit word to keep raw word alignment.
        let word = (bit_acc & 0xFFFF_FFFF) as u32;
        let shuffled = shuffle_word_quantize(word, plan);
        let off = out_words * 4;
        output[off..off + 4].copy_from_slice(&shuffled.to_le_bytes());
        out_words += 1;
    }
    if out_words * 4 < output.len() {
        output.truncate(out_words * 4);
    }
    Ok(())
}
