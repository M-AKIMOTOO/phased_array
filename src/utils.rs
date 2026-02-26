use std::f64::consts::PI;
use std::sync::Arc;

use rayon::prelude::*;
use realfft::{ComplexToReal, RealFftPlanner, RealToComplex};
use rustfft::{num_complex::Complex, Fft, FftPlanner};
use std::error::Error;

pub type DynError = Box<dyn Error + Send + Sync>;

pub struct FftHelper {
    len: usize,
    pub forward_r2c: Arc<dyn RealToComplex<f64>>,
    pub inverse_c2r: Arc<dyn ComplexToReal<f64>>,
    pub inverse_c2c: Arc<dyn Fft<f64>>,
}

impl FftHelper {
    pub fn new(len: usize) -> Self {
        let mut planner_c2c = FftPlanner::new();
        let mut planner_r2c = RealFftPlanner::new();

        let forward_r2c = planner_r2c.plan_fft_forward(len);
        let inverse_c2r = planner_r2c.plan_fft_inverse(len);
        let inverse_c2c = planner_c2c.plan_fft_inverse(len);

        Self {
            len,
            forward_r2c,
            inverse_c2r,
            inverse_c2c,
        }
    }

    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn inverse_c2c(&self, spectrum: &mut [Complex<f64>]) -> Result<(), DynError> {
        if spectrum.len() != self.len {
            return Err("Spectrum length does not match FFT configuration for C2C".into());
        }
        self.inverse_c2c.process(spectrum);
        let scale = 1.0 / self.len as f64;
        for value in spectrum.iter_mut() {
            *value *= scale;
        }
        Ok(())
    }

    pub fn forward_r2c_process(
        &self,
        input: &mut [f64],
        output: &mut [Complex<f64>],
    ) -> Result<(), DynError> {
        if input.len() != self.len {
            return Err("Input length for R2C does not match FFT configuration".into());
        }
        if output.len() != self.len / 2 + 1 {
            return Err(
                "Output length for R2C does not match expected half-spectrum length".into(),
            );
        }
        self.forward_r2c.process(input, output)?;
        Ok(())
    }

    pub fn inverse_c2r_process(
        &self,
        spectrum: &mut [Complex<f64>],
        output: &mut [f64],
    ) -> Result<(), DynError> {
        if spectrum.len() != self.len / 2 + 1 {
            return Err(
                "Input spectrum length for C2R does not match expected half-spectrum length".into(),
            );
        }
        if output.len() != self.len {
            return Err(
                "Output buffer length for C2R does not match expected time-domain length".into(),
            );
        }
        self.inverse_c2r.process(spectrum, output)?;
        let scale = 1.0 / self.len as f64;
        for value in output.iter_mut() {
            *value *= scale;
        }
        Ok(())
    }
}

#[allow(dead_code)]
pub fn apply_delay_and_rate(
    spectrum: &mut [Complex<f64>],
    sample_rate: f64,
    fft_len: usize,
    delay_seconds: f64,
    delay_rate: f64,
    geometric_accel: f64,
    frame_time: f64,
) {
    let freq_scale = sample_rate / fft_len as f64;
    let mut freqs_hz = Vec::with_capacity(spectrum.len());
    for idx in 0..spectrum.len() {
        freqs_hz.push(idx as f64 * freq_scale);
    }
    apply_delay_and_rate_with_freqs(
        spectrum,
        &freqs_hz,
        fft_len,
        delay_seconds,
        delay_rate,
        geometric_accel,
        frame_time,
    );
}

pub fn apply_delay_and_rate_with_freqs(
    spectrum: &mut [Complex<f64>],
    freqs_hz: &[f64],
    fft_len: usize,
    delay_seconds: f64,
    delay_rate: f64,
    geometric_accel: f64,
    frame_time: f64,
) {
    debug_assert_eq!(spectrum.len(), freqs_hz.len());
    let total_delay_at_time =
        delay_seconds + delay_rate * frame_time + 0.5 * geometric_accel * frame_time.powi(2);

    for (bin, &freq) in spectrum.iter_mut().zip(freqs_hz.iter()) {
        let phase_shift = -2.0 * PI * freq * total_delay_at_time;
        let rot = Complex::from_polar(1.0, phase_shift);
        *bin *= rot;
    }

    // For real output from C2R IFFT, the Nyquist frequency component (if fft_len is even)
    // must be purely real. Ensure its imaginary part is zero after phase rotation.
    if fft_len % 2 == 0 && !spectrum.is_empty() {
        let nyquist_idx = fft_len / 2;
        if nyquist_idx < spectrum.len() {
            spectrum[nyquist_idx].im = 0.0;
        }
    }
}

/// Apply delay/rate phase rotation for uniformly spaced FFT bins starting at DC.
/// This avoids per-bin sin/cos by using a complex phase recurrence.
pub fn apply_delay_and_rate_regular_bins(
    spectrum: &mut [Complex<f64>],
    fft_len: usize,
    freq_step_hz: f64,
    delay_seconds: f64,
    delay_rate: f64,
    geometric_accel: f64,
    frame_time: f64,
) {
    if spectrum.is_empty() {
        return;
    }
    let total_delay_at_time =
        delay_seconds + delay_rate * frame_time + 0.5 * geometric_accel * frame_time.powi(2);
    let phase_step = -2.0 * PI * freq_step_hz * total_delay_at_time;
    rotate_regular_bins_with_step(spectrum, fft_len, Complex::from_polar(1.0, phase_step));
}

pub fn rotate_regular_bins_with_step(
    spectrum: &mut [Complex<f64>],
    fft_len: usize,
    bin_rot_step: Complex<f64>,
) {
    if spectrum.is_empty() {
        return;
    }
    let mut rot = Complex::new(1.0, 0.0);
    for bin in spectrum.iter_mut() {
        *bin *= rot;
        rot *= bin_rot_step;
    }

    // For real output from C2R IFFT, the Nyquist frequency component (if fft_len is even)
    // must be purely real. Ensure its imaginary part is zero after phase rotation.
    if fft_len % 2 == 0 {
        let nyquist_idx = fft_len / 2;
        if nyquist_idx < spectrum.len() {
            spectrum[nyquist_idx].im = 0.0;
        }
    }
}

const PARALLEL_DECODE_THRESHOLD_BYTES: usize = 64 * 1024;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ShuffleKind {
    Identity,
    PairSwap,
    Generic,
}

pub struct DecodePlan {
    bits: usize,
    shuffle_kind: ShuffleKind,
    input_shifts: [u32; 32],
}

#[inline(always)]
fn classify_shuffle(shuffle_in: &[usize]) -> ShuffleKind {
    if shuffle_in
        .iter()
        .enumerate()
        .all(|(idx, &value)| idx == value)
    {
        ShuffleKind::Identity
    } else if shuffle_in
        .iter()
        .enumerate()
        .all(|(idx, &value)| (idx ^ 1) == value)
    {
        ShuffleKind::PairSwap
    } else {
        ShuffleKind::Generic
    }
}

#[inline(always)]
fn swap_adjacent_bits_u8(byte: u8) -> u8 {
    ((byte & 0xAA) >> 1) | ((byte & 0x55) << 1)
}

#[inline(always)]
fn swap_adjacent_bits_u32(word: u32) -> u32 {
    ((word & 0xAAAA_AAAA) >> 1) | ((word & 0x5555_5555) << 1)
}

#[inline(always)]
fn apply_shuffle_word(word: u32, shuffle_kind: ShuffleKind, input_shifts: &[u32; 32]) -> u32 {
    match shuffle_kind {
        ShuffleKind::Identity => word,
        ShuffleKind::PairSwap => swap_adjacent_bits_u32(word),
        ShuffleKind::Generic => {
            let mut shuffled = 0u32;
            for (out_bit, &in_shift) in input_shifts.iter().enumerate() {
                let bit = (word >> in_shift) & 1;
                shuffled |= bit << out_bit as u32;
            }
            shuffled
        }
    }
}

#[inline(always)]
fn decode_2bit_byte_to_levels(byte: u8, levels: &[f64], lsb_to_usb: bool, out4: &mut [f64]) {
    // gico3-compatible order for 2bit/1ch:
    // sample0 <- bits[1:0], sample1 <- bits[3:2], sample2 <- bits[5:4], sample3 <- bits[7:6]
    let code0 = (byte & 0b11) as usize;
    let mut code1 = ((byte >> 2) & 0b11) as usize;
    let code2 = ((byte >> 4) & 0b11) as usize;
    let mut code3 = ((byte >> 6) & 0b11) as usize;
    if lsb_to_usb {
        // Multiplication by (-1)^n in 2-bit code space: odd samples toggle both bits.
        code1 ^= 0b11;
        code3 ^= 0b11;
    }
    out4[0] = levels[code0];
    out4[1] = levels[code1];
    out4[2] = levels[code2];
    out4[3] = levels[code3];
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn decode_2bit_pair_swap_avx2(
    raw: &[u8],
    levels: &[f64],
    lsb_to_usb: bool,
    output: &mut [f64],
    sample_limit: usize,
    mut byte_idx: usize,
    mut out_idx: usize,
) -> (usize, usize) {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    let mask_aa = _mm256_set1_epi8(0xAAu8 as i8);
    let mask_55 = _mm256_set1_epi8(0x55u8 as i8);

    while byte_idx + 32 <= raw.len() && out_idx + 128 <= sample_limit {
        let packed = _mm256_loadu_si256(raw.as_ptr().add(byte_idx) as *const __m256i);
        let hi = _mm256_srli_epi16(_mm256_and_si256(packed, mask_aa), 1);
        let lo = _mm256_slli_epi16(_mm256_and_si256(packed, mask_55), 1);
        let swapped = _mm256_or_si256(hi, lo);

        let mut bytes = [0u8; 32];
        _mm256_storeu_si256(bytes.as_mut_ptr() as *mut __m256i, swapped);
        for &byte in &bytes {
            decode_2bit_byte_to_levels(byte, levels, lsb_to_usb, &mut output[out_idx..out_idx + 4]);
            out_idx += 4;
        }
        byte_idx += 32;
    }
    (byte_idx, out_idx)
}

#[allow(dead_code)]
pub fn decode_block(
    raw: &[u8],
    bits: usize,
    levels: &[f64],
    samples: usize,
    shuffle_in: &[usize],
    result: &mut Vec<f64>,
) -> Result<(), DynError> {
    result.resize(samples, 0.0);
    let mut buffer = Vec::new();
    decode_block_into(
        raw,
        bits,
        levels,
        samples,
        shuffle_in,
        &mut buffer,
        result.as_mut_slice(),
        false,
    )
}

pub fn decode_block_into(
    raw: &[u8],
    bits: usize,
    levels: &[f64],
    samples: usize,
    shuffle_in: &[usize],
    _bit_buffer: &mut Vec<u8>,
    output: &mut [f64],
    lsb_to_usb: bool,
) -> Result<(), DynError> {
    let plan = build_decode_plan(bits, shuffle_in)?;
    decode_block_into_with_plan(raw, levels, samples, &plan, _bit_buffer, output, lsb_to_usb)
}

pub fn build_decode_plan(bits: usize, shuffle_in: &[usize]) -> Result<DecodePlan, DynError> {
    if bits == 0 {
        return Err("Bit depth must be at least 1".into());
    }
    if bits > 32 {
        return Err("Bit depth above 32 is not supported".into());
    }
    if shuffle_in.len() != 32 {
        return Err("Shuffle map must contain 32 entries".into());
    }
    if shuffle_in.iter().any(|&idx| idx >= 32) {
        return Err("Shuffle map entries must be between 0 and 31".into());
    }

    let shuffle_kind = classify_shuffle(shuffle_in);
    let mut input_shifts = [0u32; 32];
    for (idx, &mapped) in shuffle_in.iter().enumerate() {
        input_shifts[idx] = mapped as u32;
    }

    Ok(DecodePlan {
        bits,
        shuffle_kind,
        input_shifts,
    })
}

pub fn decode_block_into_with_plan(
    raw: &[u8],
    levels: &[f64],
    samples: usize,
    plan: &DecodePlan,
    _bit_buffer: &mut Vec<u8>,
    output: &mut [f64],
    lsb_to_usb: bool,
) -> Result<(), DynError> {
    if output.len() != samples {
        return Err("Output slice length does not match expected sample count".into());
    }
    if raw.len() % 4 != 0 {
        return Err("Raw block length is not a multiple of 4 bytes".into());
    }
    let bits = plan.bits;
    let required_levels = 1usize
        .checked_shl(bits as u32)
        .ok_or("Bit depth too large for level table")?;
    if levels.len() < required_levels {
        return Err(format!(
            "Expected at least {required_levels} quantisation levels for {bits} bits, received {}",
            levels.len()
        )
        .into());
    }

    let shuffle_kind = plan.shuffle_kind;

    // 2bit + (identity | adjacent-pair swap) is by far the common hot path.
    // Decode directly from bytes, and optionally parallelise for very large buffers.
    if bits == 2 && matches!(shuffle_kind, ShuffleKind::Identity | ShuffleKind::PairSwap) {
        if samples == raw.len() * 4
            && raw.len() >= PARALLEL_DECODE_THRESHOLD_BYTES
            && rayon::current_thread_index().is_none()
        {
            if shuffle_kind == ShuffleKind::Identity {
                output
                    .par_chunks_mut(4)
                    .zip(raw.par_iter().copied())
                    .for_each(|(out4, byte)| {
                        decode_2bit_byte_to_levels(byte, levels, lsb_to_usb, out4)
                    });
            } else {
                output
                    .par_chunks_mut(4)
                    .zip(raw.par_iter().copied())
                    .for_each(|(out4, byte)| {
                        decode_2bit_byte_to_levels(
                            swap_adjacent_bits_u8(byte),
                            levels,
                            lsb_to_usb,
                            out4,
                        )
                    });
            }
            return Ok(());
        }

        let mut out_idx = 0usize;
        let mut byte_idx = 0usize;

        if shuffle_kind == ShuffleKind::PairSwap {
            #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
            {
                if std::arch::is_x86_feature_detected!("avx2") {
                    let (b, o) = unsafe {
                        decode_2bit_pair_swap_avx2(
                            raw, levels, lsb_to_usb, output, samples, byte_idx, out_idx,
                        )
                    };
                    byte_idx = b;
                    out_idx = o;
                }
            }
        }

        while byte_idx < raw.len() && out_idx + 4 <= samples {
            let byte = raw[byte_idx];
            let decoded_byte = if shuffle_kind == ShuffleKind::PairSwap {
                swap_adjacent_bits_u8(byte)
            } else {
                byte
            };
            decode_2bit_byte_to_levels(
                decoded_byte,
                levels,
                lsb_to_usb,
                &mut output[out_idx..out_idx + 4],
            );
            byte_idx += 1;
            out_idx += 4;
        }

        if out_idx != samples {
            return Err("Decoded sample count did not match expected length".into());
        }
        return Ok(());
    }

    let input_shifts = &plan.input_shifts;

    if 32 % bits == 0 {
        let mut out_idx = 0usize;
        let samples_per_word = 32 / bits;
        let code_mask = if bits == 32 {
            u32::MAX
        } else {
            (1u32 << bits as u32) - 1
        };

        for chunk in raw.chunks_exact(4) {
            if out_idx >= samples {
                break;
            }
            let word = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
            let mut shuffled = apply_shuffle_word(word, shuffle_kind, input_shifts);

            for _ in 0..samples_per_word {
                if out_idx >= samples {
                    break;
                }
                let code = if bits == 32 {
                    shuffled as usize
                } else {
                    (shuffled & code_mask) as usize
                };
                let mut value = levels[code];
                if lsb_to_usb && (out_idx & 1) == 1 {
                    value = -value;
                }
                output[out_idx] = value;
                out_idx += 1;
                if bits < 32 {
                    shuffled >>= bits as u32;
                }
            }
        }

        if out_idx != samples {
            return Err("Decoded sample count did not match expected length".into());
        }
        return Ok(());
    }

    // Non-divisor bit depths (e.g. 3,5,6bit): keep a tiny bit accumulator instead of
    // allocating/draining a dynamic bit buffer.
    let code_mask = (1u64 << bits as u32) - 1;
    let mut out_idx = 0usize;
    let mut acc = 0u64;
    let mut acc_bits = 0usize;

    for chunk in raw.chunks_exact(4) {
        let word = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
        let shuffled = apply_shuffle_word(word, shuffle_kind, input_shifts);
        // Non-divisor depths still consume bitstream from LSB to MSB.
        acc |= (shuffled as u64) << acc_bits as u32;
        acc_bits += 32;

        while acc_bits >= bits && out_idx < samples {
            let code = (acc & code_mask) as usize;
            let mut value = levels[code];
            if lsb_to_usb && (out_idx & 1) == 1 {
                value = -value;
            }
            output[out_idx] = value;
            out_idx += 1;
            acc_bits -= bits;
            acc >>= bits as u32;
        }
    }

    if out_idx != samples {
        return Err("Decoded sample count did not match expected length".into());
    }
    Ok(())
}

#[allow(dead_code)]
pub fn quantise_frame(
    samples: &[f64],
    bits: usize,
    levels: &[f64],
    shuffle_out: &[usize],
    output: &mut Vec<u8>,
) -> Result<(), DynError> {
    if bits == 0 {
        return Err("Bit depth must be at least 1".into());
    }
    if shuffle_out.len() != 32 {
        return Err("Shuffle map must contain 32 entries".into());
    }
    if levels.is_empty() {
        return Err("Quantisation levels not provided".into());
    }
    let bits_per_frame = samples
        .len()
        .checked_mul(bits)
        .ok_or("Overflow when computing bits per frame")?;
    if bits_per_frame % 32 != 0 {
        return Err("Frame bit-length must be divisible by 32".into());
    }
    output.clear();
    output.reserve(bits_per_frame / 8);

    let code_mask = if bits == 64 {
        u64::MAX
    } else {
        (1u64 << bits as u32) - 1
    };
    let mut bit_acc = 0u64;
    let mut acc_bits = 0usize;

    for &value in samples {
        let mut best_idx = 0usize;
        let mut best_error = (value - levels[0]).abs();
        for (idx, level) in levels.iter().enumerate().skip(1) {
            let error = (value - *level).abs();
            if error < best_error {
                best_error = error;
                best_idx = idx;
            }
        }
        let code = (best_idx as u64) & code_mask;
        bit_acc |= code << acc_bits as u32;
        acc_bits += bits;

        while acc_bits >= 32 {
            let word = (bit_acc & 0xFFFF_FFFF) as u32;
            let mut shuffled_word = 0u32;
            for (pos, &target) in shuffle_out.iter().enumerate() {
                let bit = (word >> pos as u32) & 1;
                shuffled_word |= bit << target as u32;
            }
            output.extend_from_slice(&shuffled_word.to_le_bytes());
            bit_acc >>= 32;
            acc_bits -= 32;
        }
    }
    debug_assert_eq!(acc_bits, 0);
    Ok(())
}

/// 累積バッファに複素スペクトルのパワー `|z|^2` を加算する。
/// 実行時に利用可能な SIMD 命令を検出して高速パスを選択する。
pub fn accumulate_power_add(dest: &mut [f64], src: &[Complex<f64>]) {
    debug_assert_eq!(dest.len(), src.len());
    let len = dest.len();
    if len == 0 {
        return;
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if len >= 2 && std::arch::is_x86_feature_detected!("avx2") {
            unsafe {
                return accumulate_power_add_avx2(dest, src);
            }
        }
    }

    // フォールバック: スカラー計算
    for (acc, value) in dest.iter_mut().zip(src.iter()) {
        *acc += value.re * value.re + value.im * value.im;
    }
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "avx2")]
unsafe fn accumulate_power_add_avx2(dest: &mut [f64], src: &[Complex<f64>]) {
    use std::arch::x86_64::*;

    let len = dest.len();
    let mut i = 0usize;
    while i + 1 < len {
        let ptr = src.as_ptr().add(i) as *const f64;
        let values = _mm256_loadu_pd(ptr);
        let squared = _mm256_mul_pd(values, values);
        let sums = _mm256_hadd_pd(squared, squared);
        let arr: [f64; 4] = std::mem::transmute(sums);
        dest[i] += arr[0];
        dest[i + 1] += arr[1];
        i += 2;
    }
    while i < len {
        let value = src[i];
        dest[i] += value.re * value.re + value.im * value.im;
        i += 1;
    }
}

pub fn unwrap_phase(phase: &[f64]) -> Vec<f64> {
    let mut unwrapped = Vec::with_capacity(phase.len());
    if let Some(&first) = phase.first() {
        unwrapped.push(first);
        let mut offset = 0.0;
        for i in 1..phase.len() {
            let diff = phase[i] - phase[i - 1];
            if diff > PI {
                offset -= 2.0 * PI;
            } else if diff < -PI {
                offset += 2.0 * PI;
            }
            unwrapped.push(phase[i] + offset);
        }
    }
    unwrapped
}

#[allow(dead_code)]
pub fn hanning_window(len: usize) -> Vec<f64> {
    let mut window = vec![0.0; len];
    for i in 0..len {
        window[i] = 0.5 * (1.0 - (2.0 * PI * i as f64 / (len as f64 - 1.0)).cos());
    }
    window
}
pub fn safe_arg(z: &Complex<f64>) -> f64 {
    if z.re == 0.0 && z.im == 0.0 {
        0.0
    } else {
        z.arg()
    }
}

#[allow(dead_code)]
pub fn inverse_permutation(p: &[usize]) -> Result<Vec<usize>, DynError> {
    if p.len() != 32 {
        return Err("Permutation must contain exactly 32 entries".into());
    }
    let mut inv = vec![0usize; 32];
    for (idx, &value) in p.iter().enumerate() {
        if value >= 32 {
            return Err("Permutation entries must be between 0 and 31".into());
        }
        inv[value] = idx;
    }
    Ok(inv)
}

#[cfg(test)]
mod tests {
    use super::decode_block_into;

    fn pack_codes_lsb_first_2bit(codes: &[usize]) -> Vec<u8> {
        assert_eq!(
            codes.len(),
            16,
            "2bit x 16 samples must form one 32-bit word"
        );
        let mut word = 0u32;
        for (idx, &code) in codes.iter().enumerate() {
            assert!(code < 4);
            word |= (code as u32) << (2 * idx);
        }
        word.to_le_bytes().to_vec()
    }

    #[test]
    fn decode_2bit_identity_keeps_00_01_10_11_order() {
        let codes = [0usize, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let raw = pack_codes_lsb_first_2bit(&codes);
        let levels = [-1.5f64, -0.5, 0.5, 1.5];
        let shuffle: Vec<usize> = (0..32).collect();

        let mut bit_buffer = Vec::new();
        let mut output = vec![0.0; 16];
        decode_block_into(
            &raw,
            2,
            &levels,
            16,
            &shuffle,
            &mut bit_buffer,
            &mut output,
            false,
        )
        .unwrap();

        let expected: Vec<f64> = codes.iter().map(|&c| levels[c]).collect();
        assert_eq!(output, expected);
    }

    #[test]
    fn level_permutation_matches_bit_swap_for_2bit() {
        let codes = [0usize, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let raw = pack_codes_lsb_first_2bit(&codes);
        let identity_shuffle: Vec<usize> = (0..32).collect();
        let bit_swap_shuffle: Vec<usize> = (0..32).step_by(2).flat_map(|i| [i + 1, i]).collect();
        let canonical_levels = [-1.5f64, -0.5, 0.5, 1.5];
        let swapped_levels = [-1.5f64, 0.5, -0.5, 1.5];

        let mut buf_a = Vec::new();
        let mut buf_b = Vec::new();
        let mut out_a = vec![0.0; 16];
        let mut out_b = vec![0.0; 16];

        decode_block_into(
            &raw,
            2,
            &swapped_levels,
            16,
            &identity_shuffle,
            &mut buf_a,
            &mut out_a,
            false,
        )
        .unwrap();
        decode_block_into(
            &raw,
            2,
            &canonical_levels,
            16,
            &bit_swap_shuffle,
            &mut buf_b,
            &mut out_b,
            false,
        )
        .unwrap();

        assert_eq!(out_a, out_b);
    }

    #[test]
    fn lsb_to_usb_flip_on_odd_samples() {
        let codes = [0usize, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let raw = pack_codes_lsb_first_2bit(&codes);
        let levels = [-1.5f64, -0.5, 0.5, 1.5];
        let shuffle: Vec<usize> = (0..32).collect();
        let mut bit_buffer = Vec::new();
        let mut output = vec![0.0; 16];

        decode_block_into(
            &raw,
            2,
            &levels,
            16,
            &shuffle,
            &mut bit_buffer,
            &mut output,
            true,
        )
        .unwrap();

        let expected: Vec<f64> = codes
            .iter()
            .enumerate()
            .map(|(idx, &c)| {
                let v = levels[c];
                if (idx & 1) == 1 {
                    -v
                } else {
                    v
                }
            })
            .collect();
        assert_eq!(output, expected);
    }
}
