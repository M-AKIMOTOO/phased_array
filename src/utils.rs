use std::f64::consts::PI;
use std::sync::Arc;

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
        Self { len, forward_r2c, inverse_c2r, inverse_c2c }
    }
    pub fn inverse_c2c(&self, spectrum: &mut [Complex<f64>]) -> Result<(), DynError> {
        if spectrum.len() != self.len { return Err("Spectrum length mismatch".into()); }
        self.inverse_c2c.process(spectrum);
        let scale = 1.0 / self.len as f64;
        for value in spectrum.iter_mut() { *value *= scale; }
        Ok(())
    }
    pub fn forward_r2c_process(&self, input: &mut [f64], output: &mut [Complex<f64>]) -> Result<(), DynError> {
        if input.len() != self.len || output.len() != self.len / 2 + 1 { return Err("Length mismatch".into()); }
        self.forward_r2c.process(input, output)?;
        Ok(())
    }
    pub fn inverse_c2r_process(&self, spectrum: &mut [Complex<f64>], output: &mut [f64]) -> Result<(), DynError> {
        if spectrum.len() != self.len / 2 + 1 || output.len() != self.len { return Err("Length mismatch".into()); }
        self.inverse_c2r.process(spectrum, output)?;
        let scale = 1.0 / self.len as f64;
        for value in output.iter_mut() { *value *= scale; }
        Ok(())
    }
}

pub fn apply_delay_and_rate_regular_bins(
    spectrum: &mut [Complex<f64>],
    fft_len: usize,
    freq_step_hz: f64,
    delay_seconds: f64,
    delay_rate: f64,
    geometric_accel: f64,
    frame_time: f64,
    is_lsb: bool,
) {
    if spectrum.is_empty() { return; }
    let total_delay_at_time = delay_seconds + delay_rate * frame_time + 0.5 * geometric_accel * frame_time.powi(2);
    // Standard VLBI sign: exp(-j * 2*PI * f * tau)
    // If we have flipped LSB to USB in time-domain, we use USB sign (-).
    // If not flipped, LSB uses positive sign (+).
    let sign = if is_lsb { 1.0 } else { -1.0 };
    let phase_step = sign * 2.0 * PI * freq_step_hz * total_delay_at_time;
    let mut rot = Complex::new(1.0, 0.0);
    let bin_rot_step = Complex::from_polar(1.0, phase_step);
    for bin in spectrum.iter_mut() { *bin *= rot; rot *= bin_rot_step; }
    if fft_len % 2 == 0 { let nyquist_idx = fft_len / 2; if nyquist_idx < spectrum.len() { spectrum[nyquist_idx].im = 0.0; } }
}

pub fn apply_integer_sample_shift_zerofill(samples: &mut [f64], shift_samples: i64) {
    if samples.is_empty() || shift_samples == 0 {
        return;
    }
    let n = samples.len();
    let s_abs = shift_samples.unsigned_abs() as usize;
    if s_abs >= n {
        samples.fill(0.0);
        return;
    }
    if shift_samples > 0 {
        // Positive delay: shift right (later in time), zero-fill head.
        samples.copy_within(..(n - s_abs), s_abs);
        samples[..s_abs].fill(0.0);
    } else {
        // Negative delay: shift left (earlier in time), zero-fill tail.
        samples.copy_within(s_abs.., 0);
        samples[(n - s_abs)..].fill(0.0);
    }
}

pub struct DecodePlan { bits: usize, shuffle_kind: ShuffleKind, input_shifts: [u32; 32] }
#[derive(Clone, Copy, Debug, PartialEq, Eq)] enum ShuffleKind { Identity, PairSwap, Generic }

pub fn build_decode_plan(bits: usize, shuffle_in: &[usize]) -> Result<DecodePlan, DynError> {
    let mut input_shifts = [0u32; 32];
    for (idx, &mapped) in shuffle_in.iter().enumerate() { input_shifts[idx] = mapped as u32; }
    let is_identity = shuffle_in.iter().enumerate().all(|(i, &v)| i == v);
    let is_pair_swap = shuffle_in.iter().enumerate().all(|(i, &v)| (i ^ 1) == v);
    let kind = if is_identity { ShuffleKind::Identity } else if is_pair_swap { ShuffleKind::PairSwap } else { ShuffleKind::Generic };
    Ok(DecodePlan { bits, shuffle_kind: kind, input_shifts })
}

pub fn decode_block_into_with_plan(raw: &[u8], levels: &[f64], samples: usize, plan: &DecodePlan, _bit_buf: &mut Vec<u8>, output: &mut [f64], lsb_to_usb: bool) -> Result<(), DynError> {
    let bits = plan.bits;
    let code_mask = (1u64 << bits) - 1;
    let mut out_idx = 0;
    let mut acc = 0u64;
    let mut acc_bits = 0;
    for chunk in raw.chunks_exact(4) {
        let mut word = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
        if plan.shuffle_kind == ShuffleKind::PairSwap { word = ((word & 0xAAAA_AAAA) >> 1) | ((word & 0x5555_5555) << 1); }
        else if plan.shuffle_kind == ShuffleKind::Generic {
            let mut shuffled = 0u32;
            for (out_bit, &in_shift) in plan.input_shifts.iter().enumerate() { shuffled |= ((word >> in_shift) & 1) << out_bit; }
            word = shuffled;
        }
        acc |= (word as u64) << acc_bits;
        acc_bits += 32;
        while acc_bits >= bits && out_idx < samples {
            let code = (acc & code_mask) as usize;
            let mut val = levels[code];
            // Spectral flip: multiply odd samples by -1
            if lsb_to_usb && (out_idx & 1) == 1 { val = -val; }
            output[out_idx] = val;
            out_idx += 1; acc_bits -= bits; acc >>= bits;
        }
    }
    if out_idx != samples {
        return Err(format!("decoded {} samples, expected {}", out_idx, samples).into());
    }
    Ok(())
}

pub fn quantise_frame(samples: &[f64], bits: usize, levels: &[f64], shuffle_out: &[usize], output: &mut Vec<u8>) -> Result<(), DynError> {
    output.clear();
    let code_mask = (1u64 << bits) - 1;
    let mut bit_acc = 0u64;
    let mut acc_bits = 0;
    for &val in samples.iter() {
        let mut best_idx = 0;
        let mut min_err = f64::MAX;
        for (idx, &lv) in levels.iter().enumerate() {
            let err = (val - lv).abs();
            if err < min_err { min_err = err; best_idx = idx; }
        }
        bit_acc |= (best_idx as u64 & code_mask) << acc_bits;
        acc_bits += bits;
        while acc_bits >= 32 {
            let word = (bit_acc & 0xFFFF_FFFF) as u32;
            let mut shuffled = 0u32;
            for (pos, &target) in shuffle_out.iter().enumerate() { shuffled |= ((word >> pos) & 1) << target; }
            output.extend_from_slice(&shuffled.to_le_bytes());
            bit_acc >>= 32; acc_bits -= 32;
        }
    }
    if acc_bits > 0 {
        // Zero-pad the tail to the next 32-bit word to keep raw word alignment.
        let word = (bit_acc & 0xFFFF_FFFF) as u32;
        let mut shuffled = 0u32;
        for (pos, &target) in shuffle_out.iter().enumerate() {
            shuffled |= ((word >> pos) & 1) << target;
        }
        output.extend_from_slice(&shuffled.to_le_bytes());
    }
    Ok(())
}

pub fn unwrap_phase(phase: &[f64]) -> Vec<f64> {
    let mut unwrapped = Vec::with_capacity(phase.len());
    if let Some(&first) = phase.first() {
        unwrapped.push(first);
        let mut offset = 0.0;
        for i in 1..phase.len() {
            let diff = phase[i] - phase[i - 1];
            if diff > PI { offset -= 2.0 * PI; } else if diff < -PI { offset += 2.0 * PI; }
            unwrapped.push(phase[i] + offset);
        }
    }
    unwrapped
}
pub fn safe_arg(z: &Complex<f64>) -> f64 { if z.re == 0.0 && z.im == 0.0 { 0.0 } else { z.arg() } }
