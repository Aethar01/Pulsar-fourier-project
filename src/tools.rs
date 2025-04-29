use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::f64::consts::PI;

use crate::ft;

pub fn read_data_file<P>(filename: P) -> io::Result<Vec<f64>>
where P: AsRef<Path> {
    let file = File::open(filename)?;
    let reader = io::BufReader::new(file);
    
    let mut data = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if let Ok(value) = line.trim().parse::<f64>() {
            data.push(value);
        }
    }
    
    Ok(data)
}

pub fn generate_synthetic_data(freq1: f64, freq2: f64) -> Vec<f64> {
    let n = 256;
    let mut data = Vec::with_capacity(n);
    
    for t in 0..n {
        let t_f64 = t as f64;
        let value = 10.0 * (2.0 * PI * freq1 * t_f64 / n as f64).cos() + 
                   5.0 * (2.0 * PI * freq2 * t_f64 / n as f64).sin();
        data.push(value);
    }
    
    data
}

pub fn analyze_pulsar_data<W: Write>(data: &[f64], mut output: W) {
    let n = data.len();
    let delta_t = 0.004; // 4ms sampling interval
    
    // Compute Fourier coefficients and power spectrum
    let (_ak, _bk, power) = ft::compute_fourier_transform(data);
    
    // Find peaks in power spectrum
    let peaks = find_peaks(&power);
    
    writeln!(output, "Fourier Analysis Results").unwrap();
    writeln!(output, "======================").unwrap();
    writeln!(output, "Data points: {}", n).unwrap();
    writeln!(output, "Sampling interval: {} s", delta_t).unwrap();
    writeln!(output, "\nFrequency (Hz)\tPower").unwrap();
    
    for (k, &p) in power.iter().enumerate() {
        let freq = k as f64 / (n as f64 * delta_t);
        writeln!(output, "{:.3}\t{:.3}", freq, p).unwrap();
    }
    
    writeln!(output, "\nSignificant Peaks:").unwrap();
    for peak in &peaks {
        let freq = peak.index as f64 / (n as f64 * delta_t);
        writeln!(output, "Frequency: {:.3} Hz, Power: {:.3}", freq, peak.value).unwrap();
    }
    
    // Phase binning analysis for the strongest peak
    if let Some(main_peak) = peaks.first() {
        let freq = main_peak.index as f64 / (n as f64 * delta_t);
        let period = 1.0 / freq;
        
        writeln!(output, "\nPhase Binning Analysis for {:.3} Hz (Period: {:.3} s)", freq, period).unwrap();
        
        let binned_data = phase_binning(data, period, delta_t, 10);
        writeln!(output, "\nPhase Bins:").unwrap();
        for (i, &bin) in binned_data.iter().enumerate() {
            writeln!(output, "Bin {}: {:.3}", i, bin).unwrap();
        }
    }
}

struct Peak {
    index: usize,
    value: f64,
}

fn find_peaks(power: &[f64]) -> Vec<Peak> {
    let mut peaks = Vec::new();
    
    // Skip DC component (k=0) and Nyquist frequency
    for k in 1..power.len() / 2 {
        if power[k] > power[k - 1] && power[k] > power[k + 1] {
            peaks.push(Peak {
                index: k,
                value: power[k],
            });
        }
    }
    
    // Sort peaks by power in descending order
    peaks.sort_by(|a, b| b.value.partial_cmp(&a.value).unwrap());
    
    peaks
}

fn phase_binning(data: &[f64], period: f64, delta_t: f64, num_bins: usize) -> Vec<f64> {
    let mut bins = vec![0.0; num_bins];
    let mut counts = vec![0; num_bins];
    
    for (t, &value) in data.iter().enumerate() {
        let time = t as f64 * delta_t;
        let phase = (time / period).fract();
        let bin_index = (phase * num_bins as f64).floor() as usize;
        
        if bin_index < num_bins {
            bins[bin_index] += value;
            counts[bin_index] += 1;
        }
    }
    
    // Normalize bins by counts
    for i in 0..num_bins {
        if counts[i] > 0 {
            bins[i] /= counts[i] as f64;
        }
    }
    
    bins
}
