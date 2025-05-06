use std::fs::File;
use std::io::{self, BufRead, Write};
use std::path::Path;
use std::f64::consts::PI;
use std::process::Command;
use std::path::PathBuf;

use crate::ft;

fn run_gnuplot(script_path: &Path) -> std::io::Result<()> {
    Command::new("gnuplot")
        .arg(script_path)
        .status()
        .map(|_| ())
}

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
    let delta_t = 0.004;
    
    for t in 0..n {
        let t_f64 = t as f64 * delta_t;
        let value = 10.0 * (2.0 * PI * freq1 * t_f64 as f64).cos() + 
                   10.0 * (2.0 * PI * freq2 * t_f64 as f64).sin();
        data.push(value);
    }
    
    data
}

pub struct Data {
    raw_data: Vec<f64>,
    ak: Vec<f64>,
    bk: Vec<f64>,
    power: Vec<f64>,
    peaks: Vec<Peak>,
    harmonics: Vec<usize>,
    n: usize,
    delta_t: f64,
    bins: Vec<f64>,
    phases_of_bins: Vec<f64>,
}

pub fn analyze_pulsar_data<W: Write>(data: Vec<f64>, mut output: W, synthetic: Option<(f64, f64)>) -> Data {
    let n = data.len();
    let delta_t = 0.004; // 4ms sampling interval
    
    // Compute Fourier coefficients and power spectrum
    let (ak, bk, power) = ft::compute_fourier_transform(&data);
    
    // Estimate noise floor
    let noise_floor = estimate_noise_floor(&power);

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
    let delta_f = 1.0 / (n as f64 * delta_t);
    let freq_error = delta_f / 2.0;
    for peak in &peaks {
        let freq = peak.index as f64 * delta_f;
        let snr = (peak.value - noise_floor) / noise_floor;
        let snr_error = (peak.value / noise_floor).sqrt();
        writeln!(output, "Frequency: {:.3} Hz, Power: {:.3}, SNR: {:.3} ± {:.3}", freq, peak.value, snr, snr_error).unwrap();
    }
    writeln!(output, "Frequency Error: +- {:.3} Hz", freq_error).unwrap();
    
    let (bins, phases_of_bins): (Vec<f64>, Vec<f64>) = {
            if let Some(main_peak) = peaks.get(0) {
                let freq = main_peak.index as f64 / (n as f64 * delta_t);
                let period = 1.0 / freq;
                
                writeln!(output, "\nPhase Binning Analysis for {:.3} Hz (Period: {:.3} s)", freq, period).unwrap();
                
                let binsnum = 11;
                let mean = data.iter().sum::<f64>() / (n as f64);
                let detrended = data.iter().map(|x| x - mean).collect::<Vec<f64>>();
                let (binned_data, _) = phase_binning(&data, period, delta_t, binsnum);
                let (binned_detrended_data, phase_of_bins) = phase_binning(&detrended, period, delta_t, binsnum);
                writeln!(output, "\nDetrended Phase Bins:").unwrap();
                for (i, &bin) in binned_detrended_data.iter().enumerate() {
                    writeln!(output, "Bin {}: {:.3}", i, bin).unwrap();
                }
                writeln!(output, "\nSum of Detrended Phase Bins: {}", binned_detrended_data.iter().sum::<f64>()).unwrap();
                (binned_data, phase_of_bins)
            }
            else {
                (Vec::new(), Vec::new())
            }
    };

    // if let Some(main_peak) = peaks.first() {
    //     let freq = main_peak.index as f64 / (n as f64 * delta_t);
    //     let period = 1.0 / freq;
        
    //     writeln!(output, "\nPhase Binning Analysis for {:.3} Hz (Period: {:.3} s)", freq, period).unwrap();
        
    //     let binned_data = phase_binning(data, period, delta_t, 10);
    //     writeln!(output, "\nPhase Bins:").unwrap();
    //     for (i, &bin) in binned_data.iter().enumerate() {
    //         writeln!(output, "Bin {}: {:.3}", i, bin).unwrap();
    //     }
    //     writeln!(output, "\nSum of Phase Bins: {}", binned_data.iter().sum::<f64>()).unwrap();
    // }

    writeln!(output, "\nHarmonic Analysis:").unwrap();
    let harmonics = {
        if let Some(main_peak) = peaks.first() {
            let harmonics = detect_harmonics(&peaks, main_peak.index);
            writeln!(output, "Fundamental Frequency: {:.3} Hz", main_peak.index as f64 / (n as f64 * delta_t)).unwrap();
            for h in &harmonics {
                let freq = *h as f64 / (n as f64 * delta_t);
                let ratio = freq / main_peak.index as f64;
                let error = (ratio.round() - ratio).abs();
                writeln!(output, "Harmonic at {:.3} Hz ({:.3}x fundamental, +-{:.3})", freq, ratio, error).unwrap();
            }
            harmonics
        } else {
            Vec::new()
        }
    };

    writeln!(output, "\nBest Period Detection:").unwrap();
    if let Some(main_peak) = peaks.first() {
        let (best_period, error) = find_best_period(&data, 1.0 / main_peak.index as f64, delta_t, 30);
        writeln!(output, "Fundamental Frequency: {:.3} Hz", main_peak.index as f64 / (n as f64 * delta_t)).unwrap();
        writeln!(output, "Best Period: {:.3} +- {:.5} s", best_period, error).unwrap();
    }

    if let Some((freq1, freq2)) = synthetic {
        writeln!(output, "\nSynthetic Data Validation:").unwrap();
        writeln!(output, "Fundamental Frequency: {:.3} Hz", freq1).unwrap();
        writeln!(output, "Second Fundamental Frequency: {:.3} Hz", freq2).unwrap();
        writeln!(output, "Validation: {}", validate_synthetic_data(freq1, freq2, &peaks, freq_error, delta_f)).unwrap();
    }

    Data {
        raw_data: data,
        ak,
        bk,
        power,
        peaks,
        harmonics,
        n,
        delta_t,
        bins,
        phases_of_bins,
    }
}

pub struct Peak {
    pub index: usize,
    pub value: f64,
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

fn phase_binning(data: &[f64], period: f64, delta_t: f64, num_bins: usize) -> (Vec<f64>, Vec<f64>) {
    let mut bins = vec![0.0; num_bins];
    let mut counts = vec![0; num_bins];
    let phase_of_bins = (0..num_bins).map(|i| (i as f64 / num_bins as f64) * 2.0 * PI).collect::<Vec<f64>>();
    
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
    bins.iter_mut()
        .zip(counts.iter())
        .for_each(|(bin, &count)| {
            *bin = if count > 0 { *bin / count as f64 } else { f64::NAN };
        });
    
    (bins, phase_of_bins)
}

pub fn validate_synthetic_data(freq1: f64, freq2: f64, peaks: &[Peak], freq_error: f64, delta_f: f64) -> bool {
    peaks.iter().any(|p| (p.index as f64 * delta_f - freq1).abs() < freq_error) &&
    peaks.iter().any(|p| (p.index as f64 * delta_f - freq2).abs() < freq_error)
}

fn detect_harmonics(peaks: &[Peak], fundamental_idx: usize) -> Vec<usize> {
    let mut harmonics = Vec::new();
    for peak in peaks {
        let ratio = peak.index as f64 / fundamental_idx as f64;
        if (ratio.round() - ratio).abs() < 0.1 && ratio.round() > 1.0 {
            harmonics.push(peak.index);
        }
    }
    harmonics
}

fn find_best_period(data: &[f64], candidate_period: f64, delta_t: f64, num_bins: usize) -> (f64, f64) {
    let mut best_period = candidate_period;
    let mut max_variation = 0.0;
    let steps = 20; // Number of steps to scan around the candidate period
    let min_period = candidate_period * 0.9;
    let max_period = candidate_period * 1.1;
    let step_size = (max_period - min_period) / steps as f64;

    for step in 0..=steps {
        let period = min_period + step as f64 * step_size;
        let (binned, _) = phase_binning(data, period, delta_t, num_bins);
        let variation = binned.iter().max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)).unwrap_or(&0.0) -
                       binned.iter().min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)).unwrap_or(&0.0);
        if variation > max_variation {
            max_variation = variation;
            best_period = period;
        }
    }

    let error = step_size / 2.0;

    (best_period, error)
}

pub fn plot_results(data: Data, temp_gnu_file: PathBuf, temp_data_file: PathBuf, is_synthetic: bool) -> std::io::Result<()> {
    let mut temp_gnu_writer = Box::new(File::create(&temp_gnu_file).expect("Unable to create temporary file"));
    let mut temp_data_writer = Box::new(File::create(&temp_data_file).expect("Unable to create temporary file"));
    let raw_data = data.raw_data;
    let ak = data.ak;
    let bk = data.bk;
    let power = data.power;
    let _peaks = data.peaks;
    let _harmonics = data.harmonics;
    let n = data.n;
    let delta_t = data.delta_t;
    let bins = data.bins;
    let phases_of_bins = data.phases_of_bins;

    for (i, &p) in power.iter().enumerate() {
        if i == 0 {
            continue;
        }
        let freq = i as f64 / (n as f64 * delta_t);
        writeln!(temp_data_writer, "{} {}", freq, p)?;
    }

    writeln!(temp_gnu_writer, "set terminal tikz tex")?;
    if is_synthetic {
        writeln!(temp_gnu_writer, "set output 'Pfp-latex/plots/ft_synthetic.tex'")?;
        writeln!(temp_gnu_writer, "set xtics add ('$k_1$' 23)")?;
        writeln!(temp_gnu_writer, "set xtics add ('$k_2$' 67)")?;
        writeln!(temp_gnu_writer, "set xtics out")?;
    } else {
        writeln!(temp_gnu_writer, "set output 'Pfp-latex/plots/ft.tex'")?;
        writeln!(temp_gnu_writer, "set xtics add ('$k_1$' 30.273)")?;
        writeln!(temp_gnu_writer, "set xtics add ('$k_2$' 59.570)")?;
        writeln!(temp_gnu_writer, "set xtics add ('$k_3$' 89.844)")?;
        writeln!(temp_gnu_writer, "set x2tics ('$k_4$' 99.609)")?;
        writeln!(temp_gnu_writer, "set xtics add ('$k_5$' 121.094)")?;
    }
    writeln!(temp_gnu_writer, "set xtics out")?;
    writeln!(temp_gnu_writer, "set xlabel 'Frequency (Hz)'")?;
    writeln!(temp_gnu_writer, "set mxtics 5")?;
    writeln!(temp_gnu_writer, "set mytics 5")?;
    writeln!(temp_gnu_writer, "set ylabel 'Power'")?;
    writeln!(temp_gnu_writer, "set arrow from {}, graph 0 to {}, graph 1 nohead dt '.'", (n-5) / 2, (n-5) / 2)?;
    writeln!(temp_gnu_writer, "plot '{}' using 1:2 notitle with lines lt black", temp_data_file.to_str().unwrap())?;
    temp_data_writer.flush()?;
    temp_gnu_writer.flush()?;
    run_gnuplot(&temp_gnu_file)?;

    let mut temp_gnu_writer = Box::new(File::create(&temp_gnu_file).expect("Unable to create temporary file"));
    let mut temp_data_writer = Box::new(File::create(&temp_data_file).expect("Unable to create temporary file"));
    let x_t = {
        (0..n).into_iter().map(|k| {
            let theta = 2.0 * PI * k as f64 / n as f64;
            let t = k as f64 * delta_t;
            ak[k] * (theta * t).cos() + bk[k] * (theta * t).sin()
        }).collect::<Vec<f64>>()
    };
    for (i, (x, dat)) in x_t.iter().zip(raw_data.iter()).enumerate() {
        if i == 0 {
            continue;
        }
        writeln!(temp_data_writer, "{} {} {}", i as f64 * delta_t, x, dat)?;
    }
    writeln!(temp_gnu_writer, "set terminal tikz tex")?;
    if is_synthetic {
        writeln!(temp_gnu_writer, "set output 'Pfp-latex/plots/reconstruction_synthetic.tex'")?;
    } else {
        writeln!(temp_gnu_writer, "set output 'Pfp-latex/plots/reconstruction.tex'")?;
    }
    writeln!(temp_gnu_writer, "set mxtics 5")?;
    writeln!(temp_gnu_writer, "set mytics 5")?;
    writeln!(temp_gnu_writer, "set xlabel 'Time (s)'")?;
    writeln!(temp_gnu_writer, "set ylabel 'Amplitude'")?;
    writeln!(temp_gnu_writer, "set xtics out")?;
    writeln!(temp_gnu_writer, "plot [0:1] [] '{}' using 1:2 notitle with lines lt black", temp_data_file.to_str().unwrap())?;

    writeln!(temp_gnu_writer, "set output 'Pfp-latex/plots/reconstruction_raw.tex'")?;
    writeln!(temp_gnu_writer, "plot [0:1] [] '{}' using 1:3 notitle with lines lt black", temp_data_file.to_str().unwrap())?;


    // writeln!(temp_gnu_writer, "set multiplot")?;
    // writeln!(temp_gnu_writer, "set size 1.0, 0.5")?;
    // writeln!(temp_gnu_writer, "set origin 0.0, 0.0")?;
    // writeln!(temp_gnu_writer, "set mxtics 5")?;
    // writeln!(temp_gnu_writer, "set mytics 5")?;
    // writeln!(temp_gnu_writer, "set xlabel 'Time (s)'")?;
    // writeln!(temp_gnu_writer, "set ylabel 'Amplitude'")?;
    // writeln!(temp_gnu_writer, "plot [0:1] [] '{}' using 1:2 notitle with lines lt black", temp_data_file.to_str().unwrap())?;
    // writeln!(temp_gnu_writer, "set origin 0.0, 0.5")?;
    // writeln!(temp_gnu_writer, "plot [0:1] [] '{}' using 1:3 notitle with lines lt black", temp_data_file.to_str().unwrap())?;
    // writeln!(temp_gnu_writer, "unset multiplot")?;

    temp_data_writer.flush()?;
    temp_gnu_writer.flush()?;
    run_gnuplot(&temp_gnu_file)?;

    // let mut temp_gnu_writer = Box::new(File::create(&temp_gnu_file).expect("Unable to create temporary file"));
    // let mut temp_data_writer = Box::new(File::create(&temp_data_file).expect("Unable to create temporary file"));
    // // find correspoteding ak and bk for the peaks
    // let mut ak_peaks = Vec::new();
    // let mut bk_peaks = Vec::new();
    // for peak in &peaks {
    //     let k = peak.index;
    //     ak_peaks.push(ak[k]);
    //     bk_peaks.push(bk[k]);
    // }
    // writeln!(temp_gnu_writer, "set terminal png")?;
    // writeln!(temp_gnu_writer, "set output 'results2.png'")?;
    // writeln!(temp_gnu_writer, "set samples 1000")?;
    // writeln!(temp_gnu_writer, "plot [t = 0:{}] {}*cos(2*pi*{}*t/{}) + {}*sin(2*pi*{}*t/{}) notitle with lines",
    //      n as f64, 
    //      ak_peaks[0], 
    //      ak.iter().position(|n| n == &ak_peaks[0]).unwrap() as f64, 
    //      n as f64,
    //      bk_peaks[0], 
    //      bk.iter().position(|n| n == &bk_peaks[0]).unwrap() as f64, 
    //      n as f64,
    // )?;
    // temp_data_writer.flush()?;
    // temp_gnu_writer.flush()?;
    // run_gnuplot(&temp_gnu_file)?;

    // let mut temp_gnu_writer = Box::new(File::create(&temp_gnu_file).expect("Unable to create temporary file"));
    // // let mut temp_data_writer = Box::new(File::create(&temp_data_file).expect("Unable to create temporary file"));
    // let dominant_peak = &peaks[0];
    // let k = dominant_peak.index;
    // let a_k = ak[k];
    // let b_k = bk[k];
    // let freq_hz = k as f64 / (n as f64 * delta_t); // Frequency in Hz
    // writeln!(temp_gnu_writer, "set terminal png")?;
    // writeln!(temp_gnu_writer, "set output 'Pfp-latex/plots/results2.png'")?;
    // writeln!(temp_gnu_writer, "set samples 1000")?;
    // writeln!(temp_gnu_writer, "set xlabel 'Time (s)'")?;
    // writeln!(temp_gnu_writer, "set ylabel 'Amplitude'")?;
    // writeln!(temp_gnu_writer, "set grid")?;
    // writeln!(
    //     temp_gnu_writer,
    //     "plot [t=0:{}] {}*cos(2*pi*{}*t) + {}*sin(2*pi*{}*t) title 'Dominant Frequency: {:.2}' with lines lw 2",
    //     (n as f64) * delta_t, // Total time in seconds
    //     a_k,
    //     freq_hz, // Use physical frequency (Hz)
    //     b_k,
    //     freq_hz,
    //     freq_hz,
    // )?;
    // run_gnuplot(&temp_gnu_file)?;
    
    let mut temp_gnu_writer = Box::new(File::create(&temp_gnu_file).expect("Unable to create temporary file"));
    let mut temp_data_writer = Box::new(File::create(&temp_data_file).expect("Unable to create temporary file"));
    // create a histogram of the phase bins
    //
    for (phase, &count) in phases_of_bins.iter().zip(bins.iter()) {
        let stddev = count.sqrt();
        writeln!(temp_data_writer, "{} {} {}", phase/(2.0*PI), count, stddev)?;
    }
    for (phase, &count) in phases_of_bins.iter().zip(bins.iter()) {
        let stddev = count.sqrt();
        writeln!(temp_data_writer, "{} {} {}", phase/(2.0*PI) + 1.0, count, stddev)?;
    }
    writeln!(temp_gnu_writer, "set terminal tikz tex")?;
    if is_synthetic {
        writeln!(temp_gnu_writer, "set output 'Pfp-latex/plots/phasebins_synthetic.tex'")?;
    } else {
        writeln!(temp_gnu_writer, "set output 'Pfp-latex/plots/phasebins.tex'")?;
    }
    writeln!(temp_gnu_writer, "set xlabel 'Pulse Phase'")?;
    writeln!(temp_gnu_writer, "set ylabel 'Count'")?;
    writeln!(temp_gnu_writer, "set mxtics 5")?;
    writeln!(temp_gnu_writer, "set mytics 5")?;
    writeln!(temp_gnu_writer, "set xtics out")?;
    writeln!(temp_gnu_writer, "set style fill empty border 0")?;
    writeln!(temp_gnu_writer, "plot [] [] '{}' using ($1+0.045):2:3 notitle with yerrorbars lt black, '{}' using 1:2 notitle with hsteps forward lt black", temp_data_file.to_str().unwrap(), temp_data_file.to_str().unwrap())?;
    temp_data_writer.flush()?;
    temp_gnu_writer.flush()?;
    run_gnuplot(&temp_gnu_file)?;

    Ok(())
}

fn estimate_noise_floor(power_spectrum: &[f64]) -> f64 {
    let mut sorted = power_spectrum.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    sorted[sorted.len() / 2]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn synthetic_data() {
        let freq1 = 21.0;
        let freq2 = 67.0;
        let synthetic_data = generate_synthetic_data(freq1, freq2);
        let data = analyze_pulsar_data(synthetic_data, Box::new(std::io::sink()), Some((freq1, freq2)));
        let delta_f = 1.0 / (data.n as f64 * data.delta_t);
        let freq_error = delta_f / 2.0;
        assert!(validate_synthetic_data(freq1, freq2, &data.peaks, freq_error, delta_f));
    }
}
