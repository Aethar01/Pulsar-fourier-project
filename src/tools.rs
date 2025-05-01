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
    
    for t in 0..n {
        let t_f64 = t as f64;
        let value = 10.0 * (2.0 * PI * freq1 * t_f64 / n as f64).cos() + 
                   5.0 * (2.0 * PI * freq2 * t_f64 / n as f64).sin();
        data.push(value);
    }
    
    data
}

pub struct Data {
    ak: Vec<f64>,
    bk: Vec<f64>,
    power: Vec<f64>,
    peaks: Vec<Peak>,
    harmonics: Vec<usize>,
    n: usize,
    delta_t: f64,
}

pub fn analyze_pulsar_data<W: Write>(data: &[f64], mut output: W, synthetic: Option<(f64, f64)>) -> Data {
    let n = data.len();
    let delta_t = 0.004; // 4ms sampling interval
    
    // Compute Fourier coefficients and power spectrum
    let (ak, bk, power) = ft::compute_fourier_transform(data);
    
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
    
    // Phase binning analysis for the strongest 10 peaks
    for i in 0..10 {
        if let Some(main_peak) = peaks.get(i) {
            let freq = main_peak.index as f64 / (n as f64 * delta_t);
            let period = 1.0 / freq;
            
            writeln!(output, "\nPhase Binning Analysis for {:.3} Hz (Period: {:.3} s)", freq, period).unwrap();
            
            let binned_data = phase_binning(data, period, delta_t, 10);
            writeln!(output, "\nPhase Bins:").unwrap();
            for (i, &bin) in binned_data.iter().enumerate() {
                writeln!(output, "Bin {}: {:.3}", i, bin).unwrap();
            }
            writeln!(output, "\nSum of Phase Bins: {}", binned_data.iter().sum::<f64>()).unwrap();
        }
    }

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
                writeln!(output, "Harmonic at {:.3} Hz ({}x fundamental)", *h as f64 / (n as f64 * delta_t), h / main_peak.index).unwrap();
            }
            harmonics
        } else {
            Vec::new()
        }
    };

    writeln!(output, "\nBest Period Detection:").unwrap();
    if let Some(main_peak) = peaks.first() {
        let best_period = find_best_period(data, 1.0 / main_peak.index as f64, delta_t, 10);
        writeln!(output, "Fundamental Frequency: {:.3} Hz", main_peak.index as f64 / (n as f64 * delta_t)).unwrap();
        writeln!(output, "Best Period: {:.3} s", best_period).unwrap();
    }

    if let Some((freq1, freq2)) = synthetic {
        writeln!(output, "\nSynthetic Data Validation:").unwrap();
        writeln!(output, "Fundamental Frequency: {:.3} Hz", freq1 / (n as f64 * delta_t)).unwrap();
        writeln!(output, "Second Fundamental Frequency: {:.3} Hz", freq2 / (n as f64 * delta_t)).unwrap();
        writeln!(output, "Validation: {}", validate_synthetic_data(freq1, freq2, &peaks)).unwrap();
    }

    Data {
        ak,
        bk,
        power,
        peaks,
        harmonics,
        n,
        delta_t,
    }
}

pub struct Peak {
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
    bins.iter_mut()
        .zip(counts.iter())
        .for_each(|(bin, &count)| {
            *bin = if count > 0 { *bin / count as f64 } else { f64::NAN };
        });
    
    bins
}

pub fn validate_synthetic_data(freq1: f64, freq2: f64, peaks: &[Peak]) -> bool {
    let n = 256;
    let delta_t = 0.004;
    let expected_freq1 = freq1 / (n as f64 * delta_t);
    let expected_freq2 = freq2 / (n as f64 * delta_t);
    peaks.iter().any(|p| (p.index as f64 - expected_freq1).abs() < 1.0) &&
    peaks.iter().any(|p| (p.index as f64 - expected_freq2).abs() < 1.0)
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

fn find_best_period(data: &[f64], candidate_period: f64, delta_t: f64, num_bins: usize) -> f64 {
    let mut best_period = candidate_period;
    let mut max_variation = 0.0;
    let steps = 20; // Number of steps to scan around the candidate period
    let min_period = candidate_period * 0.9;
    let max_period = candidate_period * 1.1;
    let step_size = (max_period - min_period) / steps as f64;

    for step in 0..=steps {
        let period = min_period + step as f64 * step_size;
        let binned = phase_binning(data, period, delta_t, num_bins);
        let variation = binned.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap() - 
                       binned.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        if variation > max_variation {
            max_variation = variation;
            best_period = period;
        }
    }
    best_period
}

pub fn plot_results(data: Data, temp_gnu_file: PathBuf, temp_data_file: PathBuf) -> std::io::Result<()> {
    let mut temp_gnu_writer = Box::new(File::create(&temp_gnu_file).expect("Unable to create temporary file"));
    let mut temp_data_writer = Box::new(File::create(&temp_data_file).expect("Unable to create temporary file"));
    let ak = data.ak;
    let bk = data.bk;
    let power = data.power;
    let peaks = data.peaks;
    let _harmonics = data.harmonics;
    let n = data.n;
    let delta_t = data.delta_t;


    for (i, &p) in power.iter().enumerate() {
        if i == 0 {
            continue;
        }
        let freq = i as f64 / (n as f64 * delta_t);
        writeln!(temp_data_writer, "{} {}", freq, p)?;
    }
    writeln!(temp_gnu_writer, "set terminal png")?;
    writeln!(temp_gnu_writer, "set output 'results.png'")?;
    writeln!(temp_gnu_writer, "set xlabel 'Frequency (Hz)'")?;
    writeln!(temp_gnu_writer, "set ylabel 'Power'")?;
    writeln!(temp_gnu_writer, "plot '{}' using 1:2 notitle with lines", temp_data_file.to_str().unwrap())?;
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
    for (i, x) in x_t.iter().enumerate() {
        if i == 0 {
            continue;
        }
        writeln!(temp_data_writer, "{} {}", i as f64 * delta_t, x)?;
    }
    writeln!(temp_gnu_writer, "set terminal png")?;
    writeln!(temp_gnu_writer, "set output 'results1.png'")?;
    writeln!(temp_gnu_writer, "set xlabel 'Time (s)'")?;
    writeln!(temp_gnu_writer, "set ylabel 'Amplitude'")?;
    writeln!(temp_gnu_writer, "plot [0:1] [] '{}' using 1:2 notitle with lines", temp_data_file.to_str().unwrap())?;
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

    let mut temp_gnu_writer = Box::new(File::create(&temp_gnu_file).expect("Unable to create temporary file"));
    // let mut temp_data_writer = Box::new(File::create(&temp_data_file).expect("Unable to create temporary file"));
    let dominant_peak = &peaks[0];
    let k = dominant_peak.index;
    let a_k = ak[k];
    let b_k = bk[k];
    let freq_hz = k as f64 / (n as f64 * delta_t); // Frequency in Hz
    writeln!(temp_gnu_writer, "set terminal pngcairo enhanced font 'Arial,12' size 800,600")?;
    writeln!(temp_gnu_writer, "set output 'results2.png'")?;
    writeln!(temp_gnu_writer, "set samples 1000")?;
    writeln!(temp_gnu_writer, "set xlabel 'Time (s)'")?;
    writeln!(temp_gnu_writer, "set ylabel 'Amplitude'")?;
    writeln!(temp_gnu_writer, "set grid")?;
    writeln!(
        temp_gnu_writer,
        "plot [t=0:{}] {}*cos(2*pi*{}*t) + {}*sin(2*pi*{}*t) title 'Dominant Frequency: {:.2}' with lines lw 2",
        (n as f64) * delta_t, // Total time in seconds
        a_k,
        freq_hz, // Use physical frequency (Hz)
        b_k,
        freq_hz,
        freq_hz,
    )?;
    run_gnuplot(&temp_gnu_file)?;

    Ok(())
}
