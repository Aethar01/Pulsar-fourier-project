use rustfft::{FftPlanner, num_complex::Complex};
use clap::Parser;
use std::path::PathBuf;
use std::fs::File;
use std::io::{self, BufRead, Write, stdout};
use std::path::Path;
use std::f64::consts::PI;



#[derive(Parser)]
#[command(name = "Pulsar Fourier Analyzer")]
#[command(about = "Analyzes pulsar data using Fourier transform and phase binning", long_about = None)]
enum Cli {
    /// Analyze pulsar data from file
    Analyze {
        /// Input data file
        #[arg(short, long, value_parser)]
        input: PathBuf,
        
        /// Output file for results (stdout if not specified)
        #[arg(short, long, value_parser)]
        output: Option<PathBuf>,
    },
    
    /// Generate and analyze synthetic data
    Synthetic {
        /// First frequency component (cosine)
        freq1: f64,
        
        /// Second frequency component (sine)
        freq2: f64,
        
        /// Output file for results (stdout if not specified)
        #[arg(short, long, value_parser)]
        output: Option<PathBuf>,
    },

    /// Compute the Fourier transform of a file and return the results to stdout or a file
    FFT {
        /// Input data file
        #[arg(short, long, value_parser)]
        input: PathBuf,

        /// Output data file
        #[arg(short, long, value_parser)]
        output: Option<PathBuf>,
    },
}

fn main() {
    let cli = Cli::parse();

    match cli {
        Cli::Analyze { input, output } => {
            let data = match read_data_file(&input) {
                Ok(data) => data,
                Err(e) => {
                    eprintln!("Error reading data file: {}", e);
                    return;
                }
            };
            
            let output_writer: Box<dyn Write> = match output {
                Some(path) => Box::new(File::create(path).expect("Unable to create output file")),
                None => Box::new(stdout()),
            };
            
            analyze_pulsar_data(&data, output_writer);
        }
        Cli::Synthetic { freq1, freq2, output } => {
            let synthetic_data = generate_synthetic_data(freq1, freq2);
            
            let output_writer: Box<dyn Write> = match output {
                Some(path) => Box::new(File::create(path).expect("Unable to create output file")),
                None => Box::new(stdout()),
            };
            
            analyze_pulsar_data(&synthetic_data, output_writer);
        }
        Cli::FFT { input, output } => {
            fft(input, output);
        }
    }
}

fn fft(input: PathBuf, output: Option<PathBuf>) {
    let input = std::fs::read_to_string(input).expect("Unable to read input file");

    let buffer = input.split("\n").into_iter()
        .map(|x| x.trim())
        .filter(|x| !x.is_empty())
        .collect::<Vec<&str>>();   
    let mut complex_buffer = Vec::new();
    buffer.iter().for_each(|x| complex_buffer.push(Complex { re: x.parse::<f32>().unwrap(), im: 0.0 }));

    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(buffer.len());

    fft.process(&mut complex_buffer);

    // complex_buffer.remove(0);

    let buff_string = complex_buffer.iter().map(|complex| format!("{}\n", complex)).collect::<String>();

    match output {
        Some(path) => {
            std::fs::write(&path, buff_string).unwrap();
            println!("Wrote to {}", path.display());
    }
        None => {
            let mut output = std::io::stdout();
            output.write_all(buff_string.as_bytes()).unwrap();
        }
    }
}

fn read_data_file<P>(filename: P) -> io::Result<Vec<f64>>
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

fn generate_synthetic_data(freq1: f64, freq2: f64) -> Vec<f64> {
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

fn analyze_pulsar_data<W: Write>(data: &[f64], mut output: W) {
    let n = data.len();
    let delta_t = 0.004; // 4ms sampling interval
    
    // Compute Fourier coefficients and power spectrum
    let (_ak, _bk, power) = compute_fourier_transform(data);
    
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

fn compute_fourier_transform(data: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = data.len();
    let mut ak = vec![0.0; n];
    let mut bk = vec![0.0; n];
    let mut power = vec![0.0; n];
    
    for k in 0..n {
        let theta = 2.0 * PI * k as f64 / n as f64;
        
        // Initialize recurrence relation
        let mut u = vec![0.0; n + 2];
        
        // Compute U_n using recurrence relation
        for i in (0..n).rev() {
            u[i] = data[i] + 2.0 * theta.cos() * u[i + 1] - u[i + 2];
        }
        
        // Compute A_k and B_k
        ak[k] = (u[0] - u[1] * theta.cos()) / n as f64;
        bk[k] = (u[1] * theta.sin()) / n as f64;
        
        // Compute power
        power[k] = ak[k].powi(2) + bk[k].powi(2);
    }
    
    (ak, bk, power)
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

