use clap::Parser;
use std::path::PathBuf;
use std::fs::File;
use std::io::{Write, stdout};

mod tools;
mod ft;

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
            let data = match tools::read_data_file(&input) {
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
            
            tools::analyze_pulsar_data(&data, output_writer, None);
        }
        Cli::Synthetic { freq1, freq2, output } => {
            let synthetic_data = tools::generate_synthetic_data(freq1, freq2);
            
            let output_writer: Box<dyn Write> = match output {
                Some(path) => Box::new(File::create(path).expect("Unable to create output file")),
                None => Box::new(stdout()),
            };
            
            tools::analyze_pulsar_data(&synthetic_data, output_writer, Some((freq1, freq2)));
        }

        Cli::FFT { input, output } => {
            ft::fft(input, output);
        }
    }
}
