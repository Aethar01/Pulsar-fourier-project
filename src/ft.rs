use rustfft::{FftPlanner, num_complex::Complex};
use std::path::PathBuf;
use std::io::Write;
use std::f64::consts::PI;


pub fn fft(input: PathBuf, output: Option<PathBuf>) {
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


pub fn compute_fourier_transform(data: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
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

