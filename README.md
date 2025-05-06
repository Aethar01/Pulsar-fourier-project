# Pulsar-fourier-project

## Dependencies
Build Dependencies:
- Rust
- Cargo

Runtime Dependencies:
- GNUPlot (for plotting, unneeded if the `--plot` flag is not used)

## Build
```bash
cargo build --release
```

## Usage
```
Usage: Pulsar-fourier-project <COMMAND>

Commands:
  analyze    Analyze pulsar data from file
  synthetic  Generate and analyze synthetic data
  fft        Compute the Fourier transform of a file and return the results to stdout or a file
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help
```

### Analyze
```
Analyze pulsar data from file

Usage: Pulsar-fourier-project analyze [OPTIONS] --input <INPUT>

Options:
  -i, --input <INPUT>    Input data file
  -o, --output <OUTPUT>  Output file for results (stdout if not specified)
  -p, --plot             Plot the results
  -h, --help             Print help
```

### Synthetic
```
Generate and analyze synthetic data

Usage: Pulsar-fourier-project synthetic [OPTIONS] <FREQ1> <FREQ2>

Arguments:
  <FREQ1>  First frequency component (cosine)
  <FREQ2>  Second frequency component (sine)

Options:
  -o, --output <OUTPUT>  Output file for results (stdout if not specified)
  -p, --plot             Plot the results
  -h, --help             Print help
```

### FFT
```
Compute the Fourier transform of a file and return the results to stdout or a file

Usage: Pulsar-fourier-project fft [OPTIONS] --input <INPUT>

Options:
  -i, --input <INPUT>    Input data file
  -o, --output <OUTPUT>  Output data file
  -h, --help             Print help
```
