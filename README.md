# best
Bam Error Stats Tool (best): analysis of error types in aligned reads.

`best` is used to assess the quality of reads after aligning them to a
reference assembly.

## Features
* Collect overall and per alignment stats
* Distribution of indel lengths
* Yield at different empirical Q-value thresholds
* Bin per read stats to easily examine the distribution of errors for certain
  types of reads
* Stats for regions specified by intervals (BED file, homopolymer regions,
  windows, matched sequences, k-mers etc.)
* Dedicated k-mer context summaries and optional position summaries for optimized k-mer runs
* Stats for quality scores vs empirical Q-values
* Optional suppression of per-feature quality score stats for large feature sets
* Multithreading for speed, with batched record processing and optional BGZF reader workers

## Usage
The [`best` Usage Guide](Usage.md) gives an overview of how to use `best`.

## Installing
1. Install [Rust](https://www.rust-lang.org/tools/install).
2. Clone this repository and navigate into the directory of this repository.
3. Run `cargo install --locked --path .`
4. Run `best input.bam reference.fasta prefix/path`

This will generate stats files with the `prefix/path` prefix.

## Development
### Running
1. Install [Rust](https://www.rust-lang.org/tools/install).
2. Clone this repository and navigate into the directory of this repository.
3. Run `cargo build --release`
4. Run `cargo run --release -- input.bam reference.fasta prefix/path` or
`target/release/best input.bam reference.fasta prefix/path`

This will generate stats files with the `prefix/path` prefix.

The built binary is located at `target/release/best`.

### Benchmarking
Use `scripts/benchmark_kmer.py` to benchmark k-mer modes across thread,
record-batch, and BAM-reader settings. The script writes raw CSV results plus
JSON and Markdown reports with repeat medians, peak memory, output-signature
checks, and candidate proof statements. The full usage guide includes an example
matrix and the current recommended high-core setting.

### Formatting
```
cargo fmt
```

### Comparing
Remember to pass the `-t 1` option to ensure that only one thread is used for
testing. Best generally tries to ensure the order of outputs is deterministic
with multiple threads, but the order of per-alignment stats is arbitrary unless
only one thread is used.

### Disclaimer

This is not an official Google product.

The code is not intended for use in any clinical settings.  It is not intended to be a medical device and is not intended for clinical use of any kind, including but not limited to diagnosis or prognosis.

No representations or warranties are made with regards to the accuracy of results generated.  User or licensee is responsible for verifying and validating accuracy when using this tool.
