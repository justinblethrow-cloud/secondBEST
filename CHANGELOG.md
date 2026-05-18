# Changelog

## Unreleased

### Added

- Add interval-local k-mer error profiling with `--intervals-kmer`, producing `summary_kmer_stats.csv`.
- Add optional per-position k-mer summaries with `--kmer-position-stats`, producing `summary_kmer_position_stats.csv`.
- Add k-mer length, sequence-context, and position-profile rollups for more compact interpretation of k-mer results.
- Add optional advanced k-mer strand and substitution spectrum summaries with `--kmer-advanced-stats`.
- Add `scripts/kmer_explore.py` for downstream k-mer enrichment, reverse-complement-collapsed summaries, reverse-complement pair diagnostics, quality calibration, homopolymer phase, motif discovery, replicate reliability, strand asymmetry, substitution spectrum, and risk-ranking reports.
- Add `scripts/kmer_explore_figures.py` for platform-labeled k-mer exploration figures and presentation-builder notes.
- Add `scripts/kmer_compare_platforms.py` for pairwise platform/chemistry k-mer error-rate hexbins and joined outlier tables.
- Add an optimized encoded k-mer path for k-mer-only interval runs, supporting k-mer sizes up to 32 and dense indexing for small k values.
- Add k-mer benchmarking support in `scripts/benchmark_kmer.py`, including repeat medians, JSON and Markdown reports, output-signature checks, record batch size matrices, and BGZF reader thread matrices.

### Changed

- Default `--bam-reader-threads` to 8 to improve BAM decompression throughput; reduce this on low-core or shared machines.
- Default `--record-batch-size` to 64 based on the downsampled BAM benchmark grid.
