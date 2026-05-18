# K-mer Exploration Reports

`scripts/kmer_explore.py` is a downstream reporting tool for BEST k-mer output
prefixes. It does not read BAMs; it consumes the CSVs produced by
`--intervals-kmer` and writes derived tables that are easier to use for
scientific interpretation.

## Recommended BEST Run

For a context-rich run:

```bash
best -t 16 \
    --no-per-aln-stats \
    --no-feature-qual-score-stats \
    --kmer-position-stats \
    --kmer-advanced-stats \
    --intervals-kmer 7 \
    -- input.bam reference.fasta output
```

`--kmer-position-stats` enables per-offset and homopolymer-phase exploration.
`--kmer-advanced-stats` enables strand asymmetry and substitution spectrum
exploration. Leave either flag off when the extra tables are not needed.

## Report Command

```bash
python3 scripts/kmer_explore.py \
    --prefix output \
    --label ont_duplex_v1 \
    --out-dir target/tmp/kmer_explore \
    --min-intervals 1000 \
    --top-n 50
```

Provide multiple `--prefix` and `--label` pairs to add replicate reliability
summaries:

```bash
python3 scripts/kmer_explore.py \
    --prefix run1/output --label replicate_1 \
    --prefix run2/output --label replicate_2 \
    --out-dir target/tmp/kmer_replicates
```

## Outputs

- `kmer_enrichment.csv`: observed errors vs expected errors under the
  k-length-specific background rate. Useful for separating common high-count
  k-mers from unusually error-enriched contexts.
- `quality_calibration.csv`: empirical k-mer QV vs mean predicted quality.
  Positive `overconfidence_qv` marks contexts where predicted quality is higher
  than empirical accuracy.
- `hp_phase_profile.csv`: error rates grouped by base, homopolymer run length
  within the k-mer, and phase such as singleton, boundary, interior, or edge.
- `position_heatmap.csv`: heatmap-ready offset/base rates for mismatch and indel
  classes.
- `motif_enrichment.csv`: model-free short motif enrichment among high-error
  k-mers, including position-specific motifs.
- `risk_ranked_kmers.csv`: contexts ranked by a support-weighted error score
  that emphasizes non-homopolymer indels and mismatches.
- `replicate_reliability.csv`: mean, variance, and spread of k-mer error rates
  across multiple prefixes.
- `strand_asymmetry.csv`: forward vs reverse error-rate differences, available
  when `--kmer-advanced-stats` was used.
- `substitution_spectrum.csv`: offset-aware substitution spectra such as `G>A`
  or `C>T`, available when `--kmer-advanced-stats` was used.
- `kmer_exploration_report.md`: compact Markdown report with top tables.

## Interpretation Notes

K-mer intervals overlap, so these reports are context summaries rather than
partitioned base counts. Enrichment and risk scores should be used to prioritize
contexts for inspection, not as formal hypothesis tests. For public claims,
prefer contexts that have high support, stable replicate behavior, and a clear
mechanistic pattern in the position, homopolymer-phase, or substitution-spectrum
tables.
