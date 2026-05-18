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

## Figure Command

After generating the exploration tables, use `scripts/kmer_explore_figures.py`
to make slide-ready figures from the output directory:

```bash
python3 scripts/kmer_explore_figures.py \
    --explore-dir target/tmp/kmer_explore \
    --out-dir target/tmp/kmer_explore/figures \
    --label "ONT simplex chr20 5x k=7" \
    --min-intervals 1000 \
    --top-n 20
```

The figure script writes PNG and SVG versions by default, plus
`figure_manifest.md` with highlights, captions, and presentation-builder notes.
Use the PNGs for dense scatter plots; point clouds are rasterized in SVG output
to keep files small enough for slide tools.

## Outputs

- `kmer_enrichment.csv`: observed errors vs expected errors under the
  k-length-specific background rate. Useful for separating common high-count
  k-mers from unusually error-enriched contexts.
- `rc_collapsed_kmers.csv`: reverse-complement-collapsed k-mer enrichment.
  Each row groups a canonical k-mer with its reverse complement, preserving the
  member list, dominant member, and member error-rate range. Use this table
  when interpreting context effects that should be independent of reported
  read orientation.
- `rc_pair_diagnostics.csv`: direct comparisons between each canonical k-mer
  and its reverse complement. This includes member-specific intervals,
  error-rate deltas, log2 error-rate ratios, and a support-weighted delta for
  ranking high-confidence orientation differences.
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
- `rc_strand_mirror.csv`: reverse-complement pair diagnostics split by forward
  and reverse alignment strand. This helps distinguish a true
  reverse-complement member difference from an ordinary strand effect.
- `substitution_spectrum.csv`: offset-aware substitution spectra such as `G>A`
  or `C>T`, available when `--kmer-advanced-stats` was used.
- `kmer_exploration_report.md`: compact Markdown report with top tables.

## Figure Outputs

- `nanopore_context_error_classes`: top reverse-complement-collapsed contexts
  with stacked mismatch, non-homopolymer indel, and homopolymer indel rates.
- `nanopore_homopolymer_run_profile`: error rate by homopolymer run length and
  base, intended to show whether the expected Nanopore homopolymer signal is
  visible.
- `nanopore_rc_pair_asymmetry`: scatter plot comparing each canonical k-mer
  with its reverse complement.
- `nanopore_rc_strand_mirror`: support-ranked mirror diagnostic that helps
  distinguish reverse-complement member differences from ordinary strand bias.
- `nanopore_quality_calibration`: predicted mean QV vs empirical context QV.
- `nanopore_motif_enrichment`: model-free motif discovery among high-error
  k-mers.
- `nanopore_substitution_spectrum`: substitution ranking derived from the
  offset-aware substitution table.
- `nanopore_offset_base_error_heatmap`: compact base-by-offset summary for
  explaining how the k-mer summaries can be decomposed into positional effects.

## Interpretation Notes

K-mer intervals overlap, so these reports are context summaries rather than
partitioned base counts. Enrichment and risk scores should be used to prioritize
contexts for inspection, not as formal hypothesis tests. BEST reports k-mers in
read orientation, which is useful for diagnosing platform and alignment
behavior, but it can split a chemistry-driven context across a k-mer and its
reverse complement. For public claims, prefer contexts that have high support,
stable replicate behavior, similar signal in the reverse-complement-collapsed
table, and a clear mechanistic pattern in the position, homopolymer-phase,
strand-mirror, or substitution-spectrum tables.

For presentation use, anchor the story around error-context discovery rather
than any single chr20 pilot result. A conservative slide sequence is:
context/error-class enrichment, homopolymer run-length profile,
reverse-complement pair asymmetry plus strand mirror, quality calibration,
motif discovery, and substitution spectrum as supporting detail.
