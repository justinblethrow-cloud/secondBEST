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
    --platform Nanopore \
    --figure-prefix nanopore \
    --min-intervals 1000 \
    --top-n 20
```

The figure script writes PNG and SVG versions by default, plus
`figure_manifest.md` with highlights, captions, and presentation-builder notes.
Use the PNGs for dense scatter plots; point clouds are rasterized in SVG output
to keep files small enough for slide tools. The default platform label and
figure prefix are Nanopore-oriented; set `--platform "PacBio Revio"` and
`--figure-prefix revio` for Revio/PacBio output.

## Platform Comparison Command

After running `scripts/kmer_explore.py` for two datasets, use
`scripts/kmer_compare_platforms.py` to make direct k-mer-by-k-mer comparison
plots:

```bash
python3 scripts/kmer_compare_platforms.py \
    --left-dir target/tmp/kmer_simplex/explore \
    --right-dir target/tmp/kmer_revio/explore \
    --left-label Simplex \
    --right-label Revio \
    --out-dir target/tmp/kmer_platform_compare \
    --min-intervals 1000 \
    --top-n 50
```

The comparison script joins shared k-mers that pass the support threshold in
both datasets. It writes one hexbin for raw k-mers and one for
reverse-complement-collapsed contexts, each with pseudo-QV axes computed as
`-10 * log10(error_rate)`. The default `--qv-error-floor 1e-6` caps zero-rate
contexts at Q60. The equality diagonal marks equal pseudo-QV; points above it
have higher pseudo-QV and lower error in the right dataset, while points below
it have higher pseudo-QV and lower error in the left dataset. Joined and outlier
CSVs are written alongside the figures so selected contexts can be annotated in
a slide or manuscript.

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

The examples below use the default `nanopore` figure prefix. If
`--figure-prefix revio` is supplied, the same figures are written with a
`revio_` prefix.

- `<prefix>_context_error_classes`: top reverse-complement-collapsed contexts
  with stacked mismatch, non-homopolymer indel, and homopolymer indel rates.
- `<prefix>_top_kmer_lollipop`: top high-risk raw k-mers ranked by error rate
  and colored by dominant mismatch/indel behavior.
- `<prefix>_top_kmer_clustered_heatmap`: top high-risk raw k-mers clustered by
  sequence and error-profile features, with a row dendrogram when SciPy is
  available.
- `<prefix>_top_kmer_sequence_heatmap`: base-by-position heatmap for top raw
  k-mers with adjacent error-rate bars.
- `<prefix>_motif_aligned_contexts`: top high-risk k-mers containing the
  leading enriched motif, aligned to the motif start.
- `<prefix>_homopolymer_run_profile`: error rate by homopolymer run length and
  base, intended to show whether the expected Nanopore homopolymer signal is
  visible.
- `<prefix>_rc_pair_asymmetry`: scatter plot comparing each canonical k-mer
  with its reverse complement.
- `<prefix>_rc_strand_mirror`: support-ranked mirror diagnostic that helps
  distinguish reverse-complement member differences from ordinary strand bias.
- `<prefix>_quality_calibration`: predicted mean QV vs empirical context QV.
- `<prefix>_motif_enrichment`: model-free motif discovery among high-error
  k-mers.
- `<prefix>_substitution_spectrum`: substitution ranking derived from the
  offset-aware substitution table.
- `<prefix>_offset_base_error_heatmap`: compact base-by-offset summary for
  explaining how the k-mer summaries can be decomposed into positional effects.

## Platform Comparison Outputs

- `platform_kmer_error_hexbin`: raw shared k-mer pseudo-QV hexbin.
- `platform_rc_collapsed_error_hexbin`: reverse-complement-collapsed shared
  context pseudo-QV hexbin.
- `platform_kmer_delta_ranked`: largest raw k-mer pseudo-QV deltas between the
  two datasets.
- `platform_rc_collapsed_delta_ranked`: largest reverse-complement-collapsed
  pseudo-QV deltas between the two datasets.
- `platform_kmer_joined.csv` and `platform_kmer_outliers.csv`: raw k-mer
  comparison tables with error-rate and pseudo-QV columns.
- `platform_rc_collapsed_joined.csv` and
  `platform_rc_collapsed_outliers.csv`: reverse-complement-collapsed comparison
  tables, including pseudo-QV and per-error-class event rates when available.
- `platform_comparison_manifest.md`: compact summary with figure names and
  presentation-builder notes.

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

See `docs/kmer_presentation_notes.md` for a presentation-builder-oriented
summary of the figure set, suggested narrative order, and caveats.

## Local 20X Platform Demonstration

A local HPRC HG002 20X demonstration compared Oxford Nanopore Simplex against
PacBio Revio with `k=7`, `--kmer-position-stats`, and
`--kmer-advanced-stats`, using `/mnt/datavault/Agentic/BEST/HPRC-HG002/HPRC-HG002.fasta`
as the reference. The generated outputs were kept under `target/tmp` rather
than committed to the repository.

Run commands:

```bash
best -t 32 \
    --bam-reader-threads 8 \
    --record-batch-size 64 \
    --no-per-aln-stats \
    --no-feature-qual-score-stats \
    --kmer-position-stats \
    --kmer-advanced-stats \
    --intervals-kmer 7 \
    -- simplex_20X.HPRC_HG002.bam HPRC-HG002.fasta simplex/run

best -t 96 \
    --bam-reader-threads 16 \
    --record-batch-size 64 \
    --no-per-aln-stats \
    --no-feature-qual-score-stats \
    --kmer-position-stats \
    --kmer-advanced-stats \
    --intervals-kmer 7 \
    -- revio_all.20X.hprc_hg002.bam HPRC-HG002.fasta revio/run
```

Observed local timings:

- Simplex 20X: 25:54 wall time, 3049% CPU, 6.7 GiB max RSS.
- Revio 20X: 4:18 wall time, 8786% CPU, 10.9 GiB max RSS.

The platform comparison joined all 16,384 raw 7-mers and 8,192
reverse-complement-collapsed contexts at `--min-intervals 1000`. In this local
run, every shared raw k-mer had a higher Simplex error rate than Revio, and the
median log2 Revio/Simplex error-rate ratio was `-4.320`. The generated
presentation outputs are:

- `target/tmp/kmer_simplex20x_hprc_k7_advanced_20260518/figures`
- `target/tmp/kmer_revio20x_hprc_k7_advanced_20260518/figures`
- `target/tmp/kmer_platform_compare_simplex20x_revio20x_20260518`
