# K-mer Presentation Notes

These notes are for building a presentation around BEST k-mer error-context
analysis, especially the local HPRC HG002 Nanopore Simplex 20X and PacBio Revio
20X comparison.

## Core Message

BEST can move from "this alignment has errors" to "these sequence contexts have
distinct, platform-specific error behavior." The strongest presentation claim
from the local 20X comparison is that the platform separation is visible for
every shared 7-mer, then remains visible after reverse-complement collapsing.

For the Simplex 20X vs Revio 20X run:

- Shared raw 7-mers: 16,384.
- Shared reverse-complement-collapsed contexts: 8,192.
- Raw 7-mers with higher Simplex error rate: 16,384.
- Raw 7-mers with higher Revio error rate: 0.
- Median Revio minus Simplex pseudo-QV delta: about +13.0.

Use those numbers as a dataset-specific demonstration, not as a universal claim
about all Nanopore or PacBio data.

## Suggested Slide Order

1. Platform overview: use `platform_kmer_error_hexbin`.
   This is the main proof figure. Axes are pseudo-QV, so higher is better. The
   diagonal is equal pseudo-QV. In the current local comparison, the density sits
   above the diagonal, showing higher Revio pseudo-QV across k-mers.

2. Orientation-robust platform overview: use
   `platform_rc_collapsed_error_hexbin`.
   This answers whether the comparison is an artifact of read orientation or
   reverse-complement splitting.

3. Concrete platform examples: use `platform_kmer_delta_ranked` and, if needed,
   `platform_rc_collapsed_delta_ranked`.
   These are label-friendly follow-ups after the hexbin establishes the global
   pattern. They provide sequences to discuss without overcrowding the hexbin.

4. Within-platform context classes: use `<prefix>_context_error_classes`.
   This shows which reverse-complement-collapsed contexts are enriched and which
   error classes drive them.

5. Top raw sequence contexts: use `<prefix>_top_kmer_lollipop` and
   `<prefix>_top_kmer_sequence_heatmap`.
   These make the actual high-risk k-mer strings visible.

6. Mechanistic grouping: use `<prefix>_top_kmer_clustered_heatmap`.
   The dendrogram groups top k-mers using sequence and error-profile features.
   Use this as an exploration slide, not as a formal phylogeny or statistical
   test.

7. Motif discovery: use `<prefix>_motif_enrichment` and
   `<prefix>_motif_aligned_contexts`.
   The enrichment plot finds compact motifs; the aligned panel makes the motif
   visually legible in top high-risk k-mers.

8. Homopolymer behavior: use `<prefix>_homopolymer_run_profile`.
   This is the clearest bridge to known long-read chemistry behavior.

9. Reverse-complement diagnostics: use `<prefix>_rc_pair_asymmetry` and
   `<prefix>_rc_strand_mirror`.
   These slides are best after the audience understands the basic context
   signal. They show whether reverse-complement pair differences persist beyond
   ordinary strand effects.

10. Calibration and secondary mechanisms: use `<prefix>_quality_calibration`,
    `<prefix>_substitution_spectrum`, and `<prefix>_offset_base_error_heatmap`
    as supporting slides.

## Figure-Specific Notes

- `platform_kmer_error_hexbin`: best single-slide platform comparison. Use the
  raw every-k-mer version for the main message.
- `platform_rc_collapsed_error_hexbin`: use immediately after the raw hexbin if
  the audience may ask about reverse-complement orientation.
- `platform_kmer_delta_ranked`: use to name concrete raw k-mers with the
  largest platform pseudo-QV differences.
- `platform_rc_collapsed_delta_ranked`: use when the story should avoid
  orientation-specific raw k-mer labels.
- `<prefix>_context_error_classes`: best single-slide within-platform context
  summary. It shows both top contexts and error class composition.
- `<prefix>_top_kmer_lollipop`: fast visual ranking of the worst raw k-mers.
  Color indicates whether the context is mismatch-heavy, non-homopolymer
  indel-heavy, or homopolymer/other-indel-heavy.
- `<prefix>_top_kmer_sequence_heatmap`: easiest way to show the audience the
  actual base strings and whether a motif or repeated base pattern is visible.
- `<prefix>_top_kmer_clustered_heatmap`: use to show structure among top
  contexts. Avoid saying the dendrogram is a biological tree; it is a clustering
  of sequence/error-profile features.
- `<prefix>_motif_enrichment`: use as a discovery result. It says which motifs
  are overrepresented among high-error k-mers.
- `<prefix>_motif_aligned_contexts`: use as the visual explanation for the
  leading motif enrichment result.
- `<prefix>_homopolymer_run_profile`: use for chemistry intuition, especially
  when discussing Nanopore homopolymer behavior.
- `<prefix>_rc_pair_asymmetry`: shows raw reverse-complement pair asymmetry.
- `<prefix>_rc_strand_mirror`: helps distinguish reverse-complement member
  differences from ordinary strand bias.
- `<prefix>_quality_calibration`: shows whether specific sequence contexts are
  overconfident relative to empirical error.
- `<prefix>_substitution_spectrum`: use when substitution biases matter; it is
  usually a support slide for long-read indel-heavy analyses.
- `<prefix>_offset_base_error_heatmap`: use as a bridge from whole-k-mer
  summaries to positional/base-level effects.

## Local Output Directories

The current local demonstration outputs are:

- Simplex figures:
  `target/tmp/kmer_simplex20x_hprc_k7_advanced_20260518/figures`
- Revio figures:
  `target/tmp/kmer_revio20x_hprc_k7_advanced_20260518/figures`
- Simplex vs Revio comparison:
  `target/tmp/kmer_platform_compare_simplex20x_revio20x_20260518`

Each output directory includes a Markdown manifest with figure files and
dataset-specific summary notes.

## Caveats

- K-mer intervals overlap, so rates are context summaries rather than
  independent base-level observations.
- The local Simplex/Revio comparison uses HPRC HG002 and a specific reference:
  `/mnt/datavault/Agentic/BEST/HPRC-HG002/HPRC-HG002.fasta`.
- Do not present one 20X run as a universal platform benchmark. Phrase it as a
  demonstration of the method and a concrete dataset comparison.
- Pseudo-QV is computed from empirical context error rate as
  `-10 * log10(error_rate)`. It is a plotting and interpretation transform, not
  a recalibrated base quality model.
- Reverse-complement-collapsed plots are preferable when the claim should be
  independent of read orientation.
