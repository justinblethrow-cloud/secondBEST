# K-mer Methods Summary

This is a compact methods note for BEST k-mer context reports and figures.

## K-mer Intervals

BEST can summarize errors over interval sets. In k-mer mode, each interval is a
reference-aligned occurrence of a read-oriented k-mer. K-mer intervals overlap,
so their rates are context summaries rather than independent base-level counts.
Use them to identify and compare sequence contexts, not as independent
hypothesis-test units.

## Error Rates and Pseudo-QV

Raw k-mer error rate is the empirical error fraction observed for intervals
assigned to that k-mer or reverse-complement-collapsed context.

Platform comparison figures transform error rate to pseudo-QV:

```text
pseudo-QV = -10 * log10(error_rate)
```

Higher pseudo-QV means lower empirical error. The transform is useful for
visualizing orders-of-magnitude differences, but it is not a recalibrated base
quality model. A configurable `--qv-error-floor` prevents zero-error contexts
from becoming infinite; the default `1e-6` caps zero-rate contexts at Q60.

## Reverse-Complement Handling

BEST reports raw k-mers in read orientation because that is useful for detecting
strand, alignment, and platform artifacts. Some chemistry-driven contexts are
better interpreted after grouping a k-mer with its reverse complement.

Use raw k-mer plots to diagnose orientation-specific behavior. Use
reverse-complement-collapsed plots when the scientific claim should be
independent of read orientation.

## Ranking and Support

`risk_ranked_kmers.csv` prioritizes k-mers using support-weighted error behavior
with emphasis on non-homopolymer indels and mismatches. Figure panels show
interval support where possible because high-error, low-support contexts are
less reliable than similarly high-error contexts with many intervals.

Typical presentation thresholds should include:

- `--min-intervals`: minimum k-mer support included in figures.
- `--top-n`: number of ranked contexts shown.
- Whether reverse-complement-collapsed or raw k-mers are being discussed.

## Clustered Heatmap

The top-k-mer clustered heatmap clusters high-risk raw k-mers using standardized
error-profile features plus one-hot sequence features. The dendrogram is an
exploratory grouping of sequence/error-profile similarity. It is not a
biological tree, a phylogeny, or a formal model of platform chemistry.

SciPy is optional. If SciPy is installed, the figure includes a hierarchical row
dendrogram. Without SciPy, BEST still writes the heatmap using ranked order.

## Plotting Dependencies

Plotting scripts require `matplotlib` and `numpy`. `scipy` is optional at
runtime but recommended for dendrogram support. The plotting requirements file
installs the full recommended plotting stack:

```bash
python3 -m pip install -r requirements-plotting.txt
```

## Provenance

Figure manifests record the command, selected plotting options, local git
commit when available, output formats, dependency versions, and whether SciPy
clustering was available. Keep the manifest with exported figures so slides can
be traced back to the exact filters and plotting transform used.
