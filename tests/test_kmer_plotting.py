from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

ROOT = Path(__file__).resolve().parents[1]


def require_plotting_stack() -> None:
    try:
        import matplotlib  # noqa: F401
        import numpy  # noqa: F401
    except ImportError as exc:
        raise unittest.SkipTest(f"plotting dependencies unavailable: {exc}") from exc


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def run_script(args: list[str], tmpdir: Path) -> None:
    env = os.environ.copy()
    env["MPLCONFIGDIR"] = str(tmpdir / "mplconfig")
    subprocess.run(
        [sys.executable, "-B", *args],
        cwd=ROOT,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )


def make_explore_fixture(explore_dir: Path) -> None:
    write_csv(
        explore_dir / "risk_ranked_kmers.csv",
        [
            {
                "sample": "toy",
                "kmer_len": 7,
                "kmer": "ACGACGA",
                "intervals": 1000,
                "identity": 0.94,
                "error_rate": 0.060,
                "errors_per_interval": 0.42,
                "mismatch_fraction": 0.20,
                "indel_fraction": 0.80,
                "non_hp_indel_fraction": 0.70,
                "risk_score": 0.90,
            },
            {
                "sample": "toy",
                "kmer_len": 7,
                "kmer": "CGACGAC",
                "intervals": 900,
                "identity": 0.95,
                "error_rate": 0.050,
                "errors_per_interval": 0.35,
                "mismatch_fraction": 0.55,
                "indel_fraction": 0.45,
                "non_hp_indel_fraction": 0.30,
                "risk_score": 0.75,
            },
            {
                "sample": "toy",
                "kmer_len": 7,
                "kmer": "TTTTTTT",
                "intervals": 800,
                "identity": 0.93,
                "error_rate": 0.070,
                "errors_per_interval": 0.49,
                "mismatch_fraction": 0.10,
                "indel_fraction": 0.90,
                "non_hp_indel_fraction": 0.10,
                "risk_score": 0.70,
            },
        ],
    )
    write_csv(
        explore_dir / "rc_collapsed_kmers.csv",
        [
            {
                "sample": "toy",
                "kmer_len": 7,
                "canonical_kmer": "AAAAAAA",
                "members": "AAAAAAA;TTTTTTT",
                "member_count": 2,
                "intervals": 1800,
                "passes_min_intervals": "true",
                "identity": 0.93,
                "error_rate": 0.070,
                "errors_per_interval": 0.49,
                "observed_errors": 882,
                "expected_errors": 440,
                "enrichment": 2.00,
                "z_score": 15.0,
                "mismatches_per_interval": 0.020,
                "non_hp_ins_per_interval": 0.030,
                "non_hp_del_per_interval": 0.040,
                "hp_ins_per_interval": 0.050,
                "hp_del_per_interval": 0.060,
                "dominant_member": "TTTTTTT",
                "member_error_rate_range": 0.010,
            },
            {
                "sample": "toy",
                "kmer_len": 7,
                "canonical_kmer": "ACGACGA",
                "members": "ACGACGA;TCGTCGT",
                "member_count": 2,
                "intervals": 1900,
                "passes_min_intervals": "true",
                "identity": 0.94,
                "error_rate": 0.060,
                "errors_per_interval": 0.42,
                "observed_errors": 798,
                "expected_errors": 500,
                "enrichment": 1.60,
                "z_score": 10.0,
                "mismatches_per_interval": 0.030,
                "non_hp_ins_per_interval": 0.020,
                "non_hp_del_per_interval": 0.010,
                "hp_ins_per_interval": 0.020,
                "hp_del_per_interval": 0.010,
                "dominant_member": "ACGACGA",
                "member_error_rate_range": 0.004,
            },
        ],
    )
    write_csv(
        explore_dir / "rc_pair_diagnostics.csv",
        [
            {
                "sample": "toy",
                "kmer_len": 7,
                "canonical_kmer": "ACGACGA",
                "member_a": "ACGACGA",
                "member_b": "TCGTCGT",
                "pair_status": "complete",
                "member_a_intervals": 1000,
                "member_b_intervals": 900,
                "min_member_intervals": 900,
                "member_a_error_rate": 0.060,
                "member_b_error_rate": 0.050,
                "delta_b_minus_a": -0.010,
                "abs_delta_error_rate": 0.010,
                "log2_error_rate_ratio": -0.25,
                "support_weighted_delta": 9.0,
            }
        ],
    )
    write_csv(
        explore_dir / "rc_strand_mirror.csv",
        [
            {
                "sample": "toy",
                "kmer_len": 7,
                "canonical_kmer": "ACGACGA",
                "member_a": "ACGACGA",
                "member_b": "TCGTCGT",
                "pair_status": "complete",
                "min_member_strand_intervals": 400,
                "member_a_forward_error_rate": 0.061,
                "member_a_reverse_error_rate": 0.059,
                "member_b_forward_error_rate": 0.050,
                "member_b_reverse_error_rate": 0.049,
                "mirror_discordance": 0.010,
                "max_within_member_delta": 0.002,
            }
        ],
    )
    write_csv(
        explore_dir / "quality_calibration.csv",
        [
            {
                "sample": "toy",
                "kmer_len": 7,
                "kmer": "ACGACGA",
                "intervals": 1000,
                "mean_qual": 28.0,
                "empirical_qv": 12.2,
                "overconfidence_qv": 15.8,
                "error_rate": 0.060,
            }
        ],
    )
    write_csv(
        explore_dir / "hp_phase_profile.csv",
        [
            {
                "sample": "toy",
                "kmer_len": 7,
                "base": base,
                "run_length_within_kmer": run_length,
                "phase": "interior",
                "matches": 100,
                "mismatches": 2,
                "non_hp_ins": 1,
                "non_hp_del": 1,
                "hp_ins": run_length,
                "hp_del": run_length,
                "error_rate": 0.05,
            }
            for base in ["A", "C", "G", "T"]
            for run_length in [1, 2]
        ],
    )
    write_csv(
        explore_dir / "position_heatmap.csv",
        [
            {
                "sample": "toy",
                "kmer_len": 7,
                "offset": offset,
                "kmer_base": base,
                "intervals": 1000,
                "matches": 950,
                "mismatches": 10,
                "non_hp_ins": 10,
                "non_hp_del": 10,
                "hp_ins": 10,
                "hp_del": 10,
                "error_rate": 0.05 + 0.001 * offset,
            }
            for offset in range(3)
            for base in ["A", "C", "G", "T"]
        ],
    )
    write_csv(
        explore_dir / "motif_enrichment.csv",
        [
            {
                "sample": "toy",
                "motif_type": "anywhere",
                "motif": "CG",
                "motif_length": 2,
                "bad_quantile_threshold": 0.05,
                "bad_with_motif": 2,
                "kmers_with_motif": 3,
                "bad_fraction_with_motif": 0.66,
                "bad_fraction_background": 0.10,
                "odds_ratio": 8.0,
            }
        ],
    )
    write_csv(
        explore_dir / "substitution_spectrum.csv",
        [
            {
                "sample": "toy",
                "kmer_len": 7,
                "offset": 0,
                "substitution": substitution,
                "count": 5,
                "intervals": 1000,
                "per_million": 5000,
            }
            for substitution in ["A>C", "C>A"]
        ],
    )


def make_compare_fixture(parent: Path) -> tuple[Path, Path]:
    left = parent / "left"
    right = parent / "right"
    base_rows = [
        ("AAAAAAA", 1000, 0.050, 1200, 0.005),
        ("ACGACGA", 900, 0.060, 1000, 0.006),
        ("CGACGAC", 800, 0.040, 900, 0.004),
    ]
    for out_dir, side in [(left, "left"), (right, "right")]:
        write_csv(
            out_dir / "kmer_enrichment.csv",
            [
                {
                    "sample": side,
                    "kmer_len": 7,
                    "kmer": kmer,
                    "intervals": left_intervals if side == "left" else right_intervals,
                    "observed_errors": 1,
                    "expected_errors": 1,
                    "enrichment": 1,
                    "z_score": 0,
                    "identity": 1 - (left_error if side == "left" else right_error),
                    "error_rate": left_error if side == "left" else right_error,
                }
                for kmer, left_intervals, left_error, right_intervals, right_error in base_rows
            ],
        )
        write_csv(
            out_dir / "rc_collapsed_kmers.csv",
            [
                {
                    "sample": side,
                    "kmer_len": 7,
                    "canonical_kmer": kmer,
                    "members": f"{kmer};{kmer[::-1]}",
                    "member_count": 2,
                    "intervals": left_intervals if side == "left" else right_intervals,
                    "passes_min_intervals": "true",
                    "identity": 1 - (left_error if side == "left" else right_error),
                    "error_rate": left_error if side == "left" else right_error,
                    "errors_per_interval": left_error * 7 if side == "left" else right_error * 7,
                    "observed_errors": 1,
                    "expected_errors": 1,
                    "enrichment": 1,
                    "z_score": 0,
                    "mismatches_per_interval": 0.01,
                    "non_hp_ins_per_interval": 0.01,
                    "non_hp_del_per_interval": 0.01,
                    "hp_ins_per_interval": 0.01,
                    "hp_del_per_interval": 0.01,
                    "dominant_member": kmer,
                    "member_error_rate_range": 0.001,
                }
                for kmer, left_intervals, left_error, right_intervals, right_error in base_rows
            ],
        )
    return left, right


class KmerPlottingSmokeTests(unittest.TestCase):
    def test_explore_figures_generate_manifest_and_top_kmer_panels(self) -> None:
        require_plotting_stack()
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            explore_dir = tmpdir / "explore"
            out_dir = tmpdir / "figures"
            make_explore_fixture(explore_dir)
            run_script(
                [
                    "scripts/kmer_explore_figures.py",
                    "--explore-dir",
                    str(explore_dir),
                    "--out-dir",
                    str(out_dir),
                    "--label",
                    "toy",
                    "--platform",
                    "Toy",
                    "--figure-prefix",
                    "toy",
                    "--min-intervals",
                    "10",
                    "--top-n",
                    "3",
                    "--formats",
                    "png",
                ],
                tmpdir,
            )
            self.assertTrue((out_dir / "toy_top_kmer_lollipop.png").exists())
            self.assertTrue((out_dir / "toy_top_kmer_clustered_heatmap.png").exists())
            self.assertTrue((out_dir / "toy_motif_aligned_contexts.png").exists())
            manifest = (out_dir / "figure_manifest.md").read_text()
            self.assertIn("## Provenance", manifest)
            self.assertIn("SciPy clustering available", manifest)

    def test_platform_compare_generates_pseudo_qv_outputs(self) -> None:
        require_plotting_stack()
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            left, right = make_compare_fixture(tmpdir)
            out_dir = tmpdir / "compare"
            run_script(
                [
                    "scripts/kmer_compare_platforms.py",
                    "--left-dir",
                    str(left),
                    "--right-dir",
                    str(right),
                    "--left-label",
                    "Simplex",
                    "--right-label",
                    "Revio",
                    "--out-dir",
                    str(out_dir),
                    "--min-intervals",
                    "10",
                    "--top-n",
                    "2",
                    "--formats",
                    "png",
                ],
                tmpdir,
            )
            self.assertTrue((out_dir / "platform_kmer_pseudo_qv_hexbin.png").exists())
            self.assertTrue((out_dir / "platform_rc_collapsed_pseudo_qv_hexbin.png").exists())
            self.assertTrue((out_dir / "platform_kmer_delta_ranked.png").exists())
            self.assertFalse((out_dir / "platform_kmer_error_hexbin.png").exists())
            manifest = (out_dir / "platform_comparison_manifest.md").read_text()
            self.assertIn("## Provenance", manifest)
            self.assertIn("QV error floor", manifest)
            joined = (out_dir / "platform_kmer_joined.csv").read_text()
            self.assertIn("Simplex_pseudo_qv", joined)
            self.assertIn("Revio_pseudo_qv", joined)


if __name__ == "__main__":
    unittest.main()
