#!/usr/bin/env python
"""
Filter real drug–target and gene–disease edges by score thresholds.

- data/real/drug_targets.csv     -> data/real/drug_targets_filtered.csv
- data/real/gene_disease.csv     -> data/real/gene_disease_filtered.csv

Thresholds are:
    drug_targets.score    >= 5.0
    gene_disease.score    >= 0.05
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd


def main() -> None:
    base = Path("data/real")

    dt_in = base / "drug_targets.csv"
    gd_in = base / "gene_disease.csv"

    dt_out = base / "drug_targets_filtered.csv"
    gd_out = base / "gene_disease_filtered.csv"

    if not dt_in.exists():
        raise SystemExit(f"Missing {dt_in}")
    if not gd_in.exists():
        raise SystemExit(f"Missing {gd_in}")

    print(f"Loading {dt_in} ...")
    dt = pd.read_csv(dt_in)

    if "score" not in dt.columns:
        raise SystemExit("drug_targets.csv has no 'score' column.")

    # Drop NA scores, then filter
    dt = dt.dropna(subset=["score"])
    dt_filtered = dt[dt["score"] >= 5.0].copy()

    print(f"Drug–target edges: {len(dt)} -> {len(dt_filtered)} after score >= 5.0")

    dt_out.parent.mkdir(parents=True, exist_ok=True)
    dt_filtered.to_csv(dt_out, index=False)
    print(f"Wrote filtered drug_targets to {dt_out.resolve()}")

    print(f"\nLoading {gd_in} ...")
    gd = pd.read_csv(gd_in)

    if "score" not in gd.columns:
        raise SystemExit("gene_disease.csv has no 'score' column.")

    gd = gd.dropna(subset=["score"])
    gd_filtered = gd[gd["score"] >= 0.05].copy()

    print(f"Gene–disease edges: {len(gd)} -> {len(gd_filtered)} after score >= 0.01")

    gd_out.parent.mkdir(parents=True, exist_ok=True)
    gd_filtered.to_csv(gd_out, index=False)
    print(f"Wrote filtered gene_disease to {gd_out.resolve()}")


if __name__ == "__main__":
    main()
