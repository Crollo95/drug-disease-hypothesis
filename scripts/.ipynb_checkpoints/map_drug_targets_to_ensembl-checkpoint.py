#!/usr/bin/env python
"""
Map ChEMBL target names in data/real/drug_targets.csv to Ensembl gene IDs
using data/real/target_mapping.csv.

Input:
    data/real/drug_targets.csv
        columns: drug_id, drug_name, gene_id, score
        where gene_id = target protein name (from ChEMBL)

    data/real/target_mapping.csv
        columns: ensembl_id, label, label_lower

Output:
    data/real/drug_targets_mapped.csv
        columns: drug_id, drug_name, gene_id, score
        where gene_id = Ensembl ID (ENSG...)

Unmapped targets are dropped (we can change this later if needed).
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd


def normalize_label(s: str) -> str:
    return str(s).strip().lower()


def main() -> None:
    base = Path("data/real")
    dt_path = base / "drug_targets.csv"
    map_path = base / "target_mapping.csv"
    out_path = base / "drug_targets_mapped.csv"

    if not dt_path.exists():
        raise SystemExit(f"Missing {dt_path}")
    if not map_path.exists():
        raise SystemExit(f"Missing {map_path}")

    dt = pd.read_csv(dt_path)
    mapping = pd.read_csv(map_path)

    if "label_lower" not in mapping.columns:
        mapping["label_lower"] = mapping["label"].astype(str).str.lower()

    # Normalize gene_id from drug_targets for matching
    dt["gene_label"] = dt["gene_id"].astype(str).str.strip()
    dt["gene_label_lower"] = dt["gene_label"].str.lower()

    # Exact case-insensitive match
    merged = dt.merge(
        mapping[["label_lower", "ensembl_id"]],
        left_on="gene_label_lower",
        right_on="label_lower",
        how="left",
    )

    total = len(merged)
    mapped = merged["ensembl_id"].notna().sum()
    print(f"Total drug–target rows: {total}")
    print(f"Successfully mapped to Ensembl: {mapped} ({mapped / total:.1%})")

    # Keep only mapped rows for now
    mapped_df = merged[merged["ensembl_id"].notna()].copy()

    # Replace gene_id with Ensembl
    mapped_df["gene_id"] = mapped_df["ensembl_id"]

    # Select canonical columns
    out_cols = ["drug_id", "drug_name", "gene_id", "score"]
    out_df = mapped_df[out_cols].drop_duplicates()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)

    print(f"Mapped drug_targets written to: {out_path.resolve()}")
    print(f"Rows in mapped drug_targets: {len(out_df)}")


if __name__ == "__main__":
    main()
