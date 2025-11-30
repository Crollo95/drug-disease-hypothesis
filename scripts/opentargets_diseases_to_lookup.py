#!/usr/bin/env python
"""
Build a disease lookup table from Open Targets Disease/Phenotype parquet files.

Input:
    data/raw/diseases/*.parquet  (Open Targets Disease/Phenotype dataset)

Output:
    data/real/disease_lookup.csv  with columns:
        - disease_id   (Open Targets diseaseId / EFO ID)
        - disease_name (human-readable label)
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, List

import pandas as pd


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert Open Targets Disease/Phenotype parquet to a simple ID→name CSV."
    )
    p.add_argument(
        "--input-dir",
        type=str,
        required=True,
        help="Directory with Disease/Phenotype parquet files.",
    )
    p.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to write disease_lookup.csv",
    )
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)

    input_dir = Path(args.input_dir)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    parquet_files = sorted(input_dir.glob("*.parquet"))
    if not parquet_files:
        raise SystemExit(f"No parquet files found in {input_dir}")

    print(f"Found {len(parquet_files)} disease parquet files. Reading...")
    dfs = [pd.read_parquet(p) for p in parquet_files]
    df = pd.concat(dfs, ignore_index=True)

    print("Columns in Disease/Phenotype dataset:")
    print(df.columns.tolist())

    # Try to identify ID and name columns robustly
    candidate_id_cols = ["id", "efoId", "diseaseId"]
    candidate_name_cols = ["name", "label", "diseaseFromSource", "diseaseName"]

    id_col = next((c for c in candidate_id_cols if c in df.columns), None)
    name_col = next((c for c in candidate_name_cols if c in df.columns), None)

    if id_col is None:
        raise SystemExit("Could not find a disease ID column in dataset.")
    if name_col is None:
        raise SystemExit("Could not find a disease name/label column in dataset.")

    print(f"Using '{id_col}' as disease_id, '{name_col}' as disease_name")

    lookup = df[[id_col, name_col]].dropna().drop_duplicates()
    lookup = lookup.rename(columns={id_col: "disease_id", name_col: "disease_name"})

    # Normalize strings
    lookup["disease_id"] = lookup["disease_id"].astype(str).str.strip()
    lookup["disease_name"] = lookup["disease_name"].astype(str).str.strip()

    lookup.to_csv(out_path, index=False)
    print(f"Disease lookup written to: {out_path.resolve()}")
    print(f"Rows: {len(lookup)}")


if __name__ == "__main__":
    main()
