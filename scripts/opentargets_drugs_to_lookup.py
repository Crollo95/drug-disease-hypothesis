#!/usr/bin/env python
"""
Build a drug lookup table from Open Targets Drug/Clinical candidates parquet files.

Input:
    data/raw/ot_drugs/*.parquet

Output:
    data/real/drug_lookup.csv with:
        - drug_id   (usually ChEMBL ID)
        - drug_name (preferred/common name if available)
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, List

import pandas as pd


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert Open Targets Drug/Clinical candidates parquet to simple ID→name CSV."
    )
    p.add_argument(
        "--input-dir",
        type=str,
        required=True,
        help="Directory with Drug/Clinical candidates parquet files.",
    )
    p.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to write drug_lookup.csv",
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

    print(f"Found {len(parquet_files)} drug parquet files. Reading...")
    dfs = [pd.read_parquet(p) for p in parquet_files]
    df = pd.concat(dfs, ignore_index=True)

    print("Columns in Drug/Clinical candidates dataset:")
    print(df.columns.tolist())

    # Try to identify ID and name columns
    candidate_id_cols = ["id", "chemblId", "drugId"]
    candidate_name_cols = ["name", "prefName", "approvedName"]

    id_col = next((c for c in candidate_id_cols if c in df.columns), None)
    name_col = next((c for c in candidate_name_cols if c in df.columns), None)

    if id_col is None:
        raise SystemExit("Could not find a drug ID column in dataset.")
    if name_col is None:
        raise SystemExit("Could not find a drug name column in dataset.")

    print(f"Using '{id_col}' as drug_id, '{name_col}' as drug_name")

    lookup = df[[id_col, name_col]].dropna().drop_duplicates()
    lookup = lookup.rename(columns={id_col: "drug_id", name_col: "drug_name"})

    lookup["drug_id"] = lookup["drug_id"].astype(str).str.strip()
    lookup["drug_name"] = lookup["drug_name"].astype(str).str.strip()

    lookup.to_csv(out_path, index=False)
    print(f"Drug lookup written to: {out_path.resolve()}")
    print(f"Rows: {len(lookup)}")


if __name__ == "__main__":
    main()
