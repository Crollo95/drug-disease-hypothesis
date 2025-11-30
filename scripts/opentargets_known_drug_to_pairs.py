#!/usr/bin/env python
"""
Convert Open Targets 'Known drug' parquet into a simple list of
(drug_id, disease_id) pairs that represent known indications.

Input:
    data/raw/known_drug/*.parquet   (Known drug dataset)

Output:
    data/real/known_indications.csv  with columns:
        - drug_id
        - disease_id
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, List

import pandas as pd


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert Open Targets Known drug parquet to drug–disease pairs."
    )
    p.add_argument(
        "--input-dir",
        type=str,
        required=True,
        help="Directory with Known drug parquet files.",
    )
    p.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to write known_indications.csv",
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

    print(f"Found {len(parquet_files)} parquet files. Reading...")
    dfs = [pd.read_parquet(p) for p in parquet_files]
    df = pd.concat(dfs, ignore_index=True)

    print("Columns in Known drug dataset:")
    print(df.columns.tolist())

    # Try to infer columns for drug and disease IDs
    candidate_drug_cols = ["drugId", "id", "chemblId", "drug_id"]
    candidate_disease_cols = ["diseaseId", "diseaseFromSourceMappedId", "efoId", "disease_id"]

    drug_col = next((c for c in candidate_drug_cols if c in df.columns), None)
    disease_col = next((c for c in candidate_disease_cols if c in df.columns), None)

    if drug_col is None:
        raise SystemExit("Could not find a drug ID column in Known drug dataset.")
    if disease_col is None:
        raise SystemExit("Could not find a disease ID column in Known drug dataset.")

    print(f"Using '{drug_col}' as drug_id, '{disease_col}' as disease_id")

    pairs = df[[drug_col, disease_col]].dropna().drop_duplicates()
    pairs = pairs.rename(columns={drug_col: "drug_id", disease_col: "disease_id"})

    pairs["drug_id"] = pairs["drug_id"].astype(str).str.strip()
    pairs["disease_id"] = pairs["disease_id"].astype(str).str.strip()

    pairs.to_csv(out_path, index=False)
    print(f"Wrote known indications to: {out_path.resolve()}")
    print(f"Unique pairs: {len(pairs)}")


if __name__ == "__main__":
    main()
