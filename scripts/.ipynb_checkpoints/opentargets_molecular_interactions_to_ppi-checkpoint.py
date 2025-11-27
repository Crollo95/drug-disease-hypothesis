#!/usr/bin/env python
"""
Convert Open Targets 'Molecular interactions' parquet files into a PPI CSV.

We treat each row as an undirected interaction between two targets (genes)
identified by Ensembl IDs.

You MUST specify which columns in the dataset contain the two target IDs,
and optionally a score/weight.

Example usage (you will adapt column names after seeing printed columns):

    python scripts/opentargets_molecular_interactions_to_ppi.py \
        --input-dir data/raw/molecular_interactions \
        --output data/real/ppi.csv \
        --gene1-col targetAId \
        --gene2-col targetBId \
        --weight-col score \
        --min-weight 0.4
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, List

import pandas as pd


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build PPI CSV from Open Targets Molecular Interactions parquet files.",
    )

    parser.add_argument(
        "--input-dir",
        type=str,
        required=True,
        help="Directory containing OT Molecular Interactions parquet files.",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to write ppi.csv.",
    )

    parser.add_argument(
        "--gene1-col",
        type=str,
        required=True,
        help="Column name for first target (e.g. targetAId).",
    )
    parser.add_argument(
        "--gene2-col",
        type=str,
        required=True,
        help="Column name for second target (e.g. targetBId).",
    )
    parser.add_argument(
        "--weight-col",
        type=str,
        default=None,
        help="Optional column name for interaction score/weight.",
    )
    parser.add_argument(
        "--min-weight",
        type=float,
        default=None,
        help="Optional minimum weight threshold (applied if weight-col is set).",
    )

    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = parse_args(argv)

    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        raise SystemExit(f"Input directory does not exist: {input_dir}")

    parquet_files: List[Path] = sorted(input_dir.glob("*.parquet"))
    if not parquet_files:
        raise SystemExit(f"No .parquet files found in {input_dir}")

    print(f"Found {len(parquet_files)} parquet files. Reading...")

    dfs = [pd.read_parquet(p) for p in parquet_files]
    df = pd.concat(dfs, ignore_index=True)

    print("Columns in Molecular Interactions dataset:")
    print(df.columns.tolist())

    for col in [args.gene1_col, args.gene2_col]:
        if col not in df.columns:
            raise SystemExit(f"Column '{col}' not found in dataset.")

    if args.weight_col is not None and args.weight_col not in df.columns:
        raise SystemExit(f"Weight column '{args.weight_col}' not found in dataset.")

    # Basic selection
    cols = [args.gene1_col, args.gene2_col]
    if args.weight_col is not None:
        cols.append(args.weight_col)

    df_ppi = df[cols].copy()

    # Rename columns to canonical names
    rename_map = {
        args.gene1_col: "gene1_id",
        args.gene2_col: "gene2_id",
    }
    if args.weight_col is not None:
        rename_map[args.weight_col] = "weight"

    df_ppi = df_ppi.rename(columns=rename_map)

    # Drop missing IDs
    df_ppi = df_ppi.dropna(subset=["gene1_id", "gene2_id"])

    # Convert to str, strip, and keep only ENSG* to focus on human genes
    df_ppi["gene1_id"] = df_ppi["gene1_id"].astype(str).str.strip()
    df_ppi["gene2_id"] = df_ppi["gene2_id"].astype(str).str.strip()

    df_ppi = df_ppi[
        df_ppi["gene1_id"].str.startswith("ENSG")
        & df_ppi["gene2_id"].str.startswith("ENSG")
    ]

    # Handle weight if present
    if "weight" in df_ppi.columns:
        df_ppi = df_ppi.dropna(subset=["weight"])
        df_ppi["weight"] = pd.to_numeric(df_ppi["weight"], errors="coerce")
        df_ppi = df_ppi.dropna(subset=["weight"])

        if args.min_weight is not None:
            df_ppi = df_ppi[df_ppi["weight"] >= args.min_weight]

    else:
        # No weight: just set all weights to 1.0
        df_ppi["weight"] = 1.0

    # Remove self-loops
    df_ppi = df_ppi[df_ppi["gene1_id"] != df_ppi["gene2_id"]]

    # Drop duplicate undirected edges (A,B) ~ (B,A)
    df_ppi["min_id"] = df_ppi[["gene1_id", "gene2_id"]].min(axis=1)
    df_ppi["max_id"] = df_ppi[["gene1_id", "gene2_id"]].max(axis=1)
    df_ppi = (
        df_ppi.groupby(["min_id", "max_id"], as_index=False)
        .agg({"weight": "max"})
        .rename(columns={"min_id": "gene1_id", "max_id": "gene2_id"})
    )

    # Write output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_ppi.to_csv(output_path, index=False)

    print(f"ppi.csv written to: {output_path.resolve()}")
    print(f"Total interactions (edges): {len(df_ppi)}")


if __name__ == "__main__":
    main()
