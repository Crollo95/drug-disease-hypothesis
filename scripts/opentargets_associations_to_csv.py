#!/usr/bin/env python
"""
Convert Open Targets 'Associations - direct (overall score)' parquet files
into a flat CSV suitable as gene_disease_raw.csv.

We assume the dataset has at least:
    - targetId   : Ensembl gene ID
    - diseaseId  : EFO disease ID
    - score      : overall association score

Some versions may also include:
    - diseaseFromSource : human-readable disease name

Usage:

    python scripts/opentargets_associations_to_csv.py \
        --input-dir data/raw/opentargets_associations_direct_overall \
        --output data/raw/gene_disease_raw.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, List

import pandas as pd


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert Open Targets association parquet files to a flat CSV.",
    )

    parser.add_argument(
        "--input-dir",
        type=str,
        required=True,
        help="Directory containing Open Targets association parquet files.",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to write the gene_disease_raw.csv file.",
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

    dfs = []
    for p in parquet_files:
        df_part = pd.read_parquet(p)
        dfs.append(df_part)

    df = pd.concat(dfs, ignore_index=True)
    print("Columns in Open Targets associations:")
    print(df.columns.tolist())

    # Try to infer relevant columns
    # gene/target ID
    if "targetId" in df.columns:
        gene_col = "targetId"
    elif "target_id" in df.columns:
        gene_col = "target_id"
    else:
        raise SystemExit("Could not find targetId or target_id in columns.")

    # disease ID
    if "diseaseId" in df.columns:
        disease_col = "diseaseId"
    elif "disease_id" in df.columns:
        disease_col = "disease_id"
    else:
        raise SystemExit("Could not find diseaseId or disease_id in columns.")

    # score
    if "score" in df.columns:
        score_col = "score"
    elif "overallScore" in df.columns:
        score_col = "overallScore"
    else:
        raise SystemExit("Could not find 'score' or 'overallScore' in columns.")

    # disease name (optional)
    if "diseaseFromSource" in df.columns:
        disease_name_col = "diseaseFromSource"
    elif "diseaseName" in df.columns:
        disease_name_col = "diseaseName"
    else:
        disease_name_col = None  # we'll fall back to disease_id

    out = {}

    out["gene_id"] = df[gene_col].astype(str).str.strip()
    out["disease_id"] = df[disease_col].astype(str).str.strip()
    out["score"] = df[score_col]

    if disease_name_col is not None:
        out["disease_name"] = df[disease_name_col].astype(str).str.strip()
    else:
        # fallback: use disease_id as name
        out["disease_name"] = out["disease_id"]

    df_out = pd.DataFrame(out)

    # Drop rows missing IDs
    df_out = df_out.dropna(subset=["gene_id", "disease_id"])

    # Deduplicate (if needed) by taking max score per pair
    df_out = (
        df_out.groupby(["gene_id", "disease_id", "disease_name"], as_index=False)
        .agg({"score": "max"})
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df_out.to_csv(output_path, index=False)
    print(f"gene_disease_raw.csv written to: {output_path.resolve()}")


if __name__ == "__main__":
    main()

