#!/usr/bin/env python
"""
Normalize a raw drug–target CSV into the canonical format used by DDH.

Canonical outputs:

    data/real/drug_targets.csv
        columns: drug_id, drug_name, gene_id, score (optional)

    data/real/drugs.csv
        columns: drug_id, drug_name

The script is intentionally generic: you tell it which input columns correspond
to drug ID, drug name, target gene/protein ID, and an optional score/confidence.

Example usage:

    python scripts/prepare_drug_targets_from_csv.py \
        --input data/raw/drug_targets_raw.csv \
        --output-drug-targets data/real/drug_targets.csv \
        --output-drugs data/real/drugs.csv \
        --drug-id-col drug_chembl_id \
        --drug-name-col pref_name \
        --gene-id-col target_gene_symbol \
        --score-col confidence_score
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import pandas as pd


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize a raw drug–target CSV into DDH canonical format.",
    )

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the raw drug–target CSV file.",
    )
    parser.add_argument(
        "--output-drug-targets",
        type=str,
        required=True,
        help="Path to write the normalized drug_targets.csv file.",
    )
    parser.add_argument(
        "--output-drugs",
        type=str,
        required=True,
        help="Path to write the normalized drugs.csv file.",
    )

    # Column mappings
    parser.add_argument(
        "--drug-id-col",
        type=str,
        required=True,
        help="Name of the column containing drug IDs in the raw CSV.",
    )
    parser.add_argument(
        "--drug-name-col",
        type=str,
        required=True,
        help="Name of the column containing drug names in the raw CSV.",
    )
    parser.add_argument(
        "--gene-id-col",
        type=str,
        required=True,
        help="Name of the column containing target gene/protein IDs in the raw CSV.",
    )
    parser.add_argument(
        "--score-col",
        type=str,
        default=None,
        help="Optional: name of the column containing a numeric score/confidence.",
    )

    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = parse_args(argv)

    input_path = Path(args.input)
    if not input_path.exists():
        raise SystemExit(f"Input file does not exist: {input_path}")

    df_raw = pd.read_csv(
    input_path,
    sep=None,            
    engine="python",     
    #on_bad_lines="skip", 
)
    print("Loaded raw file with columns:")
    print(df_raw.columns.tolist()[:50])
    
    required_cols = {args.drug_id_col, args.drug_name_col, args.gene_id_col}
    missing = required_cols - set(df_raw.columns)
    if missing:
        raise SystemExit(f"Missing required columns in input: {missing}")

    # Build normalized drug_targets table
    out_cols = {
        "drug_id": args.drug_id_col,
        "drug_name": args.drug_name_col,
        "gene_id": args.gene_id_col,
    }

    df_dt = df_raw[list(out_cols.values())].rename(columns={v: k for k, v in out_cols.items()})

    # Add score column if requested and present
    if args.score_col is not None and args.score_col in df_raw.columns:
        df_dt["score"] = df_raw[args.score_col]
    else:
        # Leave score column absent; the pipeline can handle missing scores
        pass

    # Drop rows that are missing critical IDs
    df_dt = df_dt.dropna(subset=["drug_id", "gene_id"])
    
    # Fill missing drug_name with the ID as a fallback
    df_dt["drug_name"] = df_dt["drug_name"].fillna(df_dt["drug_id"].astype(str))

    df_dt["drug_id"] = df_dt["drug_id"].astype(str).str.strip()
    df_dt["drug_name"] = df_dt["drug_name"].astype(str).str.strip()
    df_dt["gene_id"] = df_dt["gene_id"].astype(str).str.strip()


    # Optional: aggregate duplicates (same drug_id, gene_id)
    # Here, we take the max score if multiple rows exist, or just drop duplicates.
    if "score" in df_dt.columns:
        df_dt = (
            df_dt.groupby(["drug_id", "drug_name", "gene_id"], as_index=False)
            .agg({"score": "max"})
        )
    else:
        df_dt = df_dt.drop_duplicates(subset=["drug_id", "drug_name", "gene_id"])

    # Build normalized drugs table
    df_drugs = (
        df_dt[["drug_id", "drug_name"]]
        .drop_duplicates()
        .sort_values(by=["drug_id"])
        .reset_index(drop=True)
    )

    # Write outputs
    output_dt_path = Path(args.output_drug_targets)
    output_drugs_path = Path(args.output_drugs)

    output_dt_path.parent.mkdir(parents=True, exist_ok=True)
    output_drugs_path.parent.mkdir(parents=True, exist_ok=True)

    df_dt.to_csv(output_dt_path, index=False)
    df_drugs.to_csv(output_drugs_path, index=False)

    print(f"Normalized drug_targets written to: {output_dt_path.resolve()}")
    print(f"Normalized drugs written to: {output_drugs_path.resolve()}")


if __name__ == "__main__":
    main()

