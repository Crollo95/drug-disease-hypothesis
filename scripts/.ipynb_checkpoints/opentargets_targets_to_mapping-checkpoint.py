#!/usr/bin/env python
"""
Convert Open Targets Target parquet files into a mapping table.

Input:
    data/raw/targets/*.parquet  (directory with OT Target dataset)

Output:
    data/real/target_mapping.csv

Columns in output:
    ensembl_id   : Ensembl gene ID (primary key in OT target)
    label        : any label that can be used to match (symbol, name, synonym, etc.)

We will later match ChEMBL target names against 'label' (case-insensitive).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, List

import pandas as pd


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build Ensembl ID mapping from Open Targets Target parquet files.",
    )

    parser.add_argument(
        "--input-dir",
        type=str,
        required=True,
        help="Directory containing OT Target parquet files.",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to write target_mapping.csv.",
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

    print("Columns in OT Target dataset:")
    print(df.columns.tolist())

    # Identify key columns with some flexibility
    # Ensembl gene ID
    if "id" in df.columns:
        ensembl_col = "id"
    elif "targetId" in df.columns:
        ensembl_col = "targetId"
    else:
        raise SystemExit("Could not find Ensembl ID column (expected 'id' or 'targetId').")

    # Symbol
    symbol_col = None
    for cand in ["approvedSymbol", "symbol"]:
        if cand in df.columns:
            symbol_col = cand
            break

    # Name
    name_col = None
    for cand in ["approvedName", "name"]:
        if cand in df.columns:
            name_col = cand
            break

    # UniProt IDs (may be a single string or list)
    uniprot_col = None
    for cand in ["uniprotId", "uniprot_id", "uniprotIds"]:
        if cand in df.columns:
            uniprot_col = cand
            break

    # Synonyms (often a list-like column)
    synonyms_col = None
    for cand in ["synonyms", "targetSynonyms"]:
        if cand in df.columns:
            synonyms_col = cand
            break

    records = []

    records = []

    def add_label(ensembl_id: str, label: Optional[str]) -> None:
        if label is None:
            return
        label = str(label).strip()
        if not label:
            return
        records.append({"ensembl_id": ensembl_id, "label": label})

    for _, row in df.iterrows():
        ensg = str(row[ensembl_col]).strip()
        if not ensg:
            continue

        # Symbol
        if symbol_col is not None:
            val = row.get(symbol_col, None)
            if isinstance(val, str) and val.strip():
                add_label(ensg, val)

        # Approved name
        if name_col is not None:
            val = row.get(name_col, None)
            if isinstance(val, str) and val.strip():
                add_label(ensg, val)

        # UniProt IDs (may be list or string)
        if uniprot_col is not None:
            val = row.get(uniprot_col, None)
            if isinstance(val, (list, tuple, set)):
                for v in val:
                    add_label(ensg, v)
            elif isinstance(val, str) and val.strip():
                for v in val.split(","):
                    add_label(ensg, v)

        # Synonyms (list-like or string)
        if synonyms_col is not None:
            val = row.get(synonyms_col, None)
            if isinstance(val, (list, tuple, set)):
                for v in val:
                    add_label(ensg, v)
            elif isinstance(val, str) and val.strip():
                # sometimes may be a comma-separated string or similar
                add_label(ensg, val)

    if not records:
        raise SystemExit("No mapping records built; check target dataset schema.")

    df_map = pd.DataFrame(records)
    # Normalize labels for safer matching
    df_map["label"] = df_map["label"].astype(str).str.strip()
    df_map["label_lower"] = df_map["label"].str.lower()

    # Drop duplicate label–ensembl combinations
    df_map = df_map.drop_duplicates(subset=["label_lower", "ensembl_id"])

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_map.to_csv(output_path, index=False)

    print(f"target_mapping.csv written to: {output_path.resolve()}")
    print(f"Total mappings: {len(df_map)}")
    print(f"Unique labels: {df_map['label_lower'].nunique()}")
    print(f"Unique Ensembl IDs: {df_map['ensembl_id'].nunique()}")


if __name__ == "__main__":
    main()
