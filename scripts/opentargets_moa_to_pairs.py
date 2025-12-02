#!/usr/bin/env python
"""
Convert Open Targets 'Drug - mechanism of action' dataset into a simple
(drug_id, gene_id) table in Ensembl ID space.

Input:
    data/raw/drug_moa/*.parquet

Outputs:
    data/real/moa_pairs_raw.csv       (drug_id, target_raw)
    data/real/moa_pairs_ensembl.csv   (drug_id, gene_id as ENSG...)
"""

from __future__ import annotations

from pathlib import Path
import argparse
from typing import Optional, List, Any, Dict

import pandas as pd
import numpy as np

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Convert OT Drug - mechanism of action parquet to drug–gene pairs."
    )
    p.add_argument(
        "--input-dir",
        type=str,
        required=True,
        help="Directory with Drug - mechanism of action parquet files.",
    )
    p.add_argument(
        "--target-mapping",
        type=str,
        required=True,
        help="CSV with target mapping (e.g. data/real/target_mapping.csv) "
             "with at least 'ensembl_id' and some target id columns.",
    )
    p.add_argument(
        "--output-raw",
        type=str,
        required=True,
        help="Path to write moa_pairs_raw.csv",
    )
    p.add_argument(
        "--output-ensembl",
        type=str,
        required=True,
        help="Path to write moa_pairs_ensembl.csv",
    )
    return p.parse_args(argv)


def extract_target_id(t: Any) -> Optional[str]:
    """
    Extract a target identifier string from a target object.

    'targets' entries can be:
      - strings (e.g. 'ENSG00000162409')
      - dicts with fields like 'id', 'targetFromSourceId', etc.
    """
    if t is None:
        return None

    # Case 1: already a string-like ID
    if isinstance(t, str):
        val = t.strip()
        return val if val else None

    # Case 2: dict with known keys
    if isinstance(t, dict):
        candidate_keys = [
            "id",
            "targetFromSourceId",
            "targetFromSource",
            "ensemblId",
            "geneId",
        ]
        for k in candidate_keys:
            if k in t and t[k] is not None:
                return str(t[k]).strip()

    return None


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)

    input_dir = Path(args.input_dir)
    mapping_path = Path(args.target_mapping)
    out_raw = Path(args.output_raw)
    out_ens = Path(args.output_ensembl)

    out_raw.parent.mkdir(parents=True, exist_ok=True)
    out_ens.parent.mkdir(parents=True, exist_ok=True)

    parquet_files = sorted(input_dir.glob("*.parquet"))
    if not parquet_files:
        raise SystemExit(f"No parquet files found in {input_dir}")

    print(f"Found {len(parquet_files)} MoA parquet files. Reading...")
    dfs = [pd.read_parquet(p) for p in parquet_files]
    df = pd.concat(dfs, ignore_index=True)

    print("Columns in MoA dataset:")
    print(df.columns.tolist())

    if "chemblIds" not in df.columns or "targets" not in df.columns:
        raise SystemExit("Expected 'chemblIds' and 'targets' columns in MoA dataset.")

    pairs: List[Dict[str, str]] = []

    for idx, row in df.iterrows():
        chembl_ids = row["chemblIds"]
        targets = row["targets"]

        # --- Normalize chembl_ids to a list of strings ---
        if isinstance(chembl_ids, np.ndarray):
            chembl_ids_seq = chembl_ids.tolist()
        elif isinstance(chembl_ids, (list, tuple)):
            chembl_ids_seq = chembl_ids
        elif isinstance(chembl_ids, str):
            chembl_ids_seq = [chembl_ids]
        else:
            # unknown / missing format
            continue

        drug_ids = [str(d).strip() for d in chembl_ids_seq if d is not None]

        if not drug_ids:
            continue

        # --- Normalize targets to a list of items ---
        if isinstance(targets, np.ndarray):
            targets_seq = targets.tolist()
        elif isinstance(targets, (list, tuple)):
            targets_seq = targets
        else:
            # missing or unexpected format
            continue

        for tgt in targets_seq:
            tgt_id = extract_target_id(tgt)
            if tgt_id is None:
                continue
            for drug_id in drug_ids:
                pairs.append({"drug_id": drug_id, "target_raw": tgt_id})

        if (idx + 1) % 1000 == 0:
            print(f"Processed {idx+1}/{len(df)} rows...", end="\r")


    print(f"\nExtracted {len(pairs)} raw (drug_id, target_raw) pairs from MoA.")

    if not pairs:
        print("Warning: no MoA pairs extracted; check dataset structure.")
        return

    raw_df = pd.DataFrame(pairs).drop_duplicates()
    raw_df["drug_id"] = raw_df["drug_id"].astype(str).str.strip()
    raw_df["target_raw"] = raw_df["target_raw"].astype(str).str.strip()

    raw_df.to_csv(out_raw, index=False)
    print(f"Raw MoA pairs written to: {out_raw} (rows: {len(raw_df)})")

    # --------------------------
    # Map target_raw -> Ensembl
    # --------------------------
    if not mapping_path.exists():
        print(f"Target mapping {mapping_path} not found; cannot map to Ensembl.")
        return

    print(f"Loading target mapping from {mapping_path} ...")
    tm = pd.read_csv(mapping_path)

    if "ensembl_id" not in tm.columns:
        raise SystemExit("target_mapping.csv must have an 'ensembl_id' column.")

    candidate_id_cols = [
        "target_id",
        "targetFromSourceId",
        "targetFromSource",
        "approvedSymbol",
        "symbol",
        "hgnc_symbol",
        "preferredName",
        "uniprot_id",
    ]
    id_cols = [c for c in candidate_id_cols if c in tm.columns]

    if not id_cols:
        print("No obvious target identifier columns in target_mapping; "
              "writing only raw pairs.")
        return

    print(f"Using mapping ID columns: {id_cols}")

    long_parts = []
    for c in id_cols:
        sub = tm[[c, "ensembl_id"]].dropna().copy()
        sub = sub.rename(columns={c: "target_raw"})
        long_parts.append(sub)

    long_map = pd.concat(long_parts, ignore_index=True).drop_duplicates()
    long_map["target_raw"] = long_map["target_raw"].astype(str).str.strip()
    long_map["ensembl_id"] = long_map["ensembl_id"].astype(str).str.strip()

    print(f"Long mapping rows: {len(long_map)}")

    merged = raw_df.merge(long_map, on="target_raw", how="inner")
    merged = merged.rename(columns={"ensembl_id": "gene_id"})
    merged = merged[["drug_id", "gene_id"]].drop_duplicates()

    merged.to_csv(out_ens, index=False)
    print(f"Ensembl-mapped MoA pairs written to: {out_ens} (rows: {len(merged)})")


if __name__ == "__main__":
    main()
