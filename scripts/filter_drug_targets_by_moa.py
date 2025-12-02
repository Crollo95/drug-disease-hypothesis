#!/usr/bin/env python
"""
Filter drug_targets by Mechanism of Action (MoA) pairs.

Inputs:
    data/real/drug_targets_mapped.csv (or any drug_targets with ENSG gene_id)
    data/real/moa_pairs_ensembl.csv   (drug_id, gene_id from MoA)

Outputs:
    data/real/drug_targets_moa.csv    (drug_targets restricted to MoA-backed edges)
"""

from __future__ import annotations

from pathlib import Path
import argparse
from typing import Optional, List

import pandas as pd


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Filter drug_targets by MoA-supported (drug_id,gene_id) pairs."
    )
    p.add_argument(
        "--drug-targets",
        type=str,
        required=True,
        help="Input drug_targets CSV (e.g. data/real/drug_targets_mapped.csv).",
    )
    p.add_argument(
        "--moa-pairs",
        type=str,
        required=True,
        help="MoA pairs CSV (drug_id,gene_id), e.g. data/real/moa_pairs_ensembl.csv.",
    )
    p.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output CSV for MoA-filtered drug_targets.",
    )
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)

    dt_path = Path(args.drug_targets)
    moa_path = Path(args.moa_pairs)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading drug_targets from {dt_path} ...")
    dt = pd.read_csv(dt_path)
    print(f"Rows in drug_targets: {len(dt)}")

    print(f"Loading MoA pairs from {moa_path} ...")
    moa = pd.read_csv(moa_path)
    print(f"Rows in MoA pairs: {len(moa)}")

    # Normalize IDs
    for c in ["drug_id", "gene_id"]:
        if c in dt.columns:
            dt[c] = dt[c].astype(str).str.strip()
        if c in moa.columns:
            moa[c] = moa[c].astype(str).str.strip()

    # Keep only (drug_id, gene_id) that appear in MoA
    moa_pairs = moa[["drug_id", "gene_id"]].drop_duplicates()
    moa_pairs["moa_flag"] = True

    dt_moa = dt.merge(moa_pairs, on=["drug_id", "gene_id"], how="inner")
    print(f"Rows in drug_targets after MoA filter: {len(dt_moa)}")
    n_drugs = dt_moa["drug_id"].nunique()
    print(f"Unique drugs with MoA-backed targets: {n_drugs}")

    dt_moa.to_csv(out_path, index=False)
    print(f"MoA-filtered drug_targets written to: {out_path}")


if __name__ == "__main__":
    main()
