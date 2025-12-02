#!/usr/bin/env python
"""
Quick sanity checks for raw input data under data/raw/.

- Verifies that expected directories exist
- Checks that at least one file is present per dataset
- Reads one file and checks required columns (for core datasets)

Run from project root:

    python scripts/check_raw_data.py
"""

from __future__ import annotations

from pathlib import Path
import sys
from typing import List

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"

DATASETS = [
    # ----------------- CORE (required) -----------------
    {
        "name": "chembl_activities",
        "path": RAW_DIR / "chembl_activities",
        "pattern": "chembl_activities_part*.csv",
        "required_cols": [
            "Molecule ChEMBL ID",
            "Molecule Name",
            "Target Name",
            "pChEMBL Value",
        ],
        "reader": "csv",
        "required": True,
    },
    {
        "name": "association_overall_direct",
        "path": RAW_DIR / "association_overall_direct",
        "pattern": "part-*.parquet",
        "required_cols": [
            "targetId",
            "diseaseId",
            "score",
        ],
        "reader": "parquet",
        "required": True,
    },
    {
        "name": "targets",
        "path": RAW_DIR / "targets",
        "pattern": "part-*.parquet",
        "required_cols": [
            "id",
            "approvedSymbol",
        ],
        "reader": "parquet",
        "required": True,
    },
    {
        "name": "molecular_interactions",
        "path": RAW_DIR / "molecular_interactions",
        "pattern": "part-*.parquet",
        "required_cols": [
            "targetA",
            "targetB",
            "scoring",
        ],
        "reader": "parquet",
        "required": True,
    },

    # ----------------- OPTIONAL (recommended/extra) -----------------
    {
        "name": "drug_moa",
        "path": RAW_DIR / "drug_moa",
        "pattern": "part-*.parquet",
        "required_cols": [
            "chemblIds",
            "targets",
        ],
        "reader": "parquet",
        "required": False,
    },
    {
        "name": "diseases",
        "path": RAW_DIR / "diseases",
        "pattern": "disease.parquet",
        "required_cols": [
            "id",  # OT disease/phenotype id
        ],
        "reader": "parquet",
        "required": False,
    },
    {
        "name": "ot_drugs",
        "path": RAW_DIR / "ot_drugs",
        "pattern": "part-*.parquet",
        "required_cols": [
            "id",  # OT drug id
        ],
        "reader": "parquet",
        "required": False,
    },
    {
        "name": "known_drug",
        "path": RAW_DIR / "known_drug",
        "pattern": "part-*.parquet",
        "required_cols": [
            "drugId",
            "diseaseId",
        ],
        "reader": "parquet",
        "required": False,  # recommended but not strictly required
    },
    {
        "name": "drug_indications",
        "path": RAW_DIR / "drug_indications",
        "pattern": "part-*.parquet",
        "required_cols": [
            "id",  # drug id, indications are nested
        ],
        "reader": "parquet",
        "required": False,
    },
]


def read_sample_file(sample_file: Path, reader: str, dataset_name: str) -> pd.DataFrame:
    """Read a small sample of a file, with some robustness for CSV dialects."""
    if reader == "csv":
        # First try default comma-separated
        df = pd.read_csv(sample_file, nrows=5)
        # Heuristic: if we got a single giant column with semicolons in the name,
        # re-read as semicolon-separated.
        if len(df.columns) == 1 and ";" in df.columns[0]:
            print(f"  ℹ️ Detected semicolon-separated CSV for {dataset_name}, retrying with sep=';'.")
            df = pd.read_csv(sample_file, nrows=5, sep=";")
        return df
    elif reader == "parquet":
        return pd.read_parquet(sample_file)
    else:
        raise ValueError(f"Unknown reader type: {reader}")


def check_dataset(
    name: str,
    path: Path,
    pattern: str,
    required_cols: List[str],
    reader: str,
    required: bool,
) -> bool:
    print(f"\nChecking dataset: {name}")
    if not path.exists():
        msg = f"Directory not found: {path}"
        if required:
            print(f"  ❌ {msg}")
            return False
        else:
            print(f"  ⚠️ (optional) {msg}")
            return True

    files = sorted(path.glob(pattern))
    if not files:
        msg = f"No files matching pattern '{pattern}' in {path}"
        if required:
            print(f"  ❌ {msg}")
            return False
        else:
            print(f"  ⚠️ (optional) {msg}")
            return True

    sample_file = files[0]
    print(f"  ✅ Found {len(files)} file(s). Sampling: {sample_file.name}")

    df = read_sample_file(sample_file, reader=reader, dataset_name=name)

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        msg = f"Missing required columns: {missing}\n  📎 Columns present: {list(df.columns)}"
        if required:
            print(f"  ❌ {msg}")
            return False
        else:
            print(f"  ⚠️ (optional) {msg}")
            return True

    print(f"  ✅ All required columns present.")
    return True


def main() -> None:
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Raw data dir: {RAW_DIR}")

    if not RAW_DIR.exists():
        print("❌ data/raw directory does not exist. Create it and add datasets.")
        sys.exit(1)

    ok_all_required = True
    for ds in DATASETS:
        ok = check_dataset(
            name=ds["name"],
            path=ds["path"],
            pattern=ds["pattern"],
            required_cols=ds["required_cols"],
            reader=ds["reader"],
            required=ds["required"],
        )
        if ds["required"]:
            ok_all_required = ok_all_required and ok

    if ok_all_required:
        print("\n✅ All required raw dataset checks passed.")
        sys.exit(0)
    else:
        print("\n⚠️ Some required datasets are missing or malformed. See messages above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
