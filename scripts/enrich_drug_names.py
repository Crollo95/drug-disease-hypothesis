#!/usr/bin/env python
"""
Enrich drug names in data/real/drugs.csv (and optionally drug_targets.csv)
using data/real/drug_lookup.csv.

Outputs:
    Overwrites:
        data/real/drugs.csv
        data/real/drug_targets.csv (if present and has drug_name)
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd


def main() -> None:
    base = Path("data/real")
    drugs_path = base / "drugs.csv"
    lookup_path = base / "drug_lookup.csv"
    dt_path = base / "drug_targets.csv"

    if not drugs_path.exists():
        raise SystemExit(f"Missing {drugs_path}")
    if not lookup_path.exists():
        raise SystemExit(f"Missing {lookup_path}")

    print(f"Loading drugs from {drugs_path} ...")
    dr = pd.read_csv(drugs_path)
    dr["drug_id"] = dr["drug_id"].astype(str).str.strip()

    print(f"Loading lookup from {lookup_path} ...")
    lookup = pd.read_csv(lookup_path)
    lookup["drug_id"] = lookup["drug_id"].astype(str).str.strip()
    lookup["drug_name"] = lookup["drug_name"].astype(str).str.strip()

    print("Merging drug names...")
    dr = dr.merge(lookup, on="drug_id", how="left", suffixes=("", "_lookup"))

    # Prefer lookup name if available
    dr["drug_name"] = dr["drug_name_lookup"].combine_first(dr.get("drug_name"))
    dr = dr.drop(columns=[c for c in dr.columns if c.endswith("_lookup")])

    dr.to_csv(drugs_path, index=False)
    print(f"Updated drugs.csv written to: {drugs_path.resolve()}")

    # Optionally update drug_targets.csv if it exists and has drug_name
    if dt_path.exists():
        print(f"Updating drug_targets.csv at {dt_path} ...")
        dt = pd.read_csv(dt_path)
        dt["drug_id"] = dt["drug_id"].astype(str).str.strip()
        dt = dt.merge(lookup, on="drug_id", how="left", suffixes=("", "_lookup"))

        if "drug_name" in dt.columns:
            dt["drug_name"] = dt["drug_name_lookup"].combine_first(dt["drug_name"])
        else:
            dt["drug_name"] = dt["drug_name_lookup"]

        dt = dt.drop(columns=[c for c in dt.columns if c.endswith("_lookup")])
        dt.to_csv(dt_path, index=False)
        print("Updated drug_targets.csv written.")
    else:
        print("No drug_targets.csv found; skipping.")


if __name__ == "__main__":
    main()
