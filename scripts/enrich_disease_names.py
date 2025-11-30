#!/usr/bin/env python
"""
Enrich disease names in data/real/diseases.csv (and optionally gene_disease.csv)
using data/real/disease_lookup.csv.
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd


def main() -> None:
    base = Path("data/real")
    diseases_path = base / "diseases.csv"
    lookup_path = base / "disease_lookup.csv"
    gene_dis_path = base / "gene_disease.csv"

    if not diseases_path.exists():
        raise SystemExit(f"Missing {diseases_path}")
    if not lookup_path.exists():
        raise SystemExit(f"Missing {lookup_path}")

    # --- Update diseases.csv ---
    print(f"Loading diseases from {diseases_path} ...")
    dis = pd.read_csv(diseases_path)
    dis["disease_id"] = dis["disease_id"].astype(str).str.strip()

    print(f"Loading lookup from {lookup_path} ...")
    lookup = pd.read_csv(lookup_path)
    lookup["disease_id"] = lookup["disease_id"].astype(str).str.strip()
    lookup["disease_name"] = lookup["disease_name"].astype(str).str.strip()

    print("Merging disease names into diseases.csv ...")
    dis = dis.merge(lookup, on="disease_id", how="left", suffixes=("", "_lookup"))

    # If both original and lookup name exist, prefer lookup
    if "disease_name_lookup" in dis.columns and "disease_name" in dis.columns:
        dis["disease_name"] = dis["disease_name_lookup"].combine_first(dis["disease_name"])
        dis = dis.drop(columns=["disease_name_lookup"])
    elif "disease_name_lookup" in dis.columns:
        # Only lookup name exists
        dis = dis.rename(columns={"disease_name_lookup": "disease_name"})

    dis.to_csv(diseases_path, index=False)
    print(f"Updated diseases.csv written to: {diseases_path.resolve()}")

    # --- Optionally update gene_disease.csv ---
    if gene_dis_path.exists():
        print(f"Updating gene_disease.csv at {gene_dis_path} ...")
        gd = pd.read_csv(gene_dis_path)
        gd["disease_id"] = gd["disease_id"].astype(str).str.strip()

        gd = gd.merge(lookup, on="disease_id", how="left", suffixes=("", "_lookup"))

        # Cases:
        # - If gene_disease already had a disease_name column, prefer lookup where available
        # - If it did not, just use lookup's disease_name
        if "disease_name_lookup" in gd.columns and "disease_name" in gd.columns:
            gd["disease_name"] = gd["disease_name_lookup"].combine_first(gd["disease_name"])
            gd = gd.drop(columns=["disease_name_lookup"])
        elif "disease_name_lookup" in gd.columns:
            gd = gd.rename(columns={"disease_name_lookup": "disease_name"})
        # else: no name from lookup; nothing to do

        gd.to_csv(gene_dis_path, index=False)
        print("Updated gene_disease.csv written.")
    else:
        print("No gene_disease.csv found; skipping.")


if __name__ == "__main__":
    main()
