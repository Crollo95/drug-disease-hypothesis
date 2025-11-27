#!/usr/bin/env python
"""
Build genes.csv from the union of gene IDs in:
  - data/real/drug_targets.csv (gene_id column)
  - data/real/gene_disease.csv (gene_id column)

Output:
  data/real/genes.csv with columns: gene_id, symbol
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd


def main() -> None:
    base = Path("data/real")

    dt_path = base / "drug_targets_filtered.csv"
    gd_path = base / "gene_disease_filtered.csv"
    out_path = base / "genes.csv"

    if not dt_path.exists():
        raise SystemExit(f"Missing {dt_path}")
    if not gd_path.exists():
        raise SystemExit(f"Missing {gd_path}")

    dt = pd.read_csv(dt_path)
    gd = pd.read_csv(gd_path)

    genes_from_dt = set(dt["gene_id"].astype(str).str.strip())
    genes_from_gd = set(gd["gene_id"].astype(str).str.strip())

    all_genes = sorted(genes_from_dt | genes_from_gd)

    df_genes = pd.DataFrame(
        {
            "gene_id": all_genes,
            "symbol": all_genes,  # for now, symbol = ID
        }
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df_genes.to_csv(out_path, index=False)

    print(f"genes.csv written to: {out_path.resolve()}")
    print(f"Total genes: {len(df_genes)}")


if __name__ == "__main__":
    main()

