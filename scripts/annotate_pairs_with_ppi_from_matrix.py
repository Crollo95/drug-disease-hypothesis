#!/usr/bin/env python
"""
Annotate a list of drug–disease pairs with PPI-based proximity using a
precomputed gene–gene distance matrix.

Inputs:
    --pairs-csv
        CSV with at least columns:
            - drug_id
            - disease_id
        (Optionally may contain n_overlap, jaccard, etc.)

    --drug-targets-csv
        CSV with columns:
            - drug_id
            - gene_id

    --gene-disease-csv
        CSV with columns:
            - gene_id
            - disease_id

    --gene-index
        CSV produced by precompute_gene_distance_matrix.py, with columns:
            - gene_id
            - index

    --dist-matrix
        Binary file produced by precompute_gene_distance_matrix.py,
        storing an N x N uint16 matrix (row-major) with shortest path
        distances between genes.

Outputs:
    --output
        CSV = input pairs + new columns:
            - mean_distance
            - proximity_score
            - combined_score   (if n_overlap column present)

You can also limit the number of pairs processed via --max-pairs.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Dict, List, Iterable, Tuple

import numpy as np
import pandas as pd


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate drug–disease pairs with PPI proximity using a precomputed distance matrix.",
    )

    parser.add_argument(
        "--pairs-csv",
        type=str,
        required=True,
        help="CSV of drug–disease pairs (must have drug_id,disease_id).",
    )
    parser.add_argument(
        "--drug-targets-csv",
        type=str,
        required=True,
        help="CSV of drug–target associations (drug_id,gene_id).",
    )
    parser.add_argument(
        "--gene-disease-csv",
        type=str,
        required=True,
        help="CSV of gene–disease associations (gene_id,disease_id).",
    )
    parser.add_argument(
        "--gene-index",
        type=str,
        required=True,
        help="CSV mapping gene_id to index (from precompute_gene_distance_matrix.py).",
    )
    parser.add_argument(
        "--dist-matrix",
        type=str,
        required=True,
        help="Binary file with N x N uint16 distance matrix.",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to write annotated pairs CSV.",
    )

    parser.add_argument(
        "--max-pairs",
        type=int,
        default=None,
        help="Optional: max number of pairs to process (take first N rows).",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10_000,
        help="Process pairs in chunks of this size.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=1.0,
        help="Weight for normalized overlap term in combined_score.",
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=1.0,
        help="Weight for normalized proximity term in combined_score.",
    )

    return parser.parse_args(argv)


def build_gene_maps(
    dt_path: Path,
    gd_path: Path,
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Build drug_id -> list[gene_id] and disease_id -> list[gene_id] maps
    from drug_targets and gene_disease CSVs.
    """
    print(f"Loading drug_targets from {dt_path} ...")
    dt = pd.read_csv(dt_path, usecols=["drug_id", "gene_id"])
    dt["drug_id"] = dt["drug_id"].astype(str).str.strip()
    dt["gene_id"] = dt["gene_id"].astype(str).str.strip()

    print(f"Loading gene_disease from {gd_path} ...")
    gd = pd.read_csv(gd_path, usecols=["gene_id", "disease_id"])
    gd["gene_id"] = gd["gene_id"].astype(str).str.strip()
    gd["disease_id"] = gd["disease_id"].astype(str).str.strip()

    drug_to_genes: Dict[str, List[str]] = (
        dt.groupby("drug_id")["gene_id"].apply(list).to_dict()
    )
    disease_to_genes: Dict[str, List[str]] = (
        gd.groupby("disease_id")["gene_id"].apply(list).to_dict()
    )

    print(f"Drugs with targets: {len(drug_to_genes)}")
    print(f"Diseases with genes: {len(disease_to_genes)}")
    return drug_to_genes, disease_to_genes


def load_distance_matrix(
    gene_index_path: Path,
    dist_matrix_path: Path,
):
    """
    Load gene_index and the memmapped distance matrix.
    Returns:
        gene_to_idx: dict[gene_id -> int]
        dist: np.memmap of shape (N, N), dtype=uint16
        max_val: sentinel for "no path"
    """
    print(f"Loading gene_index from {gene_index_path} ...")
    idx_df = pd.read_csv(gene_index_path)
    idx_df["gene_id"] = idx_df["gene_id"].astype(str).str.strip()
    gene_to_idx = dict(zip(idx_df["gene_id"], idx_df["index"]))
    N = len(idx_df)
    print(f"Gene index loaded: {N} genes.")

    print(f"Opening distance matrix from {dist_matrix_path} ...")
    dist = np.memmap(
        dist_matrix_path,
        dtype=np.uint16,
        mode="r",
        shape=(N, N),
    )
    max_val = np.iinfo(np.uint16).max
    return gene_to_idx, dist, max_val


def mean_distance_for_pair(
    drug_id: str,
    disease_id: str,
    drug_to_genes: Dict[str, List[str]],
    disease_to_genes: Dict[str, List[str]],
    gene_to_idx: Dict[str, int],
    dist: np.memmap,
    max_val: int,
) -> float:
    """
    Compute mean shortest-path distance between all (drug_gene, disease_gene)
    pairs using the precomputed distance matrix.

    Returns inf if no valid distances.
    """
    genes_d = drug_to_genes.get(drug_id)
    genes_s = disease_to_genes.get(disease_id)
    if not genes_d or not genes_s:
        return float("inf")

    idx_d = [gene_to_idx[g] for g in set(genes_d) if g in gene_to_idx]
    idx_s = [gene_to_idx[g] for g in set(genes_s) if g in gene_to_idx]

    if not idx_d or not idx_s:
        return float("inf")

    dists: List[int] = []
    for i in idx_d:
        row = dist[i]
        for j in idx_s:
            d = int(row[j])
            if d < max_val:  # skip "no path"
                dists.append(d)

    if not dists:
        return float("inf")

    return float(sum(dists) / len(dists))


def annotate_pairs(
    pairs_df: pd.DataFrame,
    drug_to_genes: Dict[str, List[str]],
    disease_to_genes: Dict[str, List[str]],
    gene_to_idx: Dict[str, int],
    dist: np.memmap,
    max_val: int,
    chunk_size: int = 10_000,
) -> pd.DataFrame:
    """
    For each row in pairs_df (must have drug_id,disease_id), compute
    mean_distance and proximity_score and return a new DataFrame.
    """
    pairs_df = pairs_df.copy()
    pairs_df["drug_id"] = pairs_df["drug_id"].astype(str).str.strip()
    pairs_df["disease_id"] = pairs_df["disease_id"].astype(str).str.strip()

    mean_distances: List[float] = []

    n_rows = len(pairs_df)
    print(f"Annotating {n_rows} pairs with PPI proximity ...")

    # Process in chunks to allow progress reporting
    for start in range(0, n_rows, chunk_size):
        end = min(start + chunk_size, n_rows)
        chunk = pairs_df.iloc[start:end]

        for _, row in chunk.iterrows():
            d_id = row["drug_id"]
            dis_id = row["disease_id"]
            md = mean_distance_for_pair(
                d_id,
                dis_id,
                drug_to_genes,
                disease_to_genes,
                gene_to_idx,
                dist,
                max_val,
            )
            mean_distances.append(md)

        print(f"Processed rows {start}..{end-1} / {n_rows}")

    pairs_df["mean_distance"] = mean_distances

    # Convert mean_distance to a proximity_score in [0,1]
    # Simple transform: proximity = 1 / (1 + distance), with 0 for inf
    prox = []
    for d in pairs_df["mean_distance"]:
        if np.isinf(d):
            prox.append(0.0)
        else:
            prox.append(1.0 / (1.0 + d))
    pairs_df["proximity_score"] = prox

    return pairs_df


def add_combined_score(
    df: pd.DataFrame,
    alpha: float,
    beta: float,
) -> pd.DataFrame:
    """
    Optionally add combined_score as a weighted sum of normalized
    overlap (via n_overlap) and normalized proximity_score.

    If n_overlap is missing, we fall back to proximity only.
    """
    df = df.copy()

    has_overlap = "n_overlap" in df.columns

    # Normalize proximity
    prox = df["proximity_score"].astype(float)
    prox_min = prox.min()
    prox_max = prox.max()
    if prox_max > prox_min:
        prox_norm = (prox - prox_min) / (prox_max - prox_min)
    else:
        prox_norm = prox * 0.0  # all zeros

    if has_overlap:
        ov = df["n_overlap"].astype(float)
        ov_min = ov.min()
        ov_max = ov.max()
        if ov_max > ov_min:
            ov_norm = (ov - ov_min) / (ov_max - ov_min)
        else:
            ov_norm = ov * 0.0

        df["combined_score"] = alpha * ov_norm + beta * prox_norm
    else:
        # No overlap info: combined_score = proximity only
        df["combined_score"] = beta * prox_norm

    return df


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)

    pairs_path = Path(args.pairs_csv)
    dt_path = Path(args.drug_targets_csv)
    gd_path = Path(args.gene_disease_csv)
    gene_index_path = Path(args.gene_index)
    dist_matrix_path = Path(args.dist_matrix)
    output_path = Path(args.output)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Load mappings
    drug_to_genes, disease_to_genes = build_gene_maps(dt_path, gd_path)
    gene_to_idx, dist, max_val = load_distance_matrix(gene_index_path, dist_matrix_path)

    # Load pairs
    print(f"Loading pairs from {pairs_path} ...")
    pairs_df = pd.read_csv(pairs_path)

    if "drug_id" not in pairs_df.columns or "disease_id" not in pairs_df.columns:
        raise SystemExit("pairs_csv must have 'drug_id' and 'disease_id' columns.")

    if args.max_pairs is not None:
        pairs_df = pairs_df.iloc[: args.max_pairs].copy()
        print(f"Restricting to first {len(pairs_df)} pairs due to --max-pairs.")

    annotated_df = annotate_pairs(
        pairs_df,
        drug_to_genes,
        disease_to_genes,
        gene_to_idx,
        dist,
        max_val,
        chunk_size=args.chunk_size,
    )

    # Add combined_score (optional, based on n_overlap if present)
    annotated_df = add_combined_score(
        annotated_df,
        alpha=args.alpha,
        beta=args.beta,
    )

    annotated_df.to_csv(output_path, index=False)
    print(f"Annotated pairs written to: {output_path.resolve()}")


if __name__ == "__main__":
    main()
