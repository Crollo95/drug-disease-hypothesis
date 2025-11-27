#!/usr/bin/env python
"""
Precompute all-pairs shortest path distances between genes in the PPI graph.

Outputs:
    - data/real/gene_index.csv
        columns: gene_id, index

    - data/real/gene_distances.uint16.dat
        A binary file storing an N x N uint16 matrix (row-major),
        where N = number of genes in gene_index.csv.

        dist[i, j] = shortest path length between gene i and gene j in the PPI,
                     capped at max_uint16 (65535) for "no path / too far".

Usage (from project root):

    python scripts/precompute_gene_distance_matrix.py \
        --ppi-csv data/real/ppi.csv \
        --drug-targets-csv data/real/drug_targets_filtered.csv \
        --gene-disease-csv data/real/gene_disease_filtered.csv \
        --out-index data/real/gene_index.csv \
        --out-matrix data/real/gene_distances.uint16.dat
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Set, Dict

import numpy as np
import pandas as pd
import networkx as nx


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Precompute all-pairs gene–gene shortest path distances.",
    )

    parser.add_argument(
        "--ppi-csv",
        type=str,
        required=True,
        help="Path to ppi.csv (gene1_id,gene2_id,weight).",
    )
    parser.add_argument(
        "--drug-targets-csv",
        type=str,
        required=True,
        help="Path to drug_targets_filtered.csv (to include its genes).",
    )
    parser.add_argument(
        "--gene-disease-csv",
        type=str,
        required=True,
        help="Path to gene_disease_filtered.csv (to include its genes).",
    )
    parser.add_argument(
        "--out-index",
        type=str,
        required=True,
        help="Path to write gene_index.csv (gene_id,index).",
    )
    parser.add_argument(
        "--out-matrix",
        type=str,
        required=True,
        help="Path to write the uint16 distance matrix binary file.",
    )
    parser.add_argument(
        "--cutoff",
        type=int,
        default=None,
        help="Optional BFS cutoff (max distance). If None, no cutoff.",
    )

    return parser.parse_args(argv)


def build_gene_universe(
    ppi_path: Path,
    dt_path: Path,
    gd_path: Path,
) -> Set[str]:
    """Collect all gene IDs appearing in PPI, drug_targets, and gene_disease."""
    print(f"Loading PPI from {ppi_path} ...")
    ppi = pd.read_csv(ppi_path)
    genes_ppi = set(ppi["gene1_id"].astype(str).str.strip()) | set(
        ppi["gene2_id"].astype(str).str.strip()
    )

    print(f"Loading drug_targets from {dt_path} ...")
    dt = pd.read_csv(dt_path)
    genes_dt = set(dt["gene_id"].astype(str).str.strip())

    print(f"Loading gene_disease from {gd_path} ...")
    gd = pd.read_csv(gd_path)
    genes_gd = set(gd["gene_id"].astype(str).str.strip())

    all_genes = genes_ppi | genes_dt | genes_gd
    print(
        f"Genes: PPI={len(genes_ppi)}, "
        f"drug_targets={len(genes_dt)}, "
        f"gene_disease={len(genes_gd)}, "
        f"union={len(all_genes)}"
    )
    return all_genes


def build_ppi_graph(ppi_path: Path) -> nx.Graph:
    """Build an unweighted PPI graph from ppi.csv."""
    ppi = pd.read_csv(ppi_path)
    G = nx.Graph()
    # Ensure strings
    g1 = ppi["gene1_id"].astype(str).str.strip()
    g2 = ppi["gene2_id"].astype(str).str.strip()
    edges = list(zip(g1, g2))
    G.add_edges_from(edges)
    print(f"PPI graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


def main(argv: Optional[list[str]] = None) -> None:
    args = parse_args(argv)

    ppi_path = Path(args.ppi_csv)
    dt_path = Path(args.drug_targets_csv)
    gd_path = Path(args.gene_disease_csv)
    out_index_path = Path(args.out_index)
    out_matrix_path = Path(args.out_matrix)

    out_index_path.parent.mkdir(parents=True, exist_ok=True)
    out_matrix_path.parent.mkdir(parents=True, exist_ok=True)

    # 1. Build gene universe and index
    all_genes = sorted(build_gene_universe(ppi_path, dt_path, gd_path))
    N = len(all_genes)
    print(f"Total unique genes in universe: {N}")

    gene_to_idx: Dict[str, int] = {g: i for i, g in enumerate(all_genes)}

    # Save index mapping
    df_index = pd.DataFrame(
        {"gene_id": all_genes, "index": list(range(N))}
    )
    df_index.to_csv(out_index_path, index=False)
    print(f"Gene index written to: {out_index_path.resolve()}")

    # 2. Build PPI graph
    G = build_ppi_graph(ppi_path)

    # 3. Create memmap distance matrix
    max_uint16 = np.iinfo(np.uint16).max  # 65535
    print(f"Allocating distance matrix of shape ({N}, {N}) as uint16 ...")
    dist = np.memmap(
        out_matrix_path,
        dtype=np.uint16,
        mode="w+",
        shape=(N, N),
    )

    # Initialize with "infinite" distance
    dist[:] = max_uint16
    # Distance to self = 0
    for i in range(N):
        dist[i, i] = 0

    # 4. BFS from each gene that exists in the PPI graph
    genes_in_graph = [g for g in all_genes if g in G]
    print(f"Genes present in PPI graph: {len(genes_in_graph)}")

    cutoff = args.cutoff
    if cutoff is not None:
        print(f"Using BFS cutoff: {cutoff}")
    else:
        print("No BFS cutoff (full shortest paths).")

    for idx, g in enumerate(genes_in_graph, start=1):
        i = gene_to_idx[g]
        # Single-source shortest path lengths from g
        lengths = nx.single_source_shortest_path_length(G, g, cutoff=cutoff)
        for h, d in lengths.items():
            j = gene_to_idx.get(h)
            if j is None:
                continue
            if d <= max_uint16:
                dist[i, j] = d
        if idx % 100 == 0 or idx == len(genes_in_graph):
            print(f"Processed {idx}/{len(genes_in_graph)} sources")

    # Flush memmap to disk
    dist.flush()
    print(f"Distance matrix written to: {out_matrix_path.resolve()}")
    print("Done.")


if __name__ == "__main__":
    main()
