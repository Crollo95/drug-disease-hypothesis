# src/ddh/pipeline.py

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .io_handlers import load_toy_data
from .graphs import build_ppi_graph
from .scoring import (
    build_drug_target_map,
    build_disease_gene_map,
    compute_overlap_table,
    compute_overlap_table_fast,
    compute_network_proximity,
    combine_overlap_and_proximity,
    attach_entity_names,
)


def run_toy_pipeline(
    data_dir: Path,
    alpha: float = 1.0,
    beta: float = 1.0,
) -> pd.DataFrame:
    """
    Run the full hypothesis pipeline on the toy dataset.

    Steps:
        1) Load toy data
        2) Build drug->genes and disease->genes maps
        3) Compute overlap-based scores
        4) Build PPI graph and compute network proximity
        5) Combine scores into a final ranking
        6) Attach human-readable names

    Parameters
    ----------
    data_dir : Path
        Directory containing the toy CSV files (e.g., 'data/toy').
    alpha : float
        Weight for normalized overlap term in combined score.
    beta : float
        Weight for proximity term in combined score.

    Returns
    -------
    pandas.DataFrame
        Ranked drug–disease pairs with scores and names.
    """
    data_dir = Path(data_dir)

    # 1. Load
    drugs, genes, diseases, dts, gds, ppis = load_toy_data(data_dir)

    # 2. Maps
    drug_to_genes = build_drug_target_map(dts)
    disease_to_genes = build_disease_gene_map(gds)

    # 3. Overlap
    overlap_df = compute_overlap_table_fast(dts, gds)

    # 4. Network proximity
    G = build_ppi_graph(ppis)
    prox_df = compute_network_proximity(G, drug_to_genes, disease_to_genes)

    # 5. Combine
    combined_df = combine_overlap_and_proximity(
        overlap_df=overlap_df,
        proximity_df=prox_df,
        alpha=alpha,
        beta=beta,
    )

    # 6. Attach names
    combined_df = attach_entity_names(combined_df, drugs, diseases)

    return combined_df



from .config import CsvFilesConfig
from .io_handlers import load_csv_data


from .config import CsvFilesConfig
from .io_handlers import load_csv_data
from .graphs import build_ppi_graph
from .scoring import (
    build_drug_target_map,
    build_disease_gene_map,
    compute_overlap_table,
    compute_overlap_table_fast,
    compute_network_proximity,
    combine_overlap_and_proximity,
    attach_entity_names,
)
import pandas as pd


def run_csv_pipeline(
    cfg: CsvFilesConfig,
    alpha: float = 1.0,
    beta: float = 1.0,
) -> pd.DataFrame:
    """
    Run the full hypothesis pipeline on generic CSV data.

    This is analogous to run_toy_pipeline, but uses user-provided CSV files
    configured via CsvFilesConfig.

    Parameters
    ----------
    cfg : CsvFilesConfig
        Configuration with file paths and column mappings.
    alpha : float
        Weight for normalized overlap term in combined score.
    beta : float
        Weight for proximity term in combined score.

    Returns
    -------
    pandas.DataFrame
        Ranked drug–disease pairs with scores and names.
    """
    # 1. Load CSV data
    drugs, genes, diseases, dts, gds, ppis = load_csv_data(cfg)

    # 2. Maps
    drug_to_genes = build_drug_target_map(dts)
    disease_to_genes = build_disease_gene_map(gds)

    # 3. Overlap
    overlap_df = compute_overlap_table_fast(dts, gds)

    # 4. Network proximity
    if ppis:
        # Real PPI: compute actual network proximity for all pairs
        G = build_ppi_graph(ppis)
        prox_df = compute_network_proximity(G, drug_to_genes, disease_to_genes)
    else:
        # No PPI: only define proximity for pairs that have overlap
        # (since proximity adds no extra info, avoid huge cartesian product)
        if overlap_df.empty:
            # Nothing to score
            return attach_entity_names(overlap_df, drugs, diseases)

        default_distance = float("inf")
        default_prox = 0.0

        prox_df = overlap_df[["drug_id", "disease_id"]].copy()
        prox_df["mean_distance"] = default_distance
        prox_df["proximity_score"] = default_prox

    # 5. Combine
    combined_df = combine_overlap_and_proximity(
        overlap_df=overlap_df,
        proximity_df=prox_df,
        alpha=alpha,
        beta=beta,
    )

    # 6. Attach names
    combined_df = attach_entity_names(combined_df, drugs, diseases)

    return combined_df


