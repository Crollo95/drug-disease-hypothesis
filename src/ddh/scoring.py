# src/ddh/scoring.py

from __future__ import annotations

from typing import Dict, List, Set, Tuple

import pandas as pd
import numpy as np
import networkx as nx

from .data_models import DrugTargetAssoc, GeneDiseaseAssoc, Drug, Disease


def build_drug_target_map(
    drug_targets: List[DrugTargetAssoc],
) -> Dict[str, Set[str]]:
    """
    Build a mapping from drug_id -> set of target gene_ids.

    Parameters
    ----------
    drug_targets : list[DrugTargetAssoc]
        Drug–target association records.

    Returns
    -------
    dict
        Mapping: drug_id -> set of gene_ids.
    """
    mapping: Dict[str, Set[str]] = {}
    for assoc in drug_targets:
        if assoc.drug_id not in mapping:
            mapping[assoc.drug_id] = set()
        mapping[assoc.drug_id].add(assoc.gene_id)
    return mapping


def build_disease_gene_map(
    gene_diseases: List[GeneDiseaseAssoc],
) -> Dict[str, Set[str]]:
    """
    Build a mapping from disease_id -> set of associated gene_ids.

    Parameters
    ----------
    gene_diseases : list[GeneDiseaseAssoc]
        Gene–disease association records.

    Returns
    -------
    dict
        Mapping: disease_id -> set of gene_ids.
    """
    mapping: Dict[str, Set[str]] = {}
    for assoc in gene_diseases:
        if assoc.disease_id not in mapping:
            mapping[assoc.disease_id] = set()
        mapping[assoc.disease_id].add(assoc.gene_id)
    return mapping


def compute_overlap_table(
    drug_to_genes: Dict[str, Set[str]],
    disease_to_genes: Dict[str, Set[str]],
) -> pd.DataFrame:
    """
    Compute a simple overlap-based score for all drug–disease pairs.

    For each (drug_id, disease_id) pair, we compute:
        - n_overlap: number of shared genes between drug targets and disease genes
        - overlapping_genes: semicolon-separated list of overlapping gene IDs
        - jaccard: Jaccard index = |intersection| / |union|

    Only pairs with at least one overlapping gene are returned.

    Parameters
    ----------
    drug_to_genes : dict
        Mapping of drug_id -> set of target gene_ids.
    disease_to_genes : dict
        Mapping of disease_id -> set of associated gene_ids.

    Returns
    -------
    pandas.DataFrame
        Columns:
            - drug_id
            - disease_id
            - n_overlap
            - overlapping_genes
            - jaccard
    """
    rows: List[Dict[str, object]] = []

    for drug_id, drug_genes in drug_to_genes.items():
        for disease_id, disease_genes in disease_to_genes.items():
            inter = drug_genes & disease_genes
            if not inter:
                continue

            union = drug_genes | disease_genes
            n_overlap = len(inter)
            jaccard = n_overlap / len(union) if union else 0.0

            rows.append(
                {
                    "drug_id": drug_id,
                    "disease_id": disease_id,
                    "n_overlap": n_overlap,
                    "overlapping_genes": ";".join(sorted(inter)),
                    "jaccard": jaccard,
                }
            )

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(
            by=["n_overlap", "jaccard", "drug_id", "disease_id"],
            ascending=[False, False, True, True],
        ).reset_index(drop=True)
    return df


def compute_overlap_table_fast(
    drug_targets: List[DrugTargetAssoc],
    gene_diseases: List[GeneDiseaseAssoc],
) -> pd.DataFrame:
    """
    Compute overlap-based scores using a vectorized join instead of
    nested loops over all drug × disease pairs.

    Parameters
    ----------
    drug_targets : list[DrugTargetAssoc]
        Drug–target associations (drug_id, gene_id).
    gene_diseases : list[GeneDiseaseAssoc]
        Gene–disease associations (gene_id, disease_id).

    Returns
    -------
    pandas.DataFrame
        Columns:
            - drug_id
            - disease_id
            - n_overlap
            - overlapping_genes
            - jaccard
    """
    if not drug_targets or not gene_diseases:
        return pd.DataFrame(
            columns=[
                "drug_id",
                "disease_id",
                "n_overlap",
                "overlapping_genes",
                "jaccard",
            ]
        )

    # Build edge tables
    dt_df = pd.DataFrame(
        [(a.drug_id, a.gene_id) for a in drug_targets],
        columns=["drug_id", "gene_id"],
    )
    gd_df = pd.DataFrame(
        [(a.gene_id, a.disease_id) for a in gene_diseases],
        columns=["gene_id", "disease_id"],
    )

    # Drop obvious NAs
    dt_df = dt_df.dropna(subset=["drug_id", "gene_id"])
    gd_df = gd_df.dropna(subset=["gene_id", "disease_id"])

    # How many unique genes per drug / disease? (for Jaccard)
    n_genes_per_drug = (
        dt_df.groupby("drug_id")["gene_id"].nunique().rename("n_drug_genes")
    )
    n_genes_per_disease = (
        gd_df.groupby("disease_id")["gene_id"].nunique().rename("n_disease_genes")
    )

    # Join on gene_id → all overlapping triples (drug, disease, gene)
    merged = dt_df.merge(gd_df, on="gene_id", how="inner")

    if merged.empty:
        return pd.DataFrame(
            columns=[
                "drug_id",
                "disease_id",
                "n_overlap",
                "overlapping_genes",
                "jaccard",
            ]
        )

    # Aggregate overlaps
    grouped = merged.groupby(["drug_id", "disease_id"])["gene_id"]

    overlap_df = grouped.agg(
        n_overlap="nunique",
        overlapping_genes=lambda x: ";".join(sorted(set(x))),
    ).reset_index()

    # Attach gene counts for Jaccard
    overlap_df = overlap_df.merge(
        n_genes_per_drug, on="drug_id", how="left"
    ).merge(
        n_genes_per_disease, on="disease_id", how="left"
    )

    # Compute Jaccard
    union_size = (
        overlap_df["n_drug_genes"] + overlap_df["n_disease_genes"] - overlap_df["n_overlap"]
    )
    overlap_df["jaccard"] = overlap_df["n_overlap"] / union_size.replace({0: pd.NA})

    # Clean up
    overlap_df = overlap_df.drop(columns=["n_drug_genes", "n_disease_genes"])
    overlap_df = overlap_df.sort_values(
        by=["n_overlap", "jaccard", "drug_id", "disease_id"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)

    return overlap_df


def attach_entity_names(
    df: pd.DataFrame,
    drugs: List[Drug],
    diseases: List[Disease],
) -> pd.DataFrame:
    """
    Add human-readable drug and disease names to a results DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        Must contain columns 'drug_id' and 'disease_id'.
    drugs : list[Drug]
        Drug objects with id and name.
    diseases : list[Disease]
        Disease objects with id and name.

    Returns
    -------
    pandas.DataFrame
        Same as input df, with two extra columns:
            - drug_name
            - disease_name
    """
    drug_map = {d.id: d.name for d in drugs}
    disease_map = {d.id: d.name for d in diseases}

    df = df.copy()
    df["drug_name"] = df["drug_id"].map(drug_map)
    df["disease_name"] = df["disease_id"].map(disease_map)
    return df


def compute_network_proximity(
    G: nx.Graph,
    drug_to_genes: Dict[str, Set[str]],
    disease_to_genes: Dict[str, Set[str]],
    default_distance: float = 5.0,
) -> pd.DataFrame:
    """
    Compute average shortest-path distance between drug targets and disease genes.

    Smaller distances mean closer proximity in the network.

    For each (drug_id, disease_id) pair, we compute:
        - mean_distance: mean shortest-path length between all pairs
                         (drug_target, disease_gene) that are connected.
        - proximity_score: 1 / (1 + mean_distance)
                           (so larger is better; in (0, 1]).

    If no paths are found, mean_distance is set to `default_distance`
    and proximity_score is computed from that.

    Parameters
    ----------
    G : networkx.Graph
        PPI graph with gene_ids as nodes.
    drug_to_genes : dict
        Mapping drug_id -> set of gene_ids (targets).
    disease_to_genes : dict
        Mapping disease_id -> set of gene_ids (disease-associated).
    default_distance : float, optional
        Default distance when no path exists.

    Returns
    -------
    pandas.DataFrame
        Columns:
            - drug_id
            - disease_id
            - mean_distance
            - proximity_score
    """
    rows: List[Dict[str, object]] = []

    for drug_id, drug_genes in drug_to_genes.items():
        for disease_id, dis_genes in disease_to_genes.items():
            distances = []

            for t in drug_genes:
                for g in dis_genes:
                    if t in G and g in G:
                        try:
                            d = nx.shortest_path_length(G, t, g)
                            distances.append(d)
                        except nx.NetworkXNoPath:
                            continue

            if distances:
                mean_distance = float(np.mean(distances))
            else:
                mean_distance = float(default_distance)

            proximity_score = 1.0 / (1.0 + mean_distance)

            rows.append(
                {
                    "drug_id": drug_id,
                    "disease_id": disease_id,
                    "mean_distance": mean_distance,
                    "proximity_score": proximity_score,
                }
            )

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(
            by=["proximity_score", "drug_id", "disease_id"],
            ascending=[False, True, True],
        ).reset_index(drop=True)

    return df


def combine_overlap_and_proximity(
    overlap_df: pd.DataFrame,
    proximity_df: pd.DataFrame,
    alpha: float = 1.0,
    beta: float = 1.0,
) -> pd.DataFrame:
    """
    Combine overlap and network proximity into a single score.

    Combined score = alpha * normalized_overlap + beta * proximity_score

    where normalized_overlap = n_overlap / max(n_overlap) (if available).

    Parameters
    ----------
    overlap_df : pandas.DataFrame
        Result of compute_overlap_table (can be empty).
    proximity_df : pandas.DataFrame
        Result of compute_network_proximity.
    alpha : float
        Weight for normalized overlap term.
    beta : float
        Weight for proximity_score term.

    Returns
    -------
    pandas.DataFrame
        Joined table with:
            - drug_id
            - disease_id
            - n_overlap (0 if no overlap)
            - overlapping_genes (if available)
            - jaccard (if available)
            - mean_distance
            - proximity_score
            - combined_score
    """
    # Outer join: keep pairs that appear in either table
    df = pd.merge(
        proximity_df,
        overlap_df,
        on=["drug_id", "disease_id"],
        how="outer",
        suffixes=("", "_overlap"),
    )

    # If n_overlap column is missing entirely, create it as zeros
    if "n_overlap" not in df.columns:
        df["n_overlap"] = 0
    else:
        # Replace NaN with 0 for pairs that have proximity but no overlap
        df["n_overlap"] = df["n_overlap"].fillna(0)

    # Normalize overlap count to [0, 1] range
    max_overlap = df["n_overlap"].max() if not df["n_overlap"].isna().all() else 0
    if max_overlap > 0:
        df["norm_overlap"] = df["n_overlap"] / max_overlap
    else:
        df["norm_overlap"] = 0.0

    # Ensure proximity_score exists and has no NaNs
    if "proximity_score" not in df.columns:
        df["proximity_score"] = 0.0
    df["proximity_score"] = df["proximity_score"].fillna(0.0)

    # Combined score: if norm_overlap or proximity_score are numeric, this will be numeric
    df["combined_score"] = alpha * df["norm_overlap"] + beta * df["proximity_score"]

    # Sort: best combined score first
    df = df.sort_values(
        by=["combined_score", "proximity_score", "n_overlap"],
        ascending=[False, False, False],
    ).reset_index(drop=True)

    return df



