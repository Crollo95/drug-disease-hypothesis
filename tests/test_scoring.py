# tests/test_scoring.py

from pathlib import Path

import pandas as pd

from ddh.io_handlers import load_toy_data
from ddh.scoring import (
    build_drug_target_map,
    build_disease_gene_map,
    compute_overlap_table,
    compute_network_proximity,
    combine_overlap_and_proximity,
)


def _load_toy():
    data_dir = Path("data/toy")
    return load_toy_data(data_dir)


def test_build_maps_non_empty():
    drugs, genes, diseases, dts, gds, ppis = _load_toy()

    drug_to_genes = build_drug_target_map(dts)
    disease_to_genes = build_disease_gene_map(gds)

    assert len(drug_to_genes) > 0
    assert len(disease_to_genes) > 0

    # Check that a known drug has expected targets
    assert "D1" in drug_to_genes
    assert "G1" in drug_to_genes["D1"]
    assert "G2" in drug_to_genes["D1"]

    # Check that a known disease has expected genes
    assert "DIS1" in disease_to_genes
    assert "G1" in disease_to_genes["DIS1"]
    assert "G2" in disease_to_genes["DIS1"]


def test_compute_overlap_table_expected_pairs():
    drugs, genes, diseases, dts, gds, ppis = _load_toy()

    drug_to_genes = build_drug_target_map(dts)
    disease_to_genes = build_disease_gene_map(gds)

    df = compute_overlap_table(drug_to_genes, disease_to_genes)

    # There should be at least one overlap-based hypothesis
    assert not df.empty

    # D1 targets G1,G2 and DIS1 involves G1,G2 -> 2-gene overlap
    row = df[(df["drug_id"] == "D1") & (df["disease_id"] == "DIS1")]
    assert len(row) == 1

    row = row.iloc[0]
    assert row["n_overlap"] == 2
    overlapping = set(row["overlapping_genes"].split(";"))
    assert overlapping == {"G1", "G2"}

    # Jaccard should be > 0
    assert row["jaccard"] > 0.0


def test_compute_network_and_combined_scores():
    from ddh.graphs import build_ppi_graph

    drugs, genes, diseases, dts, gds, ppis = _load_toy()

    drug_to_genes = build_drug_target_map(dts)
    disease_to_genes = build_disease_gene_map(gds)

    # Build graph and proximity table
    G = build_ppi_graph(ppis)
    prox_df = compute_network_proximity(G, drug_to_genes, disease_to_genes)

    assert not prox_df.empty
    assert {"drug_id", "disease_id", "mean_distance", "proximity_score"}.issubset(
        set(prox_df.columns)
    )

    # Combined score table
    overlap_df = compute_overlap_table(drug_to_genes, disease_to_genes)
    combined_df = combine_overlap_and_proximity(
        overlap_df=overlap_df,
        proximity_df=prox_df,
        alpha=1.0,
        beta=1.0,
    )

    assert not combined_df.empty
    assert "combined_score" in combined_df.columns

    # All combined scores should be finite numbers
    assert pd.api.types.is_numeric_dtype(combined_df["combined_score"])
    assert combined_df["combined_score"].notna().all()

