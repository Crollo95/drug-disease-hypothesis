# src/ddh/io_handlers.py

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import pandas as pd

from .data_models import (
    Drug,
    Gene,
    Disease,
    DrugTargetAssoc,
    GeneDiseaseAssoc,
    GeneGeneInteraction,
)


def load_toy_data(base_path: Path) -> Tuple[
    List[Drug],
    List[Gene],
    List[Disease],
    List[DrugTargetAssoc],
    List[GeneDiseaseAssoc],
    List[GeneGeneInteraction],
]:
    """
    Load toy data from CSV files in the given directory.

    Expected files:
        - toy_drugs.csv
        - toy_genes.csv
        - toy_diseases.csv
        - toy_drug_targets.csv
        - toy_gene_disease.csv
        - toy_ppi.csv

    Parameters
    ----------
    base_path : Path
        Directory that contains the toy CSV files.

    Returns
    -------
    drugs : list[Drug]
    genes : list[Gene]
    diseases : list[Disease]
    drug_targets : list[DrugTargetAssoc]
    gene_diseases : list[GeneDiseaseAssoc]
    ppis : list[GeneGeneInteraction]
    """
    base_path = Path(base_path)

    drugs_df = pd.read_csv(base_path / "toy_drugs.csv")
    genes_df = pd.read_csv(base_path / "toy_genes.csv")
    diseases_df = pd.read_csv(base_path / "toy_diseases.csv")
    dt_df = pd.read_csv(base_path / "toy_drug_targets.csv")
    gd_df = pd.read_csv(base_path / "toy_gene_disease.csv")
    ppi_df = pd.read_csv(base_path / "toy_ppi.csv")

    drugs = [
        Drug(id=row["drug_id"], name=row["drug_name"])
        for _, row in drugs_df.iterrows()
    ]

    genes = [
        Gene(id=row["gene_id"], symbol=row["symbol"])
        for _, row in genes_df.iterrows()
    ]

    diseases = [
        Disease(id=row["disease_id"], name=row["disease_name"])
        for _, row in diseases_df.iterrows()
    ]

    drug_targets = [
        DrugTargetAssoc(
            drug_id=row["drug_id"],
            gene_id=row["gene_id"],
            source="toy",
            score=float(row["score"]) if "score" in row and pd.notna(row["score"]) else None,
        )
        for _, row in dt_df.iterrows()
    ]

    gene_diseases = [
        GeneDiseaseAssoc(
            gene_id=row["gene_id"],
            disease_id=row["disease_id"],
            source="toy",
            score=float(row["score"]) if "score" in row and pd.notna(row["score"]) else None,
        )
        for _, row in gd_df.iterrows()
    ]

    ppis = [
        GeneGeneInteraction(
            gene1_id=row["gene1_id"],
            gene2_id=row["gene2_id"],
            weight=float(row["weight"]) if "weight" in row and pd.notna(row["weight"]) else None,
            source="toy",
        )
        for _, row in ppi_df.iterrows()
    ]

    return drugs, genes, diseases, drug_targets, gene_diseases, ppis


from .config import CsvFilesConfig


def load_csv_data(cfg: CsvFilesConfig):
    """
    Load data from CSV files according to the given configuration.

    Expected minimal schemas (defaults):

    - drugs_csv:
        columns: [drug_id, drug_name]

    - genes_csv:
        columns: [gene_id, symbol]

    - diseases_csv:
        columns: [disease_id, disease_name]

    - drug_targets_csv:
        columns: [drug_id, gene_id, (optional) score]

    - gene_disease_csv:
        columns: [gene_id, disease_id, (optional) score]

    - ppi_csv (optional):
        columns: [gene1_id, gene2_id, (optional) weight]

    Returns
    -------
    drugs : list[Drug]
    genes : list[Gene]
    diseases : list[Disease]
    drug_targets : list[DrugTargetAssoc]
    gene_diseases : list[GeneDiseaseAssoc]
    ppis : list[GeneGeneInteraction]
    """
    cfg = cfg.resolve_paths()

    drugs_df = pd.read_csv(cfg.drugs_csv)
    genes_df = pd.read_csv(cfg.genes_csv)
    diseases_df = pd.read_csv(cfg.diseases_csv)
    dt_df = pd.read_csv(cfg.drug_targets_csv)
    gd_df = pd.read_csv(cfg.gene_disease_csv)

    # Drugs
    drugs = [
        Drug(
            id=row[cfg.drug_id_col],
            name=row[cfg.drug_name_col],
        )
        for _, row in drugs_df.iterrows()
    ]

    # Genes
    genes = [
        Gene(
            id=row[cfg.gene_id_col],
            symbol=row[cfg.gene_symbol_col],
        )
        for _, row in genes_df.iterrows()
    ]

    # Diseases
    diseases = [
        Disease(
            id=row[cfg.disease_id_col],
            name=row[cfg.disease_name_col],
        )
        for _, row in diseases_df.iterrows()
    ]

    # Drug-target associations
    drug_targets = []
    for _, row in dt_df.iterrows():
        score_val = None
        if cfg.dt_score_col is not None and cfg.dt_score_col in dt_df.columns:
            val = row[cfg.dt_score_col]
            if pd.notna(val):
                score_val = float(val)
        drug_targets.append(
            DrugTargetAssoc(
                drug_id=row[cfg.dt_drug_id_col],
                gene_id=row[cfg.dt_gene_id_col],
                source=str(cfg.drug_targets_csv),
                score=score_val,
            )
        )

    # Gene-disease associations
    gene_diseases = []
    for _, row in gd_df.iterrows():
        score_val = None
        if cfg.gd_score_col is not None and cfg.gd_score_col in gd_df.columns:
            val = row[cfg.gd_score_col]
            if pd.notna(val):
                score_val = float(val)
        gene_diseases.append(
            GeneDiseaseAssoc(
                gene_id=row[cfg.gd_gene_id_col],
                disease_id=row[cfg.gd_disease_id_col],
                source=str(cfg.gene_disease_csv),
                score=score_val,
            )
        )

    # PPI / gene-gene interactions
    ppis: List[GeneGeneInteraction] = []
    if cfg.ppi_csv is not None and cfg.ppi_csv.exists():
        ppi_df = pd.read_csv(cfg.ppi_csv)
        for _, row in ppi_df.iterrows():
            weight_val = None
            if cfg.ppi_weight_col is not None and cfg.ppi_weight_col in ppi_df.columns:
                val = row[cfg.ppi_weight_col]
                if pd.notna(val):
                    weight_val = float(val)
            ppis.append(
                GeneGeneInteraction(
                    gene1_id=row[cfg.ppi_gene1_col],
                    gene2_id=row[cfg.ppi_gene2_col],
                    weight=weight_val,
                    source=str(cfg.ppi_csv),
                )
            )

    return drugs, genes, diseases, drug_targets, gene_diseases, ppis

