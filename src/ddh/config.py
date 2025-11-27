# src/ddh/config.py

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class CsvFilesConfig:
    """
    Configuration for loading drug, gene, disease, and interaction data
    from CSV files with a simple schema.

    The idea:
        - You provide file paths
        - We assume standard column names, or you can override them
    """

    # File paths
    drugs_csv: Path
    genes_csv: Path
    diseases_csv: Path
    drug_targets_csv: Path
    gene_disease_csv: Path
    ppi_csv: Optional[Path] = None

    # Column names for each file (with sensible defaults)

    # Drugs: columns [drug_id_col, drug_name_col]
    drug_id_col: str = "drug_id"
    drug_name_col: str = "drug_name"

    # Genes: columns [gene_id_col, gene_symbol_col]
    gene_id_col: str = "gene_id"
    gene_symbol_col: str = "symbol"

    # Diseases: columns [disease_id_col, disease_name_col]
    disease_id_col: str = "disease_id"
    disease_name_col: str = "disease_name"

    # Drug-target associations
    dt_drug_id_col: str = "drug_id"
    dt_gene_id_col: str = "gene_id"
    dt_score_col: Optional[str] = "score"

    # Gene-disease associations
    gd_gene_id_col: str = "gene_id"
    gd_disease_id_col: str = "disease_id"
    gd_score_col: Optional[str] = "score"

    # PPI / gene-gene interactions
    ppi_gene1_col: str = "gene1_id"
    ppi_gene2_col: str = "gene2_id"
    ppi_weight_col: Optional[str] = "weight"

    def resolve_paths(self, base_dir: Path | None = None) -> "CsvFilesConfig":
        """
        Return a copy of this config with all paths resolved (absolute).
        If base_dir is provided, relative paths are interpreted relative to it.
        """
        base_dir = Path(base_dir) if base_dir is not None else Path(".")
        return CsvFilesConfig(
            drugs_csv=(base_dir / self.drugs_csv).resolve(),
            genes_csv=(base_dir / self.genes_csv).resolve(),
            diseases_csv=(base_dir / self.diseases_csv).resolve(),
            drug_targets_csv=(base_dir / self.drug_targets_csv).resolve(),
            gene_disease_csv=(base_dir / self.gene_disease_csv).resolve(),
            ppi_csv=(base_dir / self.ppi_csv).resolve() if self.ppi_csv is not None else None,
            drug_id_col=self.drug_id_col,
            drug_name_col=self.drug_name_col,
            gene_id_col=self.gene_id_col,
            gene_symbol_col=self.gene_symbol_col,
            disease_id_col=self.disease_id_col,
            disease_name_col=self.disease_name_col,
            dt_drug_id_col=self.dt_drug_id_col,
            dt_gene_id_col=self.dt_gene_id_col,
            dt_score_col=self.dt_score_col,
            gd_gene_id_col=self.gd_gene_id_col,
            gd_disease_id_col=self.gd_disease_id_col,
            gd_score_col=self.gd_score_col,
            ppi_gene1_col=self.ppi_gene1_col,
            ppi_gene2_col=self.ppi_gene2_col,
            ppi_weight_col=self.ppi_weight_col,
        )

