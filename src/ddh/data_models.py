# src/ddh/data_models.py

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Drug:
    """
    Representation of a drug or compound.

    Attributes
    ----------
    id : str
        Stable identifier (e.g. 'DRUGBANK:DB0001', 'CHEMBL:CHEMBL25', or internal 'D123').
    name : str
        Human-readable drug name.
    """
    id: str
    name: str


@dataclass(frozen=True)
class Gene:
    """
    Representation of a gene.

    Attributes
    ----------
    id : str
        Stable gene identifier (e.g. 'ENTREZ:7157', 'HGNC:TP53', etc.).
    symbol : str
        Standard gene symbol.
    """
    id: str
    symbol: str


@dataclass(frozen=True)
class Disease:
    """
    Representation of a disease or phenotype.

    Attributes
    ----------
    id : str
        Stable disease identifier (e.g. 'MONDO:0005148', 'DOID:9352').
    name : str
        Human-readable disease name.
    """
    id: str
    name: str


@dataclass(frozen=True)
class DrugTargetAssoc:
    """
    Association between a drug and a target gene.

    Attributes
    ----------
    drug_id : str
        Identifier of the drug (should match Drug.id).
    gene_id : str
        Identifier of the target gene (should match Gene.id).
    source : str
        Source database or file name (e.g. 'toy', 'DrugBank', 'ChEMBL').
    score : Optional[float]
        Optional score (e.g., confidence, binding affinity, etc.).
    """
    drug_id: str
    gene_id: str
    source: str
    score: Optional[float] = None


@dataclass(frozen=True)
class GeneDiseaseAssoc:
    """
    Association between a gene and a disease.

    Attributes
    ----------
    gene_id : str
        Identifier of the gene (should match Gene.id).
    disease_id : str
        Identifier of the disease (should match Disease.id).
    source : str
        Source database or file name (e.g. 'toy', 'DisGeNET').
    score : Optional[float]
        Optional association score (e.g., evidence score).
    """
    gene_id: str
    disease_id: str
    source: str
    score: Optional[float] = None


@dataclass(frozen=True)
class GeneGeneInteraction:
    """
    Interaction between two genes in a network (e.g., PPI).

    Attributes
    ----------
    gene1_id : str
        Identifier of the first gene.
    gene2_id : str
        Identifier of the second gene.
    weight : Optional[float]
        Optional weight representing confidence/strength of interaction.
    source : str
        Source database or network name (default 'ppi').
    """
    gene1_id: str
    gene2_id: str
    weight: Optional[float] = None
    source: str = "ppi"

