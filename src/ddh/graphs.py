# src/ddh/graphs.py

from __future__ import annotations

from typing import List

import networkx as nx

from .data_models import GeneGeneInteraction


def build_ppi_graph(interactions: List[GeneGeneInteraction]) -> nx.Graph:
    """
    Build an undirected PPI graph from a list of gene–gene interactions.

    Nodes
    -----
    gene_id

    Edges
    -----
    (gene1_id, gene2_id) with attribute 'weight'.

    Parameters
    ----------
    interactions : list[GeneGeneInteraction]
        Interaction records.

    Returns
    -------
    networkx.Graph
        Graph with genes as nodes and interactions as edges.
    """
    G = nx.Graph()
    for inter in interactions:
        w = inter.weight if inter.weight is not None else 1.0
        G.add_edge(inter.gene1_id, inter.gene2_id, weight=w)
    return G

