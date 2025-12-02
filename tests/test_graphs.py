# tests/test_graphs.py

from pathlib import Path

from ddh.io_handlers import load_toy_data
from ddh.graphs import build_ppi_graph


def test_build_ppi_graph_nodes_and_edges():
    data_dir = Path("data/toy")
    drugs, genes, diseases, dts, gds, ppis = load_toy_data(data_dir)

    G = build_ppi_graph(ppis)

    # The graph should not be empty
    assert G.number_of_nodes() > 0
    assert G.number_of_edges() > 0

    # Check that known nodes exist
    assert "G1" in G.nodes
    assert "G2" in G.nodes
    assert "G3" in G.nodes
    assert "G4" in G.nodes

    # Check known edges exist from toy_ppi.csv
    assert G.has_edge("G1", "G2")
    assert G.has_edge("G1", "G3")
    assert G.has_edge("G2", "G3")
    assert G.has_edge("G3", "G4")


def test_ppi_edge_weights():
    data_dir = Path("data/toy")
    _, _, _, _, _, ppis = load_toy_data(data_dir)

    G = build_ppi_graph(ppis)

    # Check that weights were set
    w_12 = G["G1"]["G2"]["weight"]
    w_13 = G["G1"]["G3"]["weight"]
    w_23 = G["G2"]["G3"]["weight"]
    w_34 = G["G3"]["G4"]["weight"]

    assert w_12 == 0.8
    assert w_13 == 0.6
    assert w_23 == 0.5
    assert w_34 == 0.4

