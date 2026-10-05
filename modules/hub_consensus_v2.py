from __future__ import annotations

import networkx as nx
import pandas as pd


def build_graph(network: pd.DataFrame) -> nx.Graph:
    """Build an unweighted graph for visualization only.

    Centrality calculations are performed in R/igraph via GitHub Actions.
    NetworkX is retained here only for local network rendering/layout and for
    identifying which mapped targets have at least one retained PPI edge.
    """
    graph = nx.Graph()
    if network.empty:
        return graph

    for _, row in network.iterrows():
        a = str(row.get("preferredName_A", row.get("stringId_A", ""))).strip()
        b = str(row.get("preferredName_B", row.get("stringId_B", ""))).strip()
        if not a or not b:
            continue
        score = float(row.get("score", 0.0) or 0.0)
        graph.add_edge(a, b, string_confidence=score)

    return graph


def connected_and_isolated_targets(
    mapping: pd.DataFrame,
    graph: nx.Graph,
) -> tuple[list[str], list[str]]:
    if mapping.empty:
        return [], []

    if "preferredName" in mapping.columns:
        mapped = mapping["preferredName"].dropna().astype(str).drop_duplicates().tolist()
    else:
        mapped = mapping["stringId"].dropna().astype(str).drop_duplicates().tolist()

    connected_set = set(graph.nodes)
    connected = [gene for gene in mapped if gene in connected_set]
    isolated = [gene for gene in mapped if gene not in connected_set]
    return connected, isolated
