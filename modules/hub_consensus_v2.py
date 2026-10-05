from __future__ import annotations

import networkx as nx
import pandas as pd


def build_graph(network: pd.DataFrame) -> nx.Graph:
    """Build an unweighted graph for visualization and connectivity checks.

    Regulatory STRING output is represented as a directed graph. Functional and
    physical networks remain undirected. Centrality calculations are performed
    in R/igraph for the validated undirected workflow.
    """
    if network.empty:
        return nx.Graph()

    regulatory = {
        "source_preferred_name",
        "target_preferred_name",
    }.issubset(network.columns)
    graph: nx.Graph = nx.DiGraph() if regulatory else nx.Graph()

    for _, row in network.iterrows():
        a = str(
            row.get(
                "preferredName_A",
                row.get("source_preferred_name", row.get("stringId_A", "")),
            )
        ).strip()
        b = str(
            row.get(
                "preferredName_B",
                row.get("target_preferred_name", row.get("stringId_B", "")),
            )
        ).strip()
        if not a or not b:
            continue

        score = float(
            row.get(
                "analysis_score",
                row.get("score", row.get("combined_score", 0.0)),
            )
            or 0.0
        )
        attrs = {"string_confidence": score}
        sign = str(row.get("sign", "")).strip()
        if sign:
            attrs["sign"] = sign
        graph.add_edge(a, b, **attrs)

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
