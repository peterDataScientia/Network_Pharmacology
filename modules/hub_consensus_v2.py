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



def build_display_graph(
    network: pd.DataFrame,
    mapping: pd.DataFrame,
    *,
    hide_disconnected_nodes: bool = False,
    use_query_labels: bool = False,
) -> nx.Graph:
    """Build the graph used only for plotting.

    The analysis graph intentionally contains only retained edges. For display,
    users may choose to keep mapped query proteins that have no retained edge.
    This separation prevents a visualization preference from changing centrality
    or the connected/isolated-target classification.
    """
    graph = build_graph(network)
    if hide_disconnected_nodes or mapping.empty:
        return graph

    if use_query_labels and "queryItem" in mapping.columns:
        names = mapping["queryItem"]
    elif "preferredName" in mapping.columns:
        names = mapping["preferredName"]
    elif "stringId" in mapping.columns:
        names = mapping["stringId"]
    else:
        return graph

    for value in names.dropna().astype(str):
        name = value.strip()
        if name:
            graph.add_node(name)

    return graph
