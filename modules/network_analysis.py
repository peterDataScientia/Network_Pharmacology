from __future__ import annotations

import networkx as nx
import pandas as pd

PRIMARY_METRICS = ["Degree", "Betweenness", "Closeness", "Eigenvector"]


def build_graph(network: pd.DataFrame) -> nx.Graph:
    """Build the connected PPI graph returned by STRING.

    STRING confidence is retained as edge metadata for reporting/visualization,
    but the primary hub analysis below is intentionally unweighted. The selected
    STRING confidence threshold determines which edges enter the graph.
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


def centrality_table(graph: nx.Graph) -> pd.DataFrame:
    """Calculate the four unweighted topology metrics used for consensus hubs."""
    if graph.number_of_nodes() == 0:
        return pd.DataFrame()

    degree = dict(graph.degree())
    betweenness = nx.betweenness_centrality(graph, weight=None, normalized=True)
    closeness = nx.closeness_centrality(graph)

    try:
        eigenvector = nx.eigenvector_centrality(
            graph,
            max_iter=5000,
            tol=1e-10,
            weight=None,
        )
    except (nx.PowerIterationFailedConvergence, nx.NetworkXException):
        eigenvector = {node: float("nan") for node in graph.nodes}

    df = pd.DataFrame(
        {
            "Gene": list(graph.nodes),
            "Degree": [degree[n] for n in graph.nodes],
            "Betweenness": [betweenness[n] for n in graph.nodes],
            "Closeness": [closeness[n] for n in graph.nodes],
            "Eigenvector": [eigenvector[n] for n in graph.nodes],
        }
    )

    for metric in PRIMARY_METRICS:
        df[f"{metric} rank"] = (
            df[metric]
            .rank(method="min", ascending=False, na_option="bottom")
            .astype(int)
        )

    return (
        df.sort_values(
            ["Degree", "Betweenness", "Closeness", "Eigenvector", "Gene"],
            ascending=[False, False, False, False, True],
        )
        .reset_index(drop=True)
    )


def consensus_hub_analysis(
    centrality: pd.DataFrame,
    top_n: int = 10,
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Identify hubs by intersection of Top-N lists from four centralities.

    The procedure mirrors the manuscript workflow:
      Degree ∩ Betweenness ∩ Closeness ∩ Eigenvector.

    Returns:
      all_ranked: every connected gene with metric ranks and consensus membership
      consensus_hubs: genes appearing in Top-N for all four metrics
      effective_top_n: min(requested Top-N, connected node count)
    """
    if centrality.empty:
        return pd.DataFrame(), pd.DataFrame(), 0

    effective_top_n = min(int(top_n), len(centrality))
    out = centrality.copy()

    membership_columns = []
    for metric in PRIMARY_METRICS:
        ordered = (
            out.sort_values(
                [metric, "Degree", "Gene"],
                ascending=[False, False, True],
                na_position="last",
            )
            .head(effective_top_n)["Gene"]
            .tolist()
        )
        selected = set(ordered)
        col = f"Top {effective_top_n} {metric}"
        out[col] = out["Gene"].isin(selected)
        membership_columns.append(col)

    out["Consensus count"] = out[membership_columns].sum(axis=1).astype(int)
    out["Consensus"] = out["Consensus count"].astype(str) + "/4"
    out["4/4 consensus hub"] = out["Consensus count"] == 4

    rank_cols = [f"{m} rank" for m in PRIMARY_METRICS]
    out["Mean rank"] = out[rank_cols].mean(axis=1)

    out = (
        out.sort_values(
            ["Consensus count", "Mean rank", "Degree rank", "Gene"],
            ascending=[False, True, True, True],
        )
        .reset_index(drop=True)
    )

    hubs = out[out["4/4 consensus hub"]].copy().reset_index(drop=True)
    return out, hubs, effective_top_n


def connected_and_isolated_targets(
    mapping: pd.DataFrame,
    graph: nx.Graph,
) -> tuple[list[str], list[str]]:
    """Return mapped target names split into connected and isolated-at-threshold."""
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
