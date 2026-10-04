from __future__ import annotations

import networkx as nx
import numpy as np
import pandas as pd


def build_graph(network: pd.DataFrame) -> nx.Graph:
    graph = nx.Graph()
    if network.empty:
        return graph

    for _, row in network.iterrows():
        a = str(row.get("preferredName_A", row.get("stringId_A", ""))).strip()
        b = str(row.get("preferredName_B", row.get("stringId_B", ""))).strip()
        if not a or not b:
            continue
        score = float(row.get("score", 0.0) or 0.0)
        score = max(score, 1e-9)
        graph.add_edge(a, b, weight=score, distance=1.0 / score)
    return graph


def pagerank_numpy(
    graph: nx.Graph,
    alpha: float = 0.85,
    max_iter: int = 200,
    tol: float = 1e-10,
) -> dict:
    """Weighted PageRank implemented with NumPy only.

    NetworkX's PageRank implementation depends on SciPy. Streamlit Community
    Cloud may not have SciPy installed, so this lightweight implementation keeps
    the app deployable without adding that heavy dependency.
    """
    nodes = list(graph.nodes)
    n = len(nodes)
    if n == 0:
        return {}

    index = {node: i for i, node in enumerate(nodes)}
    matrix = np.zeros((n, n), dtype=float)

    for u, v, data in graph.edges(data=True):
        weight = float(data.get("weight", 1.0) or 1.0)
        i, j = index[u], index[v]
        matrix[i, j] += weight
        matrix[j, i] += weight

    row_sums = matrix.sum(axis=1)
    dangling = row_sums == 0

    transition = np.zeros_like(matrix)
    active = ~dangling
    transition[active] = matrix[active] / row_sums[active, None]

    rank = np.full(n, 1.0 / n)
    teleport = np.full(n, (1.0 - alpha) / n)

    for _ in range(max_iter):
        dangling_mass = rank[dangling].sum() / n
        new_rank = alpha * (rank @ transition + dangling_mass) + teleport

        if np.abs(new_rank - rank).sum() < tol:
            rank = new_rank
            break
        rank = new_rank

    total = rank.sum()
    if total > 0:
        rank = rank / total

    return {node: float(rank[index[node]]) for node in nodes}


def centrality_table(graph: nx.Graph) -> pd.DataFrame:
    if graph.number_of_nodes() == 0:
        return pd.DataFrame()

    degree = dict(graph.degree())
    degree_cent = nx.degree_centrality(graph)
    strength = dict(graph.degree(weight="weight"))
    betweenness = nx.betweenness_centrality(graph, weight="distance", normalized=True)
    closeness = nx.closeness_centrality(graph, distance="distance")
    pagerank = pagerank_numpy(graph)

    try:
        eigenvector = nx.eigenvector_centrality(graph, max_iter=2000, weight="weight")
    except (nx.PowerIterationFailedConvergence, nx.NetworkXException):
        eigenvector = {node: np.nan for node in graph.nodes}

    df = pd.DataFrame(
        {
            "Gene": list(graph.nodes),
            "Degree": [degree[n] for n in graph.nodes],
            "Degree centrality": [degree_cent[n] for n in graph.nodes],
            "Weighted degree": [strength[n] for n in graph.nodes],
            "Betweenness": [betweenness[n] for n in graph.nodes],
            "Closeness": [closeness[n] for n in graph.nodes],
            "Eigenvector": [eigenvector[n] for n in graph.nodes],
            "PageRank": [pagerank[n] for n in graph.nodes],
        }
    )

    metrics = ["Degree", "Betweenness", "Closeness", "Eigenvector", "PageRank"]
    ranks = []
    for metric in metrics:
        ranks.append(df[metric].rank(method="average", ascending=False, pct=True))
    df["Composite score"] = 1.0 - pd.concat(ranks, axis=1).mean(axis=1)
    return df.sort_values(["Composite score", "Degree"], ascending=[False, False]).reset_index(drop=True)


def choose_hubs(centrality: pd.DataFrame, metric: str, top_n: int) -> pd.DataFrame:
    if centrality.empty:
        return centrality
    metric = metric if metric in centrality.columns else "Degree"
    return (
        centrality.sort_values([metric, "Degree"], ascending=[False, False])
        .head(min(top_n, len(centrality)))
        .reset_index(drop=True)
    )
