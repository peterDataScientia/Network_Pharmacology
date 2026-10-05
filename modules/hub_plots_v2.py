from __future__ import annotations

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

FIGSIZE = (6.66, 3.90)


def network_figure(
    graph: nx.Graph,
    hubs: set[str] | None = None,
    *,
    layout: str = "force_directed",
    show_labels: bool = True,
    label_font_size: int | None = None,
    colorblind_friendly: bool = True,
    center_node_labels: bool = False,
    show_regulatory_signs: bool = True,
):
    hubs = hubs or set()
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    if graph.number_of_nodes() == 0:
        ax.text(
            0.5,
            0.5,
            "No interactions passed the selected settings.",
            ha="center",
            va="center",
        )
        ax.axis("off")
        return fig

    if layout == "circular":
        pos = nx.circular_layout(graph)
    else:
        pos = nx.spring_layout(graph, seed=42, weight=None)

    degree = dict(graph.degree())
    sizes = [90 + 30 * degree[n] for n in graph.nodes]

    if colorblind_friendly:
        hub_color = "#0072B2"
        other_color = "#B8B8B8"
        edge_color = "#5F6368"
    else:
        hub_color = "#2E9B4D"
        other_color = "#4D7FB8"
        edge_color = "#7C8793"

    node_colors = [hub_color if n in hubs else other_color for n in graph.nodes]
    widths = [
        0.55 + 1.6 * float(graph[u][v].get("string_confidence", 0.0) or 0.0)
        for u, v in graph.edges
    ]

    nx.draw_networkx_edges(
        graph,
        pos,
        width=widths,
        alpha=0.42,
        edge_color=edge_color,
        arrows=graph.is_directed(),
        arrowsize=12 if graph.is_directed() else 10,
        connectionstyle="arc3,rad=0.05" if graph.is_directed() else "arc3",
        ax=ax,
    )
    nx.draw_networkx_nodes(
        graph,
        pos,
        node_size=sizes,
        node_color=node_colors,
        edgecolors="white",
        linewidths=0.7,
        alpha=0.95,
        ax=ax,
    )

    if show_labels:
        if label_font_size is None:
            label_size = 7 if graph.number_of_nodes() <= 40 else 5
        else:
            label_size = max(5, min(int(label_font_size), 50))

        if center_node_labels:
            label_pos = pos
            valign = "center"
        else:
            ys = [float(value[1]) for value in pos.values()]
            span = (max(ys) - min(ys)) if ys else 1.0
            offset = max(span * 0.035, 0.025)
            label_pos = {
                node: (float(x), float(y) - offset)
                for node, (x, y) in pos.items()
            }
            valign = "top"

        nx.draw_networkx_labels(
            graph,
            label_pos,
            font_size=label_size,
            font_weight="bold",
            verticalalignment=valign,
            ax=ax,
        )

    if graph.is_directed() and show_regulatory_signs:
        edge_labels = {}
        for u, v, attrs in graph.edges(data=True):
            sign = str(attrs.get("sign", "")).strip().lower()
            if sign == "pos":
                edge_labels[(u, v)] = "+"
            elif sign == "neg":
                edge_labels[(u, v)] = "−"
        if edge_labels:
            nx.draw_networkx_edge_labels(
                graph,
                pos,
                edge_labels=edge_labels,
                font_size=8,
                rotate=False,
                label_pos=0.55,
                bbox={"alpha": 0.0, "edgecolor": "none"},
                ax=ax,
            )

    title = (
        "STRING regulatory network"
        if graph.is_directed()
        else "STRING PPI with 4/4 consensus hubs"
    )
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.axis("off")
    fig.tight_layout()
    return fig


def consensus_centrality_figure(ranked: pd.DataFrame, top_n: int):
    metrics = ["Degree", "Betweenness", "Closeness", "Eigenvector"]
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.2), dpi=150)
    axes = axes.ravel()

    if ranked.empty:
        for ax in axes:
            ax.axis("off")
        axes[0].text(0.5, 0.5, "No connected PPI nodes.", ha="center", va="center")
        fig.tight_layout()
        return fig

    consensus = set(
        ranked.loc[ranked["4/4 consensus hub"], "Gene"].astype(str)
        if "4/4 consensus hub" in ranked.columns
        else []
    )

    for idx, (ax, metric) in enumerate(zip(axes, metrics)):
        data = (
            ranked.sort_values(
                [metric, "Degree", "Gene"],
                ascending=[False, False, True],
            )
            .head(min(top_n, len(ranked)))
            .sort_values(metric, ascending=True)
        )
        colors = ["#2E9B4D" if gene in consensus else "#D94A45" for gene in data["Gene"]]
        ax.barh(data["Gene"], data[metric], color=colors)
        ax.set_xlabel(metric, fontsize=9, fontweight="bold")
        ax.set_ylabel("")
        ax.set_title(
            f"({chr(65 + idx)}) {metric}",
            fontsize=10,
            fontweight="bold",
            loc="left",
        )
        ax.tick_params(
            axis="both",
            labelsize=8,
            width=1.2,
            length=4,
            direction="out",
            pad=3,
        )
        ax.spines["left"].set_linewidth(1.5)
        ax.spines["bottom"].set_linewidth(1.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.text(
        0.5,
        0.01,
        f"Green = 4/4 consensus hub; red = not 4/4 consensus · Top {top_n} per metric",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    return fig
