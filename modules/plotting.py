from __future__ import annotations

from io import BytesIO

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

FIGSIZE = (6.66, 3.90)
EXPORT_DPI = 600


def figure_bytes(fig, fmt: str) -> bytes:
    buf = BytesIO()
    kwargs = {"format": fmt, "bbox_inches": "tight", "facecolor": "white"}
    if fmt.lower() == "png":
        kwargs["dpi"] = EXPORT_DPI
    fig.savefig(buf, **kwargs)
    buf.seek(0)
    return buf.getvalue()


def network_figure(graph: nx.Graph, hubs: set[str] | None = None):
    hubs = hubs or set()
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=150)
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    if graph.number_of_nodes() == 0:
        ax.text(0.5, 0.5, "No interactions passed the selected threshold.", ha="center", va="center")
        ax.axis("off")
        return fig

    pos = nx.spring_layout(graph, seed=42, weight=None)
    degrees = dict(graph.degree())
    sizes = [90 + 30 * degrees[n] for n in graph.nodes]
    node_colors = ["#C43D3D" if n in hubs else "#4D7FB8" for n in graph.nodes]
    widths = [0.55 + 1.6 * graph[u][v].get("string_confidence", 0.0) for u, v in graph.edges]

    nx.draw_networkx_edges(graph, pos, width=widths, alpha=0.35, edge_color="#7C8793", ax=ax)
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
    label_size = 7 if graph.number_of_nodes() <= 40 else 5
    nx.draw_networkx_labels(graph, pos, font_size=label_size, font_weight="bold", ax=ax)
    ax.set_title("STRING protein–protein interaction network", fontsize=11, fontweight="bold")
    ax.axis("off")
    fig.tight_layout()
    return fig


def centrality_figure(hubs: pd.DataFrame, metric: str):
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=150)
    data = hubs.sort_values(metric, ascending=True)
    ax.barh(data["Gene"], data[metric])
    ax.set_xlabel(metric, fontsize=11, fontweight="bold")
    ax.set_ylabel("Gene", fontsize=11, fontweight="bold")
    ax.set_title(f"Top hub genes by {metric}", fontsize=11, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9, width=1.4, length=5, direction="out", pad=4)
    ax.spines["left"].set_linewidth(1.8)
    ax.spines["bottom"].set_linewidth(1.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig



def consensus_centrality_figure(
    ranked: pd.DataFrame,
    top_n: int,
):
    """Four-panel Top-N centrality plot with 4/4 consensus hubs highlighted."""
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
            ranked.sort_values([metric, "Degree", "Gene"], ascending=[False, False, True])
            .head(min(top_n, len(ranked)))
            .sort_values(metric, ascending=True)
        )
        colors = ["#2E9B4D" if gene in consensus else "#D94A45" for gene in data["Gene"]]
        ax.barh(data["Gene"], data[metric], color=colors)
        ax.set_xlabel(metric, fontsize=9, fontweight="bold")
        ax.set_ylabel("")
        ax.set_title(f"({chr(65 + idx)}) {metric}", fontsize=10, fontweight="bold", loc="left")
        ax.tick_params(axis="both", labelsize=8, width=1.2, length=4, direction="out", pad=3)
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

def enrichment_figure(data: pd.DataFrame, title: str, top_n: int = 15):
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=150)
    if data.empty:
        ax.text(0.5, 0.5, "No significant terms.", ha="center", va="center")
        ax.axis("off")
        return fig

    plot = data.copy()
    plot["fdr"] = pd.to_numeric(plot["fdr"], errors="coerce")
    plot["number_of_genes"] = pd.to_numeric(plot["number_of_genes"], errors="coerce").fillna(1)
    plot = plot.dropna(subset=["fdr"]).sort_values("fdr").head(top_n).copy()
    plot["-log10(FDR)"] = -np.log10(plot["fdr"].clip(lower=1e-300))
    plot["label"] = plot["description"].astype(str).str.slice(0, 52)
    plot = plot.sort_values("-log10(FDR)", ascending=True)

    sizes = 20 + 18 * np.sqrt(plot["number_of_genes"].to_numpy())
    sc = ax.scatter(
        plot["-log10(FDR)"],
        plot["label"],
        s=sizes,
        c=plot["-log10(FDR)"],
        cmap="viridis",
        alpha=0.85,
        edgecolors="black",
        linewidths=0.35,
    )
    ax.set_xlabel("−log10(FDR)", fontsize=11, fontweight="bold")
    ax.set_ylabel("")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.tick_params(axis="both", labelsize=8, width=1.4, length=5, direction="out", pad=4)
    ax.spines["left"].set_linewidth(1.8)
    ax.spines["bottom"].set_linewidth(1.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("−log10(FDR)", fontsize=9)
    fig.tight_layout()
    return fig
