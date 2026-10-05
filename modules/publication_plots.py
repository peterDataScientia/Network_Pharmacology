from __future__ import annotations

import math
import textwrap

import matplotlib.pyplot as plt
import pandas as pd


CATEGORY_SPECS = [
    ("go_mf_reduced", "GO-MF"),
    ("go_cc_reduced", "GO-CC"),
    ("go_bp_reduced", "GO-BP"),
    ("reactome_reduced", "Reactome"),
]

CATEGORY_COLORS = {
    "GO-MF": "#8C6BB1",
    "GO-CC": "#41AB5D",
    "GO-BP": "#2B8CBE",
    "Reactome": "#E34A33",
}


def _gene_ratio(value) -> float:
    text = str(value or "").strip()
    if "/" in text:
        left, right = text.split("/", 1)
        try:
            numerator = float(left)
            denominator = float(right)
            return numerator / denominator if denominator else 0.0
        except ValueError:
            return 0.0
    try:
        return float(text)
    except ValueError:
        return 0.0


def prepare_publication_figure_data(
    reduced_tables: dict[str, pd.DataFrame],
    *,
    categories: list[str] | None = None,
    top_n_by_category: dict[str, int] | None = None,
    sort_by: str = "Adjusted P-value",
    display_cutoff: float | None = None,
    selected_term_ids: dict[str, list[str]] | None = None,
    wrap_width: int = 48,
) -> pd.DataFrame:
    categories = categories or [label for _, label in CATEGORY_SPECS]
    top_n_by_category = top_n_by_category or {}
    selected_term_ids = selected_term_ids or {}
    frames: list[pd.DataFrame] = []

    for key, category in CATEGORY_SPECS:
        if category not in categories:
            continue
        df = reduced_tables.get(key, pd.DataFrame()).copy()
        if df.empty or "p.adjust" not in df.columns or "Count" not in df.columns:
            continue

        df["p.adjust"] = pd.to_numeric(df["p.adjust"], errors="coerce")
        df["Count"] = pd.to_numeric(df["Count"], errors="coerce")
        df = df.dropna(subset=["p.adjust", "Count"])
        df = df[df["p.adjust"] > 0]
        if display_cutoff is not None:
            df = df[df["p.adjust"] < float(display_cutoff)]

        ids = [str(x) for x in selected_term_ids.get(category, []) if str(x)]
        if ids and "ID" in df.columns:
            df = df[df["ID"].astype(str).isin(ids)]

        if df.empty:
            continue

        df["GeneRatioValue"] = (
            df["GeneRatio"].map(_gene_ratio)
            if "GeneRatio" in df.columns
            else 0.0
        )
        if sort_by == "Gene count":
            df = df.sort_values(
                ["Count", "p.adjust"],
                ascending=[False, True],
            )
        elif sort_by == "Gene ratio":
            df = df.sort_values(
                ["GeneRatioValue", "p.adjust"],
                ascending=[False, True],
            )
        else:
            df = df.sort_values(
                ["p.adjust", "Count"],
                ascending=[True, False],
            )

        if not ids:
            limit = int(top_n_by_category.get(category, 10))
            df = df.head(max(1, limit))

        df["Category"] = category
        raw_labels = df.get("Description", df.get("ID", "")).astype(str)
        df["Label"] = raw_labels.map(
            lambda value: textwrap.fill(value, width=max(20, int(wrap_width)))
        )
        frames.append(df)

    if not frames:
        return pd.DataFrame()

    data = pd.concat(frames, ignore_index=True)
    data["Significance"] = -data["p.adjust"].map(math.log10)
    order = {"GO-MF": 0, "GO-CC": 1, "GO-BP": 2, "Reactome": 3}
    data["_cat"] = data["Category"].map(order)
    return data


def _empty_figure(message: str = "No reduced enrichment terms available."):
    fig, ax = plt.subplots(figsize=(8, 4), dpi=150)
    ax.text(0.5, 0.5, message, ha="center", va="center")
    ax.axis("off")
    return fig


def mirrored_enrichment_figure(
    reduced_tables: dict[str, pd.DataFrame],
    top_n_per_category: int = 10,
    *,
    categories: list[str] | None = None,
    top_n_by_category: dict[str, int] | None = None,
    sort_by: str = "Adjusted P-value",
    display_cutoff: float | None = None,
    selected_term_ids: dict[str, list[str]] | None = None,
    wrap_width: int = 48,
):
    if top_n_by_category is None:
        top_n_by_category = {
            label: int(top_n_per_category)
            for _, label in CATEGORY_SPECS
        }

    data = prepare_publication_figure_data(
        reduced_tables,
        categories=categories,
        top_n_by_category=top_n_by_category,
        sort_by=sort_by,
        display_cutoff=display_cutoff,
        selected_term_ids=selected_term_ids,
        wrap_width=wrap_width,
    )
    if data.empty:
        return _empty_figure()

    data = data.sort_values(
        ["_cat", "p.adjust", "Count"],
        ascending=[True, True, False],
    )
    data = data.iloc[::-1].reset_index(drop=True)

    n = len(data)
    height = max(6.0, 0.32 * n + 1.8)
    fig, ax = plt.subplots(figsize=(12, height), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    data["LeftValue"] = -data["Significance"]
    data["RightValue"] = data["Count"]
    y = list(range(n))
    ax.barh(
        y,
        data["LeftValue"],
        height=0.70,
        color=[CATEGORY_COLORS[c] for c in data["Category"]],
        edgecolor="black",
        linewidth=0.45,
    )
    ax.barh(
        y,
        data["RightValue"],
        height=0.70,
        color="#404040",
        edgecolor="black",
        linewidth=0.45,
    )
    ax.axvline(0, linewidth=1.2, color="black")

    max_sig = max(float(data["Significance"].max()), 1.0)
    max_count = max(float(data["RightValue"].max()), 1.0)
    ax.set_xlim(-(max_sig * 2.20), max_count * 1.35)
    ax.set_ylim(-1, n - 0.3)

    label_x = -(max_sig * 1.08)
    count_offset = max_count * 0.018
    for i, row in data.iterrows():
        ax.text(
            label_x,
            i,
            row["Label"],
            ha="right",
            va="center",
            fontsize=8.5,
            fontweight="bold",
            clip_on=False,
        )
        ax.text(
            row["RightValue"] + count_offset,
            i,
            f"{int(row['Count'])}",
            ha="left",
            va="center",
            fontsize=8.2,
            fontweight="bold",
        )

    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=9, width=1.0, length=5, direction="out")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.1)
    ax.grid(False)

    ax.text(
        0.25,
        -0.035,
        r"$-\log_{10}$ (adjusted P-value)",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10,
        fontweight="bold",
    )
    ax.text(
        0.75,
        -0.035,
        "Gene Count",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10,
        fontweight="bold",
    )

    handles = [
        plt.Rectangle(
            (0, 0),
            1,
            1,
            facecolor=CATEGORY_COLORS[c],
            edgecolor="black",
            label=c,
        )
        for c in ["GO-MF", "GO-CC", "GO-BP", "Reactome"]
        if c in set(data["Category"])
    ]
    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        fontsize=9,
    )
    fig.subplots_adjust(left=0.50, right=0.84, top=0.98, bottom=0.10)
    return fig


def horizontal_bar_enrichment_figure(
    reduced_tables: dict[str, pd.DataFrame],
    **kwargs,
):
    data = prepare_publication_figure_data(reduced_tables, **kwargs)
    if data.empty:
        return _empty_figure()

    data = data.sort_values(
        ["_cat", "p.adjust", "Count"],
        ascending=[True, True, False],
    ).iloc[::-1].reset_index(drop=True)
    n = len(data)
    fig, ax = plt.subplots(figsize=(10, max(5.0, 0.34 * n + 1.5)), dpi=150)
    y = list(range(n))
    ax.barh(
        y,
        data["Significance"],
        color=[CATEGORY_COLORS[c] for c in data["Category"]],
        edgecolor="black",
        linewidth=0.45,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(data["Label"], fontsize=8.5)
    ax.set_xlabel(r"$-\log_{10}$ (adjusted P-value)", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", direction="out", labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def dot_enrichment_figure(
    reduced_tables: dict[str, pd.DataFrame],
    **kwargs,
):
    data = prepare_publication_figure_data(reduced_tables, **kwargs)
    if data.empty:
        return _empty_figure()

    data = data.sort_values(
        ["_cat", "p.adjust", "Count"],
        ascending=[True, True, False],
    ).iloc[::-1].reset_index(drop=True)
    n = len(data)
    fig, ax = plt.subplots(figsize=(10, max(5.0, 0.34 * n + 1.5)), dpi=150)
    y = list(range(n))
    sizes = 24 + 18 * data["Count"].astype(float)
    sc = ax.scatter(
        data["GeneRatioValue"],
        y,
        s=sizes,
        c=data["Significance"],
        alpha=0.85,
        edgecolors="black",
        linewidths=0.4,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(data["Label"], fontsize=8.5)
    ax.set_xlabel("Gene ratio", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", direction="out", labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label(r"$-\log_{10}$ (adjusted P-value)", fontsize=9)
    fig.tight_layout()
    return fig


def publication_enrichment_figure(
    reduced_tables: dict[str, pd.DataFrame],
    *,
    plot_type: str = "Mirrored",
    categories: list[str] | None = None,
    top_n_by_category: dict[str, int] | None = None,
    sort_by: str = "Adjusted P-value",
    display_cutoff: float | None = None,
    selected_term_ids: dict[str, list[str]] | None = None,
    wrap_width: int = 48,
):
    kwargs = {
        "categories": categories,
        "top_n_by_category": top_n_by_category,
        "sort_by": sort_by,
        "display_cutoff": display_cutoff,
        "selected_term_ids": selected_term_ids,
        "wrap_width": wrap_width,
    }
    if plot_type == "Dot plot":
        return dot_enrichment_figure(reduced_tables, **kwargs)
    if plot_type == "Horizontal bar":
        return horizontal_bar_enrichment_figure(reduced_tables, **kwargs)
    return mirrored_enrichment_figure(reduced_tables, **kwargs)
