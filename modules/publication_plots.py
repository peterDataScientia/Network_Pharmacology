from __future__ import annotations

import math

import matplotlib.pyplot as plt
import pandas as pd


def mirrored_enrichment_figure(
    reduced_tables: dict[str, pd.DataFrame],
    top_n_per_category: int = 10,
):
    """Publication-style mirrored enrichment figure.

    Left: -log10(BH-adjusted P)
    Right: gene Count
    Categories: GO-MF, GO-CC, GO-BP, Reactome
    """
    specs = [
        ("go_mf_reduced", "GO-MF"),
        ("go_cc_reduced", "GO-CC"),
        ("go_bp_reduced", "GO-BP"),
        ("reactome_reduced", "Reactome"),
    ]

    frames = []
    for key, category in specs:
        df = reduced_tables.get(key, pd.DataFrame()).copy()
        if df.empty or "p.adjust" not in df.columns or "Count" not in df.columns:
            continue
        df["p.adjust"] = pd.to_numeric(df["p.adjust"], errors="coerce")
        df["Count"] = pd.to_numeric(df["Count"], errors="coerce")
        df = df.dropna(subset=["p.adjust", "Count"])
        df = df[df["p.adjust"] > 0]
        df = df.sort_values(["p.adjust", "Count"], ascending=[True, False]).head(
            int(top_n_per_category)
        )
        if df.empty:
            continue
        df["Category"] = category
        df["Label"] = df.get("Description", df.get("ID", "")).astype(str)
        frames.append(df)

    if not frames:
        fig, ax = plt.subplots(figsize=(8, 4), dpi=150)
        ax.text(0.5, 0.5, "No reduced enrichment terms available.", ha="center", va="center")
        ax.axis("off")
        return fig

    data = pd.concat(frames, ignore_index=True)
    data["Significance"] = -data["p.adjust"].map(math.log10)
    data["LeftValue"] = -data["Significance"]
    data["RightValue"] = data["Count"]

    # Keep categories grouped and strongest terms visually prominent.
    cat_order = {"GO-MF": 0, "GO-CC": 1, "GO-BP": 2, "Reactome": 3}
    data["_cat"] = data["Category"].map(cat_order)
    data = data.sort_values(
        ["_cat", "p.adjust", "Count"],
        ascending=[True, True, False],
    ).reset_index(drop=True)
    data = data.iloc[::-1].reset_index(drop=True)

    n = len(data)
    height = max(6.0, 0.32 * n + 1.8)
    fig, ax = plt.subplots(figsize=(12, height), dpi=150)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    category_colors = {
        "GO-MF": "#8C6BB1",
        "GO-CC": "#41AB5D",
        "GO-BP": "#2B8CBE",
        "Reactome": "#E34A33",
    }

    y = list(range(n))
    ax.barh(
        y,
        data["LeftValue"],
        height=0.70,
        color=[category_colors[c] for c in data["Category"]],
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
    left_limit = -(max_sig * 2.15)
    right_limit = max_count * 1.35
    ax.set_xlim(left_limit, right_limit)
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
        r"$-\log_{10}$ (BH-adjusted P-value)",
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
        plt.Rectangle((0, 0), 1, 1, facecolor=category_colors[c], edgecolor="black", label=c)
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
