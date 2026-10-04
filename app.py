from __future__ import annotations

import json

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from modules.io_utils import build_results_zip, dataframe_tsv, normalize_targets, targets_from_upload
from modules.network_analysis import build_graph, centrality_table, choose_hubs
from modules.plotting import centrality_figure, enrichment_figure, figure_bytes, network_figure
from modules.string_api import StringAPIError, run_string_workflow

st.set_page_config(page_title="Network Pharmacology Analyzer", page_icon="🧬", layout="wide")

SPECIES = {
    "Homo sapiens (Human)": 9606,
    "Mus musculus (Mouse)": 10090,
    "Rattus norvegicus (Rat)": 10116,
}
CATEGORY_LABELS = {
    "Process": "GO Biological Process",
    "Function": "GO Molecular Function",
    "Component": "GO Cellular Component",
    "KEGG": "KEGG Pathways",
    "RCTM": "Reactome Pathways",
}


def enrichment_subset(df: pd.DataFrame, category: str, fdr_cutoff: float) -> pd.DataFrame:
    if df.empty or "category" not in df.columns:
        return pd.DataFrame()
    out = df[df["category"] == category].copy()
    if "fdr" in out.columns:
        out["fdr"] = pd.to_numeric(out["fdr"], errors="coerce")
        out = out[out["fdr"] <= fdr_cutoff]
        out = out.sort_values("fdr")
    return out.reset_index(drop=True)


def render_downloads(prefix: str, fig):
    c1, c2, c3 = st.columns(3)
    c1.download_button(
        "PNG (600 dpi)",
        figure_bytes(fig, "png"),
        f"{prefix}.png",
        "image/png",
        key=f"{prefix}_png",
    )
    c2.download_button(
        "PDF",
        figure_bytes(fig, "pdf"),
        f"{prefix}.pdf",
        "application/pdf",
        key=f"{prefix}_pdf",
    )
    c3.download_button(
        "SVG",
        figure_bytes(fig, "svg"),
        f"{prefix}.svg",
        "image/svg+xml",
        key=f"{prefix}_svg",
    )


st.title("🧬 Network Pharmacology & Toxicology Analyzer")
st.caption(
    "From a common-target list to STRING PPI, hub genes, GO, KEGG and Reactome results — "
    "with downloadable publication-ready outputs."
)

with st.expander("What this app does", expanded=False):
    st.markdown(
        """
        1. Validates and maps submitted targets with STRING.
        2. Retrieves a STRING protein–protein interaction (PPI) network at your chosen confidence threshold.
        3. Calculates network centrality and ranks candidate hub genes.
        4. Retrieves GO Biological Process, Molecular Function, Cellular Component, KEGG and Reactome enrichment.
        5. Exports raw TSV tables, high-resolution PNG/PDF/SVG figures and a complete ZIP package.
        """
    )
    st.info(
        "Your submitted target identifiers are sent to the STRING API to perform mapping, "
        "network retrieval and enrichment."
    )

with st.sidebar:
    st.header("Analysis settings")
    species_label = st.selectbox("Organism", list(SPECIES), index=0)
    species = SPECIES[species_label]
    network_type = st.selectbox("STRING network type", ["functional", "physical"], index=0)
    flavor_options = ["evidence", "confidence", "actions"]
    if network_type == "functional":
        flavor_options.append("typed")
    network_flavor = st.selectbox(
        "STRING native figure style",
        flavor_options,
        index=0,
        help=(
            "Evidence shows evidence-channel edge colors; confidence emphasizes combined "
            "interaction confidence; actions shows predicted molecular actions; typed is "
            "available for functional networks."
        ),
    )
    score_label = st.select_slider(
        "Minimum STRING interaction score",
        options=[150, 400, 700, 900],
        value=700,
        format_func=lambda x: {
            150: "Low (0.15)",
            400: "Medium (0.40)",
            700: "High (0.70)",
            900: "Highest (0.90)",
        }[x],
    )
    fdr_cutoff = st.select_slider(
        "Enrichment FDR cutoff",
        options=[0.001, 0.01, 0.05, 0.10],
        value=0.05,
    )
    hub_metric = st.selectbox(
        "Hub ranking metric",
        ["Degree", "Betweenness", "Closeness", "Eigenvector", "PageRank", "Composite score"],
        index=0,
    )
    top_n = st.slider("Number of hub genes", 5, 30, 10)
    enrichment_top_n = st.slider("Terms per enrichment plot", 5, 25, 15)

left, right = st.columns([1.3, 1])
with left:
    target_text = st.text_area(
        "Paste common targets",
        placeholder="AKT1\nTP53\nEGFR\nMAPK1\n...",
        height=220,
        help=(
            "Gene symbols, UniProt IDs and other identifiers recognized by STRING are accepted. "
            "Separate by lines, commas, spaces or semicolons."
        ),
    )
with right:
    uploaded = st.file_uploader("Or upload targets", type=["txt", "csv", "tsv"])
    st.markdown("**Example**")
    st.code("AKT1\nTP53\nEGFR\nMAPK1\nSTAT3\nJUN", language=None)

manual_targets = normalize_targets(target_text)
uploaded_targets = targets_from_upload(uploaded)
targets = []
seen = set()
for gene in manual_targets + uploaded_targets:
    key = gene.upper()
    if key not in seen:
        seen.add(key)
        targets.append(gene)

st.write(f"**Detected targets:** {len(targets)}")
if targets:
    st.caption(", ".join(targets[:30]) + (" …" if len(targets) > 30 else ""))

run = st.button("Run complete analysis", type="primary", use_container_width=True)

if run:
    if not targets:
        st.error("Please paste or upload at least one target.")
        st.stop()
    if len(targets) > 500:
        st.error("For this first web version, please submit 500 targets or fewer per analysis.")
        st.stop()

    with st.spinner("Mapping targets and querying STRING…"):
        try:
            mapping, network, enrichment, native_media = run_string_workflow(
                targets,
                species,
                score_label,
                network_type,
                network_flavor=network_flavor,
            )
        except StringAPIError as exc:
            st.error(str(exc))
            st.stop()

    graph = build_graph(network)
    centrality = centrality_table(graph)
    hubs = choose_hubs(centrality, hub_metric, top_n)

    mapped_queries = set(mapping.get("queryItem", pd.Series(dtype=str)).astype(str).str.upper())
    unresolved = [g for g in targets if g.upper() not in mapped_queries]

    st.session_state["analysis"] = {
        "mapping": mapping,
        "network": network,
        "enrichment": enrichment,
        "native_media": native_media,
        "centrality": centrality,
        "hubs": hubs,
        "unresolved": unresolved,
        "settings": {
            "species": species_label,
            "taxon_id": species,
            "network_type": network_type,
            "network_flavor": network_flavor,
            "required_score": score_label,
            "fdr_cutoff": fdr_cutoff,
            "hub_metric": hub_metric,
            "top_n": top_n,
            "enrichment_top_n": enrichment_top_n,
            "submitted_targets": targets,
        },
    }

analysis = st.session_state.get("analysis")
if analysis:
    mapping = analysis["mapping"]
    network = analysis["network"]
    enrichment = analysis["enrichment"]
    native_media = analysis.get(
        "native_media",
        {"highres_png": None, "svg": None, "link": None, "errors": []},
    )
    centrality = analysis["centrality"]
    hubs = analysis["hubs"]
    unresolved = analysis["unresolved"]
    settings = analysis["settings"]
    graph = build_graph(network)

    st.divider()
    st.header("Results")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Submitted", len(settings["submitted_targets"]))
    m2.metric("Mapped by STRING", len(mapping))
    m3.metric("Network nodes", graph.number_of_nodes())
    m4.metric("Network edges", graph.number_of_edges())

    if unresolved:
        st.warning(
            f"{len(unresolved)} identifier(s) were not mapped: "
            + ", ".join(unresolved[:20])
            + (" …" if len(unresolved) > 20 else "")
        )

    tab_map, tab_net, tab_hub, tab_enrich, tab_export = st.tabs(
        ["Target mapping", "PPI network", "Hub genes", "Enrichment", "Download package"]
    )

    with tab_map:
        st.subheader("STRING identifier mapping")
        show_cols = [
            c
            for c in ["queryItem", "preferredName", "stringId", "taxonName", "annotation"]
            if c in mapping.columns
        ]
        st.dataframe(
            mapping[show_cols] if show_cols else mapping,
            use_container_width=True,
            hide_index=True,
        )
        st.download_button(
            "Download mapping TSV",
            dataframe_tsv(mapping),
            "string_mapping.tsv",
            "text/tab-separated-values",
        )

    with tab_net:
        st.subheader("Protein–protein interaction network")
        st.caption(
            "The STRING-native figure is shown first. The second view is generated locally "
            "to highlight the hub genes calculated by this app."
        )

        native_tab, local_tab = st.tabs(
            ["STRING native figure", "Hub-highlighted figure"]
        )

        with native_tab:
            if native_media.get("highres_png"):
                st.image(
                    native_media["highres_png"],
                    caption=(
                        "Official STRING-rendered network · "
                        f"{settings.get('network_flavor', 'evidence')} style"
                    ),
                    use_container_width=True,
                )

                d1, d2, d3 = st.columns(3)
                d1.download_button(
                    "Download STRING high-res PNG",
                    native_media["highres_png"],
                    "string_native_network_highres.png",
                    "image/png",
                    key="string_native_png",
                    use_container_width=True,
                )
                if native_media.get("svg"):
                    d2.download_button(
                        "Download STRING SVG",
                        native_media["svg"],
                        "string_native_network.svg",
                        "image/svg+xml",
                        key="string_native_svg",
                        use_container_width=True,
                    )
                else:
                    d2.caption("SVG unavailable for this run.")

                if native_media.get("link"):
                    d3.link_button(
                        "Open this network in STRING",
                        native_media["link"],
                        use_container_width=True,
                    )
                else:
                    d3.caption("STRING webpage link unavailable.")
            else:
                st.warning(
                    "STRING did not return the native high-resolution image for this run. "
                    "The network data and local figure are still available."
                )

            if native_media.get("errors"):
                with st.expander("STRING media retrieval notes"):
                    for error in native_media["errors"]:
                        st.write(f"• {error}")

        with local_tab:
            if network.empty:
                st.warning(
                    "No interactions passed the selected STRING score threshold, so a "
                    "hub-highlighted local network cannot be constructed."
                )
            else:
                fig_net = network_figure(
                    graph,
                    set(hubs["Gene"]) if not hubs.empty else set(),
                )
                st.pyplot(fig_net, use_container_width=True)
                render_downloads("hub_highlighted_ppi_network", fig_net)
                plt.close(fig_net)
                st.caption(
                    "Red nodes are the current top-ranked hub genes; blue nodes are other "
                    "network proteins. This is an app-generated figure, not STRING's native rendering."
                )

        if network.empty:
            st.warning(
                "No PPI edges passed the selected STRING score threshold. Try a lower "
                "threshold or review the submitted targets."
            )
        else:
            st.markdown("#### STRING interaction table")
            st.dataframe(network, use_container_width=True, hide_index=True)
            st.download_button(
                "Download STRING edges TSV",
                dataframe_tsv(network),
                "string_network_edges.tsv",
                "text/tab-separated-values",
            )

    with tab_hub:
        st.subheader(f"Hub genes ranked by {settings['hub_metric']}")
        if centrality.empty:
            st.warning(
                "Centrality cannot be calculated because no PPI edges passed the selected threshold."
            )
        else:
            st.dataframe(hubs, use_container_width=True, hide_index=True)
            fig_hub = centrality_figure(hubs, settings["hub_metric"])
            st.pyplot(fig_hub, use_container_width=True)
            render_downloads("hub_genes_centrality", fig_hub)
            plt.close(fig_hub)
            st.download_button(
                "Download all centrality metrics",
                dataframe_tsv(centrality),
                "centrality_all_genes.tsv",
                "text/tab-separated-values",
            )
            st.download_button(
                "Download top hub genes",
                dataframe_tsv(hubs),
                "hub_genes.tsv",
                "text/tab-separated-values",
            )
            st.info(
                "Hub status is network- and threshold-dependent. Degree, betweenness, closeness, "
                "eigenvector and PageRank quantify different aspects of centrality; they should "
                "not be interpreted as direct biological causality."
            )

    with tab_enrich:
        st.subheader("Functional enrichment")
        if len(mapping) < 2:
            st.warning(
                "STRING enrichment is not interpreted for a single mapped target. "
                "Submit at least two targets."
            )
        elif enrichment.empty:
            st.warning("No enrichment results were returned.")
        else:
            for category, label in CATEGORY_LABELS.items():
                subset = enrichment_subset(
                    enrichment,
                    category,
                    settings["fdr_cutoff"],
                )
                with st.expander(
                    f"{label} — {len(subset)} significant terms",
                    expanded=category in {"Process", "KEGG", "RCTM"},
                ):
                    if subset.empty:
                        st.caption(f"No terms at FDR ≤ {settings['fdr_cutoff']}.")
                        continue
                    cols = [
                        c
                        for c in [
                            "term",
                            "description",
                            "number_of_genes",
                            "number_of_genes_in_background",
                            "strength",
                            "signal",
                            "fdr",
                            "preferredNames",
                        ]
                        if c in subset.columns
                    ]
                    st.dataframe(
                        subset[cols] if cols else subset,
                        use_container_width=True,
                        hide_index=True,
                    )
                    fig = enrichment_figure(
                        subset,
                        label,
                        settings["enrichment_top_n"],
                    )
                    st.pyplot(fig, use_container_width=True)
                    safe = category.lower()
                    render_downloads(f"enrichment_{safe}", fig)
                    plt.close(fig)
                    st.download_button(
                        f"Download {label} TSV",
                        dataframe_tsv(subset),
                        f"enrichment_{safe}.tsv",
                        "text/tab-separated-values",
                        key=f"tsv_{safe}",
                    )

    with tab_export:
        st.subheader("Complete reproducible results package")
        files = {
            "tables/string_mapping.tsv": dataframe_tsv(mapping),
            "tables/string_network_edges.tsv": dataframe_tsv(network),
            "tables/centrality_all_genes.tsv": dataframe_tsv(centrality),
            "tables/hub_genes.tsv": dataframe_tsv(hubs),
            "tables/enrichment_all.tsv": dataframe_tsv(enrichment),
        }

        if native_media.get("highres_png"):
            files["figures/string_native_network_highres.png"] = native_media["highres_png"]
        if native_media.get("svg"):
            files["figures/string_native_network.svg"] = native_media["svg"]
        if native_media.get("link"):
            files["string_native_network_link.txt"] = native_media["link"].encode("utf-8")

        if graph.number_of_nodes() > 0:
            fig = network_figure(
                graph,
                set(hubs["Gene"]) if not hubs.empty else set(),
            )
            for fmt in ["png", "pdf", "svg"]:
                files[f"figures/string_ppi_network.{fmt}"] = figure_bytes(fig, fmt)
            plt.close(fig)

        if not hubs.empty:
            fig = centrality_figure(hubs, settings["hub_metric"])
            for fmt in ["png", "pdf", "svg"]:
                files[f"figures/hub_genes_centrality.{fmt}"] = figure_bytes(fig, fmt)
            plt.close(fig)

        for category, label in CATEGORY_LABELS.items():
            subset = enrichment_subset(
                enrichment,
                category,
                settings["fdr_cutoff"],
            )
            files[f"tables/enrichment_{category.lower()}.tsv"] = dataframe_tsv(subset)
            if not subset.empty:
                fig = enrichment_figure(
                    subset,
                    label,
                    settings["enrichment_top_n"],
                )
                for fmt in ["png", "pdf", "svg"]:
                    files[f"figures/enrichment_{category.lower()}.{fmt}"] = figure_bytes(fig, fmt)
                plt.close(fig)

        methods = (
            "NETWORK PHARMACOLOGY / TOXICOLOGY ANALYSIS\n\n"
            f"Organism: {settings['species']} (NCBI taxon {settings['taxon_id']})\n"
            f"STRING network type: {settings['network_type']}\n"
            f"STRING native figure style: {settings.get('network_flavor', 'evidence')}\n"
            f"Minimum STRING interaction score: {settings['required_score']}/1000\n"
            f"Enrichment significance threshold: FDR <= {settings['fdr_cutoff']}\n"
            f"Hub ranking metric: {settings['hub_metric']}\n"
            f"Top hub genes requested: {settings['top_n']}\n\n"
            "STRING identifiers were mapped using get_string_ids. The PPI edge list was obtained "
            "from the STRING network API with no added neighbor nodes. The official STRING-native "
            "network was also retrieved as a high-resolution PNG and SVG when available. Network "
            "centrality was calculated locally with NetworkX. Confidence scores were used as edge strengths; "
            "inverse confidence was used as distance for shortest-path-based metrics. Functional "
            "enrichment was obtained from the STRING enrichment API and filtered by false discovery "
            "rate (FDR).\n"
        )
        files["settings.json"] = json.dumps(settings, indent=2).encode("utf-8")
        zip_bytes = build_results_zip(files, methods)
        st.download_button(
            "Download complete results ZIP",
            zip_bytes,
            "network_pharmacology_results.zip",
            "application/zip",
            type="primary",
            use_container_width=True,
        )
        st.code(methods, language=None)

st.divider()
st.caption(
    "This tool supports exploratory/research analysis. Results depend on target identity, organism, "
    "STRING version, interaction threshold and enrichment background; biological conclusions require "
    "domain-specific validation."
)
