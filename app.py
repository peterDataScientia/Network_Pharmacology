from __future__ import annotations

import json

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from modules.io_utils import build_results_zip, dataframe_tsv, normalize_targets, targets_from_upload
from modules.hub_consensus_v2 import (
    build_graph,
    connected_and_isolated_targets,
)
from modules.centrality_github import (
    CentralityGitHubError,
    run_r_igraph_centrality,
)
from modules.hub_plots_v2 import consensus_centrality_figure, network_figure
from modules.plotting import enrichment_figure, figure_bytes
from modules.publication_ui import (
    publication_export_files,
    publication_methods_text,
    render_publication_enrichment,
)
from modules.string_api import StringAPIError, run_string_workflow
from modules.string_filters import SOURCE_LABELS
from modules.string_settings_ui import render_string_settings

st.set_page_config(page_title="Network Pharmacology Analyzer", page_icon="🧬", layout="wide")

APP_STATE_VERSION = 8
if st.session_state.get("_app_state_version") != APP_STATE_VERSION:
    st.session_state.pop("analysis", None)
    st.session_state.pop("publication_result", None)
    st.session_state.pop("publication_result_signature", None)
    st.session_state["_app_state_version"] = APP_STATE_VERSION

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
        3. Calculates Degree, Betweenness, Closeness and Eigenvector centrality in R/igraph and identifies 4/4 consensus hubs.
        4. Provides Quick STRING enrichment plus a validated Publication Enrichment workflow using clusterProfiler and ReactomePA.
        5. Reduces redundant GO terms with Wang semantic similarity and Reactome pathways with Jaccard similarity.
        6. Exports raw tables, non-redundant tables, high-resolution PNG/PDF/SVG figures and reproducibility metadata.
        """
    )
    st.info(
        "Your submitted target identifiers are sent to the STRING API to perform mapping, "
        "network retrieval and enrichment."
    )

with st.sidebar:
    string_options = render_string_settings(st)
    st.divider()
    st.header("Analysis settings")
    fdr_cutoff = st.select_slider(
        "Quick STRING enrichment FDR cutoff",
        options=[0.001, 0.01, 0.05, 0.10],
        value=0.05,
    )
    top_n = st.slider(
        "Top N per centrality metric",
        5,
        30,
        10,
        help=(
            "For functional/physical networks, the app takes the Top N genes from "
            "Degree, Betweenness, Closeness and Eigenvector centrality. The existing "
            "4/4 consensus workflow is not applied to directed regulatory networks."
        ),
    )
    enrichment_top_n = st.slider("Terms per enrichment plot", 5, 25, 15)

species_label = string_options["species_label"]
species = string_options["species"]
string_version = string_options["string_version"]
network_type = string_options["network_type"]
network_flavor = string_options["network_flavor"]
active_sources = string_options["active_sources"]
score_label = string_options["required_score"]

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

previous_analysis = st.session_state.get("analysis")
current_signature = {
    "submitted_targets": list(targets),
    "taxon_id": species,
    "string_version": string_version,
    "network_type": network_type,
    "network_flavor": network_flavor,
    "active_sources": list(active_sources),
    "required_score": score_label,
    "first_shell": string_options["first_shell"],
    "second_shell": string_options["second_shell"],
    "typed_physical_edges": string_options["typed_physical_edges"],
    "typed_regulatory_edges": string_options["typed_regulatory_edges"],
    "show_regulatory_signs": string_options["show_regulatory_signs"],
    "fdr_cutoff": fdr_cutoff,
    "top_n": top_n,
    "enrichment_top_n": enrichment_top_n,
}
settings_changed = False
if previous_analysis:
    prior = previous_analysis.get("settings", {})
    settings_changed = any(
        prior.get(key) != value
        for key, value in current_signature.items()
    )
    if settings_changed:
        st.info(
            "Settings or targets have changed. Results below still belong to the previous "
            "run until you click **Update analysis**."
        )

run_label = "Update analysis" if previous_analysis else "Run complete analysis"
run = st.button(run_label, type="primary", use_container_width=True)

if run:
    if not active_sources:
        st.error("Select at least one active STRING interaction source.")
        st.stop()
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
                active_sources=active_sources,
                first_shell=string_options["first_shell"],
                second_shell=string_options["second_shell"],
                typed_physical_edges=string_options["typed_physical_edges"],
                typed_regulatory_edges=string_options["typed_regulatory_edges"],
                show_regulatory_signs=string_options["show_regulatory_signs"],
                bubble_3d=string_options["bubble_3d"],
                block_structure_pics=string_options["block_structure_pics"],
                center_node_labels=string_options["center_node_labels"],
                show_query_node_labels=string_options["show_query_node_labels"],
                hide_disconnected_nodes=string_options["hide_disconnected_nodes"],
                hide_node_labels=string_options["hide_node_labels"],
                label_font_size=string_options["label_font_size"],
                string_version=string_version,
            )
        except StringAPIError as exc:
            st.error(str(exc))
            st.stop()

    graph = build_graph(network)
    connected_targets, isolated_targets = connected_and_isolated_targets(mapping, graph)

    if network.empty:
        centrality = pd.DataFrame()
        consensus_ranked = pd.DataFrame()
        hubs = pd.DataFrame()
        effective_top_n = 0
        centrality_provenance = {
            "centrality_engine": "R/igraph",
            "nodes": 0,
            "edges": 0,
        }
    elif network_type == "regulatory":
        centrality = pd.DataFrame()
        consensus_ranked = pd.DataFrame()
        hubs = pd.DataFrame()
        effective_top_n = 0
        centrality_provenance = {
            "centrality_engine": "not-run",
            "reason": (
                "Directed regulatory STRING networks are shown and exported as directed "
                "graphs; the validated undirected 4/4 consensus workflow is not applied."
            ),
            "nodes": graph.number_of_nodes(),
            "edges": graph.number_of_edges(),
        }
    else:
        centrality_status = st.status(
            "Submitting centrality analysis to R/igraph…",
            expanded=True,
        )

        def centrality_status_update(message: str) -> None:
            centrality_status.write(message)

        try:
            (
                centrality,
                consensus_ranked,
                hubs,
                effective_top_n,
                centrality_provenance,
            ) = run_r_igraph_centrality(
                mapping=mapping,
                network=network,
                taxon_id=species,
                required_score=score_label,
                network_type=network_type,
                string_version=string_version,
                top_n=top_n,
                network_flavor=network_flavor,
                active_sources=active_sources,
                add_nodes=string_options["add_nodes"],
                status_callback=centrality_status_update,
            )
        except CentralityGitHubError as exc:
            centrality_status.update(
                label="R/igraph centrality analysis failed",
                state="error",
                expanded=True,
            )
            st.error(str(exc))
            st.stop()
        else:
            centrality_status.update(
                label="R/igraph centrality analysis completed",
                state="complete",
                expanded=False,
            )

    mapped_queries = set(mapping.get("queryItem", pd.Series(dtype=str)).astype(str).str.upper())
    unresolved = [g for g in targets if g.upper() not in mapped_queries]

    st.session_state.pop("publication_result", None)
    st.session_state.pop("publication_result_signature", None)

    st.session_state["analysis"] = {
        "mapping": mapping,
        "network": network,
        "enrichment": enrichment,
        "native_media": native_media,
        "centrality": centrality,
        "consensus_ranked": consensus_ranked,
        "hubs": hubs,
        "connected_targets": connected_targets,
        "isolated_targets": isolated_targets,
        "effective_top_n": effective_top_n,
        "centrality_provenance": centrality_provenance,
        "unresolved": unresolved,
        "settings": {
            "species": species_label,
            "taxon_id": species,
            "network_type": network_type,
            "network_type_label": string_options["network_type_label"],
            "network_flavor": network_flavor,
            "string_version": string_version,
            "required_score": score_label,
            "active_sources": list(active_sources),
            "evidence_transfer": True,
            "first_shell": string_options["first_shell"],
            "second_shell": string_options["second_shell"],
            "add_nodes": string_options["add_nodes"],
            "layout": string_options["layout"],
            "colorblind_friendly": string_options["colorblind_friendly"],
            "bubble_3d": string_options["bubble_3d"],
            "block_structure_pics": string_options["block_structure_pics"],
            "center_node_labels": string_options["center_node_labels"],
            "show_query_node_labels": string_options["show_query_node_labels"],
            "hide_disconnected_nodes": string_options["hide_disconnected_nodes"],
            "hide_node_labels": string_options["hide_node_labels"],
            "label_font_size": string_options["label_font_size"],
            "typed_physical_edges": string_options["typed_physical_edges"],
            "typed_regulatory_edges": string_options["typed_regulatory_edges"],
            "show_regulatory_signs": string_options["show_regulatory_signs"],
            "fdr_cutoff": fdr_cutoff,
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
        {"highres_png": None, "svg": None, "link": None, "errors": [], "warnings": []},
    )
    centrality = analysis["centrality"]
    consensus_ranked = analysis["consensus_ranked"]
    hubs = analysis["hubs"]
    connected_targets = analysis["connected_targets"]
    isolated_targets = analysis["isolated_targets"]
    effective_top_n = analysis["effective_top_n"]
    centrality_provenance = analysis.get("centrality_provenance", {})
    unresolved = analysis["unresolved"]
    settings = analysis["settings"]
    graph = build_graph(network)

    st.divider()
    st.header("Results")
    source_names = [
        SOURCE_LABELS.get(source, source)
        for source in settings.get("active_sources", [])
    ]
    centrality_label = (
        "not applied to directed regulatory network"
        if settings["network_type"] == "regulatory"
        else "R/igraph"
    )
    st.caption(
        f"STRING v{settings.get('string_version', '12.0')} · "
        f"{settings['network_type']} network · "
        f"{settings.get('network_flavor', 'evidence')} edges · "
        f"score ≥ {settings['required_score']/1000:.3f} · "
        f"added interactors: {settings.get('add_nodes', 0)} · "
        f"centrality: {centrality_label}"
    )
    if source_names:
        st.caption("Active evidence sources · " + " · ".join(source_names))
    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Submitted", len(settings["submitted_targets"]))
    m2.metric("Mapped by STRING", len(mapping))
    m3.metric("Connected", len(connected_targets))
    m4.metric("Network edges", graph.number_of_edges())
    m5.metric("4/4 hubs", len(hubs))

    if unresolved:
        st.warning(
            f"{len(unresolved)} identifier(s) were not mapped: "
            + ", ".join(unresolved[:20])
            + (" …" if len(unresolved) > 20 else "")
        )

    tab_map, tab_net, tab_hub, tab_enrich, tab_export = st.tabs(
        ["Target mapping", "Network", "Hub genes", "Enrichment", "Download package"]
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
        st.subheader("STRING interaction network")
        st.caption(
            "The STRING-native rendering is shown first. The app-generated view uses the "
            "selected evidence-source filter and your local layout/display settings."
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

            if native_media.get("warnings"):
                for warning in native_media["warnings"]:
                    st.warning(warning)
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
                    layout=settings.get("layout", "force_directed"),
                    show_labels=not settings.get("hide_node_labels", False),
                    label_font_size=settings.get("label_font_size", 12),
                    colorblind_friendly=settings.get("colorblind_friendly", True),
                )
                st.pyplot(fig_net, use_container_width=True)
                render_downloads("consensus_hub_highlighted_ppi_network", fig_net)
                plt.close(fig_net)
                if settings["network_type"] == "regulatory":
                    st.caption(
                        "Directed arrows represent STRING regulatory relationships. "
                        "The validated undirected 4/4 hub workflow is not overlaid on this view."
                    )
                else:
                    st.caption(
                        "Highlighted nodes are 4/4 consensus hubs from the validated "
                        "R/igraph workflow. This app-generated view follows the selected "
                        "source filter and layout."
                    )

        if network.empty:
            st.warning(
                "No interactions passed the selected STRING network settings. Try a lower "
                "score threshold, enable additional evidence sources, or review the targets."
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
        st.subheader("Consensus hub-target analysis")
        if settings["network_type"] == "regulatory":
            st.info(
                "This is a directed regulatory network. The existing 4/4 consensus method "
                "was validated for undirected functional/physical topology, so the app does "
                "not silently reuse it here. The directed network remains available in the "
                "Network tab and downloads."
            )
        else:
            st.caption(
                "Primary hub definition: a gene must rank within the Top "
                f"{effective_top_n} connected targets for all four unweighted topology metrics "
                "(Degree, Betweenness, Closeness and Eigenvector), calculated in R/igraph."
            )

        if settings["network_type"] == "regulatory":
            pass
        elif centrality.empty:
            st.warning(
                "Centrality cannot be calculated because no interactions passed the selected settings."
            )
        else:
            if len(centrality) <= effective_top_n:
                st.warning(
                    "The number of connected targets is less than or equal to the selected Top-N. "
                    "Every connected target can therefore enter every Top-N list, so 4/4 consensus "
                    "has little discriminatory value. Reduce Top N or use a larger network."
                )

            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Connected targets", len(centrality))
            c2.metric("Top N per metric", effective_top_n)
            c3.metric("4/4 consensus hubs", len(hubs))
            c4.metric("Isolated at threshold", len(isolated_targets))

            st.markdown("#### 4/4 consensus hubs")
            if hubs.empty:
                st.info("No gene appeared in the Top-N lists of all four centrality measures.")
            else:
                hub_cols = [
                    "Gene",
                    "Degree", "Degree rank",
                    "Betweenness", "Betweenness rank",
                    "Closeness", "Closeness rank",
                    "Eigenvector", "Eigenvector rank",
                    "Consensus", "Mean rank",
                ]
                st.dataframe(hubs[hub_cols], use_container_width=True, hide_index=True)

            fig_hub = consensus_centrality_figure(consensus_ranked, effective_top_n)
            st.pyplot(fig_hub, use_container_width=True)
            render_downloads("consensus_hub_centrality_4panel", fig_hub)
            plt.close(fig_hub)

            st.markdown("#### All connected targets and consensus membership")
            st.dataframe(consensus_ranked, use_container_width=True, hide_index=True)

            if isolated_targets:
                with st.expander(
                    f"Mapped targets with no retained PPI edge at this threshold ({len(isolated_targets)})"
                ):
                    st.write(", ".join(isolated_targets))
                    st.caption(
                        "These targets are retained in the analysis record but are excluded from "
                        "topology-based hub ranking because they have no edge in the filtered PPI network."
                    )

            d1, d2, d3 = st.columns(3)
            d1.download_button(
                "Download consensus hubs",
                dataframe_tsv(hubs),
                "consensus_hubs_4of4.tsv",
                "text/tab-separated-values",
                use_container_width=True,
            )
            d2.download_button(
                "Download all rankings",
                dataframe_tsv(consensus_ranked),
                "centrality_consensus_rankings.tsv",
                "text/tab-separated-values",
                use_container_width=True,
            )
            isolated_df = pd.DataFrame({"Gene": isolated_targets})
            d3.download_button(
                "Download isolated targets",
                dataframe_tsv(isolated_df),
                "isolated_targets_at_threshold.tsv",
                "text/tab-separated-values",
                use_container_width=True,
            )

            engine = centrality_provenance.get("summary", {}).get("engine", {})
            st.info(
                "STRING confidence is used here to decide which PPI edges are retained. "
                "The four primary centrality measures are calculated in R/igraph on the "
                "unweighted filtered topology. "
                f"Engine: {engine.get('R', 'R 4.6.1')} · "
                f"igraph {engine.get('igraph', '2.3.4')}. "
                "A 4/4 hub is a network-topology consensus candidate, not proof of causality "
                "or therapeutic importance."
            )

    with tab_enrich:
        quick_enrich_tab, publication_enrich_tab = st.tabs(
            ["Quick enrichment (STRING)", "Publication enrichment (R/Bioconductor)"]
        )

        with quick_enrich_tab:
            st.subheader("Quick STRING enrichment")
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



        with publication_enrich_tab:
            render_publication_enrichment(mapping, settings)

    with tab_export:
        st.subheader("Complete reproducible results package")
        files = {
            "tables/string_mapping.tsv": dataframe_tsv(mapping),
            "tables/string_network_edges.tsv": dataframe_tsv(network),
            "tables/centrality_all_genes.tsv": dataframe_tsv(centrality),
            "tables/centrality_consensus_rankings.tsv": dataframe_tsv(consensus_ranked),
            "tables/consensus_hubs_4of4.tsv": dataframe_tsv(hubs),
            "tables/isolated_targets_at_threshold.tsv": dataframe_tsv(
                pd.DataFrame({"Gene": isolated_targets})
            ),
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
                layout=settings.get("layout", "force_directed"),
                show_labels=not settings.get("hide_node_labels", False),
                label_font_size=settings.get("label_font_size", 12),
                colorblind_friendly=settings.get("colorblind_friendly", True),
            )
            for fmt in ["png", "pdf", "svg"]:
                files[f"figures/string_interaction_network.{fmt}"] = figure_bytes(fig, fmt)
            plt.close(fig)

        if not centrality.empty:
            fig = consensus_centrality_figure(consensus_ranked, effective_top_n)
            for fmt in ["png", "pdf", "svg"]:
                files[f"figures/consensus_hub_centrality_4panel.{fmt}"] = figure_bytes(fig, fmt)
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

        publication_result = st.session_state.get("publication_result")
        publication_signature = st.session_state.get("publication_result_signature", {})
        publication_matches_current = (
            publication_result is not None
            and publication_signature.get("taxon_id") == settings["taxon_id"]
            and publication_signature.get("foreground")
            == list(
                dict.fromkeys(
                    str(x).strip()
                    for x in settings.get("submitted_targets", [])
                    if str(x).strip()
                )
            )
            and publication_signature.get("background_label")
            == st.session_state.get("publication_background_mode")
            and publication_signature.get("source_type")
            == st.session_state.get("publication_source_type")
        )
        if publication_matches_current:
            publication_top_n = int(st.session_state.get("publication_figure_top_n", 10))
            files.update(
                publication_export_files(
                    publication_result,
                    top_n_per_category=publication_top_n,
                )
            )

        if settings["network_type"] == "regulatory":
            hub_method_line = (
                "Hub analysis: not applied; directed regulatory network retained as directed topology\n"
            )
            engine_line = "Centrality engine: not run for regulatory network\n"
            hub_prose = (
                "Because this analysis used STRING's directed regulatory network, the existing "
                "undirected 4/4 consensus-centrality workflow was not applied. Regulatory edges "
                "and their directions were preserved in the network table and app-generated graph. "
            )
        else:
            hub_method_line = (
                "Hub analysis: unweighted Degree, Betweenness, Closeness and Eigenvector "
                "centrality in R/igraph\n"
            )
            engine = centrality_provenance.get("summary", {}).get("engine", {})
            engine_line = (
                f"Centrality engine: {engine.get('R', 'R 4.6.1')} / "
                f"igraph {engine.get('igraph', '2.3.4')}\n"
            )
            hub_prose = (
                "Degree, betweenness, closeness and eigenvector centrality were calculated "
                "in R using igraph on the resulting unweighted topology. The GitHub Actions "
                "centrality job independently re-fetched the same version-pinned STRING network, "
                "applied the same evidence-source and neighborhood settings, and required an exact "
                "edge-set hash match before accepting the R results. For each metric, the Top-N "
                "connected targets were selected, and targets present in all four Top-N lists were "
                "designated 4/4 consensus hubs. "
            )

        methods = (
            "NETWORK PHARMACOLOGY / TOXICOLOGY ANALYSIS\n\n"
            f"Organism: {settings['species']} (NCBI taxon {settings['taxon_id']})\n"
            f"STRING version: {settings.get('string_version', '12.0')}\n"
            f"STRING network type: {settings['network_type']}\n"
            f"STRING edge meaning: {settings.get('network_flavor', 'evidence')}\n"
            f"Active evidence sources: {', '.join(source_names) if source_names else 'none'}\n"
            f"Evidence transfer: included (public API does not separate direct/transferred channel scores)\n"
            f"Minimum STRING interaction score: {settings['required_score']}/1000\n"
            f"Added interactors: first shell {settings.get('first_shell', 0)}, "
            f"second shell {settings.get('second_shell', 0)} "
            f"(tabular API total add_nodes={settings.get('add_nodes', 0)})\n"
            f"Local network layout: {settings.get('layout', 'force_directed')}\n"
            f"Enrichment significance threshold: FDR <= {settings['fdr_cutoff']}\n"
            + hub_method_line
            + engine_line
            + f"Top N per centrality metric: {effective_top_n}\n"
            + f"4/4 consensus hubs identified: {len(hubs)}\n\n"
            + "STRING identifiers were mapped using get_string_ids. Network interactions were "
            "retrieved from the version-pinned STRING network API using the selected network type, "
            "score threshold and neighborhood size. When the user disabled evidence channels, "
            "the retained edge score was recomputed from the selected STRING channel scores using "
            "STRING's documented prior-corrected probabilistic combination rule, and edges below "
            "the requested confidence threshold were removed. STRING's public image/link API does "
            "not expose evidence-channel filtering, so when a subset of channels was selected the "
            "app-generated network is the authoritative filtered topology and the native STRING "
            "image is presented with that limitation. STRING confidence represents evidence support, "
            "not biochemical interaction strength or binding affinity. "
            + hub_prose
            + "Mapped query targets with no retained interaction were reported separately. "
            "Functional enrichment was obtained from the STRING enrichment API and filtered by "
            "false discovery rate (FDR).\n"
        )
        if publication_matches_current:
            methods += "\n\n" + publication_methods_text(publication_result.summary)

        files["settings.json"] = json.dumps(settings, indent=2).encode("utf-8")
        files["centrality_provenance.json"] = json.dumps(
            centrality_provenance,
            indent=2,
        ).encode("utf-8")
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
