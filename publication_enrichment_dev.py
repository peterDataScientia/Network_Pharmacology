from __future__ import annotations

import json
from io import BytesIO
from zipfile import ZIP_DEFLATED, ZipFile

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from modules.io_utils import normalize_targets, targets_from_upload
from modules.plotting import figure_bytes
from modules.publication_enrichment import (
    PublicationEnrichmentError,
    rscript_available,
    run_publication_enrichment,
)
from modules.publication_plots import mirrored_enrichment_figure

st.set_page_config(
    page_title="Publication Enrichment DEV",
    page_icon="🧪",
    layout="wide",
)

SPECIES = {
    "Homo sapiens (Human)": 9606,
    "Mus musculus (Mouse)": 10090,
    "Rattus norvegicus (Rat)": 10116,
}

SOURCE_HELP = {
    "Experimental omics": (
        "Recommended background: genes/proteins that passed assay detection and QC and therefore "
        "could have entered the foreground."
    ),
    "Target-prediction / database intersection": (
        "There is usually no assay-specific detectable universe. An explicit all-annotated organism "
        "background is a transparent general reference, but database/prediction ascertainment bias remains."
    ),
    "Candidate / panel study": (
        "Recommended background: the full candidate/panel gene set from which the selected targets arose."
    ),
    "Other / custom": (
        "Use a custom background when you know which genes had a realistic opportunity to be selected."
    ),
}

st.title("🧪 Publication Enrichment — Development Prototype")
st.caption(
    "Isolated R/Bioconductor prototype. This branch does not replace the production STRING enrichment."
)

if not rscript_available():
    st.warning(
        "Rscript is not currently available in this deployment. The interface is ready, but the "
        "Bioconductor runtime must be installed before Publication Enrichment can execute."
    )

with st.sidebar:
    st.header("Scientific design")
    species_label = st.selectbox("Organism", list(SPECIES), index=0)
    taxon_id = SPECIES[species_label]

    source_type = st.selectbox(
        "Where did this target list come from?",
        list(SOURCE_HELP),
        index=1,
    )
    st.info(SOURCE_HELP[source_type])

    if source_type == "Experimental omics":
        default_bg = "Custom study background"
    elif source_type == "Candidate / panel study":
        default_bg = "Custom study background"
    elif source_type == "Target-prediction / database intersection":
        default_bg = "All annotated genes for organism"
    else:
        default_bg = "All annotated genes for organism"

    bg_options = [
        "Custom study background",
        "All annotated genes for organism",
        "Package/default background",
    ]
    background_label = st.radio(
        "Enrichment background",
        bg_options,
        index=bg_options.index(default_bg),
        help=(
            "The background is part of the statistical question. Different backgrounds can change "
            "enrichment ratios, P-values and term rankings."
        ),
    )

    top_n = st.slider("Terms per category in publication figure", 3, 20, 10)

left, right = st.columns([1.3, 1])
with left:
    target_text = st.text_area(
        "Foreground target genes",
        placeholder="AKT1\nTP53\nSTAT3\n...",
        height=240,
    )
    uploaded = st.file_uploader(
        "Or upload foreground targets",
        type=["txt", "csv", "tsv"],
        key="foreground_upload",
    )

with right:
    st.markdown("### Background")
    custom_bg_upload = None
    custom_bg_text = ""
    if background_label == "Custom study background":
        custom_bg_text = st.text_area(
            "Paste custom background genes",
            placeholder="All genes/proteins that could have been selected...",
            height=170,
        )
        custom_bg_upload = st.file_uploader(
            "Or upload custom background",
            type=["txt", "csv", "tsv"],
            key="background_upload",
        )
    elif background_label == "All annotated genes for organism":
        st.success(
            "The R engine will explicitly construct the universe from all Entrez IDs represented "
            "in the selected organism's OrgDb annotation package."
        )
    else:
        st.warning(
            "Package/default background is provided mainly for reproducing older analyses. "
            "For new studies, an explicit scientifically justified background is preferable."
        )

foreground = []
seen = set()
for gene in normalize_targets(target_text) + targets_from_upload(uploaded):
    key = gene.upper()
    if key not in seen:
        seen.add(key)
        foreground.append(gene)

custom_background = []
if background_label == "Custom study background":
    seen_bg = set()
    for gene in normalize_targets(custom_bg_text) + targets_from_upload(custom_bg_upload):
        key = gene.upper()
        if key not in seen_bg:
            seen_bg.add(key)
            custom_background.append(gene)

st.write(f"**Foreground genes detected:** {len(foreground)}")
if foreground:
    st.caption(", ".join(foreground[:40]) + (" …" if len(foreground) > 40 else ""))

if background_label == "Custom study background":
    st.write(f"**Background genes detected:** {len(custom_background)}")

bg_mode = {
    "Custom study background": "custom",
    "All annotated genes for organism": "annotated",
    "Package/default background": "default",
}[background_label]

run = st.button("Run publication enrichment", type="primary", use_container_width=True)

if run:
    if len(foreground) < 2:
        st.error("Provide at least two foreground genes.")
        st.stop()
    if bg_mode == "custom" and not custom_background:
        st.error("Custom background is selected, but no background genes were provided.")
        st.stop()

    with st.spinner("Running clusterProfiler + ReactomePA in R…"):
        try:
            result = run_publication_enrichment(
                targets=foreground,
                taxon_id=taxon_id,
                background_mode=bg_mode,
                custom_background=custom_background,
            )
        except PublicationEnrichmentError as exc:
            st.error(str(exc))
            st.stop()

    st.session_state["publication_result"] = result

result = st.session_state.get("publication_result")
if result:
    summary = result.summary
    tables = result.tables

    st.divider()
    st.header("Publication enrichment results")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Submitted", summary.get("foreground_symbols_submitted", 0))
    c2.metric("Mapped Entrez", summary.get("foreground_entrez_mapped", 0))
    c3.metric("GO-BP reduced", summary.get("counts", {}).get("go_bp_reduced", 0))
    c4.metric("Reactome reduced", summary.get("counts", {}).get("reactome_reduced", 0))

    st.info(
        f"Background: {summary.get('background_description', 'unknown')} · "
        f"Significance: {summary.get('significance', 'BH-adjusted P < 0.05')}"
    )

    tab_qc, tab_raw, tab_reduced, tab_figure, tab_methods = st.tabs(
        ["Mapping & QC", "Raw significant", "Non-redundant", "Publication figure", "Methods & export"]
    )

    with tab_qc:
        st.subheader("Identifier mapping")
        st.dataframe(
            tables.get("mapping_foreground", pd.DataFrame()),
            use_container_width=True,
            hide_index=True,
        )
        unmapped = tables.get("unmapped_foreground", pd.DataFrame())
        if not unmapped.empty:
            st.warning(f"{len(unmapped)} foreground symbol(s) were not mapped.")
            st.dataframe(unmapped, use_container_width=True, hide_index=True)

        if "background_universe" in tables:
            st.caption(
                f"Explicit background universe contains {len(tables['background_universe'])} Entrez IDs."
            )

    with tab_raw:
        for key, label in [
            ("go_bp_raw", "GO Biological Process"),
            ("go_cc_raw", "GO Cellular Component"),
            ("go_mf_raw", "GO Molecular Function"),
            ("reactome_raw", "Reactome"),
        ]:
            df = tables.get(key, pd.DataFrame())
            with st.expander(f"{label} — {len(df)} significant terms", expanded=key == "go_bp_raw"):
                st.dataframe(df, use_container_width=True, hide_index=True)

    with tab_reduced:
        st.caption(
            "GO: Wang semantic similarity cutoff 0.70. "
            "Reactome: Jaccard gene-set similarity cutoff 0.70."
        )
        for key, label in [
            ("go_bp_reduced", "GO Biological Process"),
            ("go_cc_reduced", "GO Cellular Component"),
            ("go_mf_reduced", "GO Molecular Function"),
            ("reactome_reduced", "Reactome"),
        ]:
            df = tables.get(key, pd.DataFrame())
            with st.expander(f"{label} — {len(df)} non-redundant terms", expanded=True):
                st.dataframe(df, use_container_width=True, hide_index=True)

    with tab_figure:
        fig = mirrored_enrichment_figure(tables, top_n_per_category=top_n)
        st.pyplot(fig, use_container_width=True)
        d1, d2, d3 = st.columns(3)
        d1.download_button(
            "PNG (600 dpi)",
            figure_bytes(fig, "png"),
            "publication_enrichment_mirrored.png",
            "image/png",
            use_container_width=True,
        )
        d2.download_button(
            "PDF",
            figure_bytes(fig, "pdf"),
            "publication_enrichment_mirrored.pdf",
            "application/pdf",
            use_container_width=True,
        )
        d3.download_button(
            "SVG",
            figure_bytes(fig, "svg"),
            "publication_enrichment_mirrored.svg",
            "image/svg+xml",
            use_container_width=True,
        )
        plt.close(fig)

    with tab_methods:
        st.json(summary)

        methods = (
            "PUBLICATION ENRICHMENT\n\n"
            f"Organism: {summary.get('organism')}\n"
            f"Annotation database: {summary.get('organism_db')}\n"
            f"Background: {summary.get('background_description')}\n"
            "Identifier mapping: gene symbol to Entrez ID using AnnotationDbi/OrgDb\n"
            "GO: clusterProfiler::enrichGO for BP, CC and MF\n"
            "Reactome: ReactomePA::enrichPathway\n"
            "Multiple testing: Benjamini-Hochberg\n"
            "Significance: adjusted P < 0.05\n"
            "GO redundancy reduction: Wang semantic similarity cutoff 0.70; "
            "representative selected by minimum adjusted P\n"
            "Reactome redundancy reduction: Jaccard gene-set similarity cutoff 0.70, "
            "average-linkage clustering cut at distance 0.30; representative selected by "
            "minimum adjusted P then maximum gene count\n"
        )
        st.code(methods, language=None)

        memory = BytesIO()
        with ZipFile(memory, "w", ZIP_DEFLATED) as zf:
            zf.writestr("summary.json", json.dumps(summary, indent=2))
            zf.writestr("METHODS.txt", methods)
            for key, df in tables.items():
                zf.writestr(f"tables/{key}.csv", df.to_csv(index=False))
            fig = mirrored_enrichment_figure(tables, top_n_per_category=top_n)
            zf.writestr("figures/publication_enrichment_mirrored.png", figure_bytes(fig, "png"))
            zf.writestr("figures/publication_enrichment_mirrored.pdf", figure_bytes(fig, "pdf"))
            zf.writestr("figures/publication_enrichment_mirrored.svg", figure_bytes(fig, "svg"))
            plt.close(fig)

        st.download_button(
            "Download publication enrichment ZIP",
            memory.getvalue(),
            "publication_enrichment_results.zip",
            "application/zip",
            type="primary",
            use_container_width=True,
        )
