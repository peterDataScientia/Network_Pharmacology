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
    PublicationEnrichmentResult,
    publication_environment_status,
    publication_execution_mode,
    rscript_available,
    run_publication_enrichment,
)
from modules.publication_plots import mirrored_enrichment_figure
from modules.reference_validation import (
    EGCG_RISI_ENRICHMENT_DEFAULT_EXPECTED,
    is_egcg_risi_reference,
)


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


def publication_methods_text(summary: dict) -> str:
    versions = summary.get("versions", {})
    return (
        "PUBLICATION ENRICHMENT\n\n"
        f"Organism: {summary.get('organism')}\n"
        f"Annotation database: {summary.get('organism_db')}\n"
        f"Background: {summary.get('background_description')}\n"
        "Identifier mapping: gene symbol to Entrez ID using AnnotationDbi/OrgDb\n"
        "GO: clusterProfiler::enrichGO for Biological Process, Cellular Component and Molecular Function\n"
        "Reactome: ReactomePA::enrichPathway\n"
        "Multiple testing: Benjamini-Hochberg\n"
        "Significance: adjusted P < 0.05\n"
        "GO redundancy reduction: Wang semantic similarity cutoff 0.70; "
        "representative selected by minimum adjusted P\n"
        "GO-BP preprocessing: strongest 100 significant BP terms before semantic reduction\n"
        "Reactome redundancy reduction: Jaccard gene-set similarity cutoff 0.70, "
        "average-linkage clustering cut at distance 0.30; representative selected by "
        "minimum adjusted P then maximum gene count\n"
        f"R: {versions.get('R', 'not recorded')}\n"
        f"clusterProfiler: {versions.get('clusterProfiler', 'not recorded')}\n"
        f"ReactomePA: {versions.get('ReactomePA', 'not recorded')}\n"
        f"AnnotationDbi: {versions.get('AnnotationDbi', 'not recorded')}\n"
        f"GOSemSim: {versions.get('GOSemSim', 'not recorded')}\n"
        f"Organism annotation package version: {versions.get('organism_db', 'not recorded')}\n"
        f"Execution backend: {summary.get('execution', {}).get('executor', 'local-r')}\n"
        f"Execution revision: {summary.get('execution', {}).get('ref_sha', 'not recorded')}\n"
    )


def publication_export_files(
    result: PublicationEnrichmentResult,
    top_n_per_category: int = 10,
) -> dict[str, bytes]:
    files: dict[str, bytes] = {
        "publication_enrichment/summary.json": json.dumps(
            result.summary, indent=2
        ).encode("utf-8"),
        "publication_enrichment/METHODS.txt": publication_methods_text(
            result.summary
        ).encode("utf-8"),
    }
    for key, df in result.tables.items():
        files[f"publication_enrichment/tables/{key}.csv"] = df.to_csv(
            index=False
        ).encode("utf-8")

    fig = mirrored_enrichment_figure(
        result.tables,
        top_n_per_category=top_n_per_category,
    )
    for fmt in ["png", "pdf", "svg"]:
        files[
            f"publication_enrichment/figures/publication_enrichment_mirrored.{fmt}"
        ] = figure_bytes(fig, fmt)
    plt.close(fig)
    return files


def _mapped_symbols(mapping: pd.DataFrame) -> list[str]:
    if mapping.empty:
        return []
    if "preferredName" in mapping.columns:
        values = mapping["preferredName"].dropna().astype(str).tolist()
    elif "queryItem" in mapping.columns:
        values = mapping["queryItem"].dropna().astype(str).tolist()
    else:
        values = []
    return list(dict.fromkeys(v.strip() for v in values if v.strip()))


def _submitted_symbols(settings: dict) -> list[str]:
    values = settings.get("submitted_targets", []) or []
    return list(
        dict.fromkeys(
            str(v).strip()
            for v in values
            if str(v).strip()
        )
    )


def render_publication_enrichment(
    mapping: pd.DataFrame,
    settings: dict,
) -> None:
    st.markdown("### Publication enrichment")
    st.caption(
        "Validated R/Bioconductor workflow: clusterProfiler + ReactomePA, "
        "BH-adjusted P < 0.05, GO Wang reduction at 0.70 and Reactome Jaccard reduction at 0.70."
    )

    # Publication enrichment must be independent of STRING alias/preferred-name
    # changes. The validated manuscript workflow used the original submitted
    # gene-symbol list and mapped SYMBOL -> ENTREZ directly with AnnotationDbi.
    foreground = _submitted_symbols(settings)
    foreground_source = "submitted target symbols"
    if not foreground:
        foreground = _mapped_symbols(mapping)
        foreground_source = "STRING-mapped preferred gene symbols (fallback)"

    if len(foreground) < 2:
        st.warning("At least two target symbols are required.")
        return

    st.info(
        f"This analysis uses {len(foreground)} {foreground_source}. "
        "Publication enrichment maps these symbols directly to Entrez IDs with AnnotationDbi; "
        "STRING preferred-name substitutions are not used when submitted targets are available."
    )

    taxon_id = int(settings["taxon_id"])
    executor = publication_execution_mode()

    if executor == "unavailable":
        st.warning(
            "Publication Enrichment compute is not configured and Rscript is not available locally. "
            "Quick STRING enrichment remains fully available."
        )
        return

    ready, missing = publication_environment_status(taxon_id)
    if ready:
        if executor == "github-actions":
            st.success(
                "GitHub Actions publication runner is ready. "
                "Heavy R/Bioconductor analysis will run on a temporary hosted runner."
            )
        else:
            st.success("Local R/Bioconductor publication environment is ready.")
    else:
        if executor == "github-actions":
            st.error(
                "The GitHub Actions publication runner is not ready. "
                "Quick STRING enrichment remains available."
            )
        else:
            st.error(
                "Publication Enrichment is not ready on this Streamlit server. "
                "Automatic Bioconductor installation from the analysis button is disabled."
            )
        if missing:
            st.caption("Status: " + ", ".join(missing))
        return

    source_type = st.selectbox(
        "Where did this target list come from?",
        list(SOURCE_HELP),
        index=1,
        key="publication_source_type",
    )
    st.info(SOURCE_HELP[source_type])

    if source_type in {"Experimental omics", "Candidate / panel study"}:
        default_bg = "Custom study background"
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
        key="publication_background_mode",
        help=(
            "The background is part of the statistical question. Different backgrounds can change "
            "enrichment ratios, P-values and term rankings."
        ),
    )

    custom_background: list[str] = []
    if background_label == "Custom study background":
        left, right = st.columns([1.2, 1])
        with left:
            bg_text = st.text_area(
                "Paste custom background genes",
                placeholder="All genes/proteins that could realistically have entered the foreground...",
                height=160,
                key="publication_background_text",
            )
        with right:
            bg_upload = st.file_uploader(
                "Or upload custom background",
                type=["txt", "csv", "tsv"],
                key="publication_background_upload",
            )

        seen = set()
        for gene in normalize_targets(bg_text) + targets_from_upload(bg_upload):
            key = gene.upper()
            if key not in seen:
                seen.add(key)
                custom_background.append(gene)
        st.caption(f"Custom background genes detected: {len(custom_background)}")

    elif background_label == "All annotated genes for organism":
        st.success(
            "The universe will be explicitly constructed from all Entrez IDs represented "
            "in the selected organism annotation package."
        )
    else:
        if is_egcg_risi_reference(foreground):
            st.success(
                "Manuscript reproduction mode: this is the validated EGCG/RISI background setting. "
                "Expected reference counts are GO-BP 1432→49, GO-CC 23→17, "
                "GO-MF 54→26 and Reactome 276→91."
            )
        else:
            st.warning(
                "Package/default background is retained mainly for reproducing older analyses. "
                "For a new study, prefer a scientifically justified explicit background when possible."
            )

    figure_top_n = st.slider(
        "Terms per category in publication figure",
        3,
        20,
        10,
        key="publication_figure_top_n",
    )

    bg_mode = {
        "Custom study background": "custom",
        "All annotated genes for organism": "annotated",
        "Package/default background": "default",
    }[background_label]

    if st.button(
        "Run publication enrichment",
        type="primary",
        use_container_width=True,
        key="run_publication_enrichment",
    ):
        if bg_mode == "custom" and not custom_background:
            st.error("Custom background is selected, but no background genes were provided.")
        else:
            status_box = st.status(
                (
                    "Submitting publication enrichment to GitHub Actions…"
                    if executor == "github-actions"
                    else "Running publication enrichment…"
                ),
                expanded=executor == "github-actions",
            )

            def status_update(message: str) -> None:
                status_box.write(message)

            try:
                result = run_publication_enrichment(
                    targets=foreground,
                    taxon_id=taxon_id,
                    background_mode=bg_mode,
                    custom_background=custom_background,
                    auto_install=False,
                    status_callback=status_update,
                )
            except PublicationEnrichmentError as exc:
                status_box.update(
                    label="Publication enrichment failed",
                    state="error",
                    expanded=True,
                )
                st.error(str(exc))
            else:
                status_box.update(
                    label="Publication enrichment completed",
                    state="complete",
                    expanded=False,
                )
                st.session_state["publication_result"] = result
                st.session_state["publication_result_signature"] = {
                    "taxon_id": taxon_id,
                    "foreground": foreground,
                    "background_label": background_label,
                    "background_mode": bg_mode,
                    "custom_background": custom_background,
                    "source_type": source_type,
                }
                st.success("Publication enrichment completed.")

    result = st.session_state.get("publication_result")
    signature = st.session_state.get("publication_result_signature", {})
    if result is None:
        return

    if (
        signature.get("taxon_id") != taxon_id
        or signature.get("foreground") != foreground
        or signature.get("background_mode") != bg_mode
        or signature.get("custom_background", []) != custom_background
        or signature.get("source_type") != source_type
    ):
        st.warning(
            "The displayed Publication Enrichment settings differ from the stored result. "
            "Run Publication Enrichment again before interpreting or exporting these results."
        )
        return

    summary = result.summary
    tables = result.tables

    st.divider()
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Submitted", summary.get("foreground_symbols_submitted", 0))
    c2.metric("Mapped Entrez", summary.get("foreground_entrez_mapped", 0))
    c3.metric("GO-BP reduced", summary.get("counts", {}).get("go_bp_reduced", 0))
    c4.metric("Reactome reduced", summary.get("counts", {}).get("reactome_reduced", 0))

    st.info(
        f"Background: {summary.get('background_description', 'unknown')} · "
        f"Significance: {summary.get('significance', 'BH-adjusted P < 0.05')}"
    )

    if is_egcg_risi_reference(foreground) and bg_mode == "default":
        observed = summary.get("counts", {})
        expected = EGCG_RISI_ENRICHMENT_DEFAULT_EXPECTED
        mismatches = {
            key: (observed.get(key), value)
            for key, value in expected.items()
            if observed.get(key) != value
        }
        if mismatches:
            detail = "; ".join(
                f"{key}: observed {got}, expected {want}"
                for key, (got, want) in mismatches.items()
            )
            st.error(
                "EGCG/RISI manuscript reproducibility check FAILED. "
                "The validated package/default reference is "
                "GO-BP 1432→49, GO-CC 23→17, GO-MF 54→26, Reactome 276→91. "
                + detail
            )
        else:
            st.success(
                "EGCG/RISI manuscript reproducibility check PASSED: "
                "GO-BP 1432→49 · GO-CC 23→17 · GO-MF 54→26 · Reactome 276→91."
            )

    tab_qc, tab_raw, tab_reduced, tab_figure, tab_methods = st.tabs(
        [
            "Mapping & QC",
            "Raw significant",
            "Non-redundant",
            "Publication figure",
            "Methods & export",
        ]
    )

    with tab_qc:
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
            with st.expander(
                f"{label} — {len(df)} significant terms",
                expanded=key == "go_bp_raw",
            ):
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
            with st.expander(
                f"{label} — {len(df)} non-redundant terms",
                expanded=True,
            ):
                st.dataframe(df, use_container_width=True, hide_index=True)

    with tab_figure:
        fig = mirrored_enrichment_figure(
            tables,
            top_n_per_category=figure_top_n,
        )
        st.pyplot(fig, use_container_width=True)
        d1, d2, d3 = st.columns(3)
        d1.download_button(
            "PNG (600 dpi)",
            figure_bytes(fig, "png"),
            "publication_enrichment_mirrored.png",
            "image/png",
            use_container_width=True,
            key="publication_plot_png",
        )
        d2.download_button(
            "PDF",
            figure_bytes(fig, "pdf"),
            "publication_enrichment_mirrored.pdf",
            "application/pdf",
            use_container_width=True,
            key="publication_plot_pdf",
        )
        d3.download_button(
            "SVG",
            figure_bytes(fig, "svg"),
            "publication_enrichment_mirrored.svg",
            "image/svg+xml",
            use_container_width=True,
            key="publication_plot_svg",
        )
        plt.close(fig)

    with tab_methods:
        methods = publication_methods_text(summary)
        st.code(methods, language=None)

        memory = BytesIO()
        with ZipFile(memory, "w", ZIP_DEFLATED) as zf:
            for name, payload in publication_export_files(
                result,
                top_n_per_category=figure_top_n,
            ).items():
                zf.writestr(name, payload)

        st.download_button(
            "Download publication enrichment ZIP",
            memory.getvalue(),
            "publication_enrichment_results.zip",
            "application/zip",
            type="primary",
            use_container_width=True,
            key="publication_enrichment_zip",
        )
