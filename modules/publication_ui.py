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
    run_publication_enrichment,
)
from modules.publication_plots import (
    CATEGORY_SPECS,
    publication_enrichment_figure,
)
from modules.publication_settings import (
    DEFAULT_PUBLICATION_ANALYSIS_SETTINGS,
    PublicationSettingsError,
    normalize_publication_analysis_settings,
)
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

DATABASE_OPTIONS = {
    "GO Biological Process": "BP",
    "GO Cellular Component": "CC",
    "GO Molecular Function": "MF",
    "Reactome": "REACTOME",
}

P_ADJUST_OPTIONS = {
    "Benjamini–Hochberg (BH)": "BH",
    "Bonferroni": "bonferroni",
    "Holm": "holm",
    "Benjamini–Yekutieli (BY)": "BY",
}

TABLE_CATEGORY_KEYS = {
    "GO-BP": ("go_bp_all", "go_bp_raw", "go_bp_reduced"),
    "GO-CC": ("go_cc_all", "go_cc_raw", "go_cc_reduced"),
    "GO-MF": ("go_mf_all", "go_mf_raw", "go_mf_reduced"),
    "Reactome": ("reactome_all", "reactome_raw", "reactome_reduced"),
}


def publication_methods_text(summary: dict) -> str:
    versions = summary.get("versions", {})
    settings = summary.get("analysis_settings", {})
    ontologies = settings.get("ontologies", ["BP", "CC", "MF"])
    include_reactome = bool(settings.get("include_reactome", True))
    databases = [f"GO-{x}" for x in ontologies]
    if include_reactome:
        databases.append("Reactome")

    bp_preselect = settings.get("go_bp_preselect_n", 100)
    bp_preselect_text = (
        "all significant GO-BP terms"
        if bp_preselect in (None, "all")
        else f"strongest {bp_preselect} significant GO-BP terms"
    )

    go_reduction = (
        f"Wang semantic similarity cutoff {settings.get('go_similarity_cutoff', 0.70):.2f}; "
        "representative selected by minimum adjusted P"
        if settings.get("reduce_go", True)
        else "disabled"
    )
    reactome_reduction = (
        f"Jaccard gene-set similarity cutoff "
        f"{settings.get('reactome_jaccard_cutoff', 0.70):.2f}; "
        "average-linkage clustering; representative selected by minimum adjusted P then maximum gene count"
        if settings.get("reduce_reactome", True)
        else "disabled"
    )

    return (
        "PUBLICATION ENRICHMENT\n\n"
        f"Organism: {summary.get('organism')}\n"
        f"Annotation database: {summary.get('organism_db')}\n"
        f"Background: {summary.get('background_description')}\n"
        f"Databases/ontologies: {', '.join(databases) if databases else 'none'}\n"
        "Identifier mapping: gene symbol to Entrez ID using AnnotationDbi/OrgDb\n"
        f"Multiple testing: {settings.get('p_adjust_method', 'BH')}\n"
        f"Significance: adjusted P < {settings.get('p_adjust_cutoff', 0.05)}\n"
        f"Gene-set size: {settings.get('min_gs_size', 10)} to "
        f"{settings.get('max_gs_size', 500)} genes\n"
        "GO enrichment: clusterProfiler::enrichGO for enabled ontologies\n"
        "Reactome enrichment: ReactomePA::enrichPathway when enabled\n"
        f"GO redundancy reduction: {go_reduction}\n"
        f"GO-BP semantic-reduction input: {bp_preselect_text}\n"
        f"Reactome redundancy reduction: {reactome_reduction}\n"
        f"R: {versions.get('R', 'not recorded')}\n"
        f"clusterProfiler: {versions.get('clusterProfiler', 'not recorded')}\n"
        f"ReactomePA: {versions.get('ReactomePA', 'not recorded')}\n"
        f"AnnotationDbi: {versions.get('AnnotationDbi', 'not recorded')}\n"
        f"GOSemSim: {versions.get('GOSemSim', 'not recorded')}\n"
        f"Organism annotation package version: {versions.get('organism_db', 'not recorded')}\n"
        f"Execution backend: {summary.get('execution', {}).get('executor', 'local-r')}\n"
        f"Execution revision: {summary.get('execution', {}).get('ref_sha', 'not recorded')}\n"
    )


def _default_figure_settings(top_n_per_category: int = 10) -> dict:
    return {
        "plot_type": "Mirrored",
        "categories": ["GO-MF", "GO-CC", "GO-BP", "Reactome"],
        "top_n_by_category": {
            "GO-MF": int(top_n_per_category),
            "GO-CC": int(top_n_per_category),
            "GO-BP": int(top_n_per_category),
            "Reactome": int(top_n_per_category),
        },
        "sort_by": "Adjusted P-value",
        "display_cutoff": None,
        "selected_term_ids": {},
        "wrap_width": 48,
        "manual_selection": False,
    }


def publication_export_files(
    result: PublicationEnrichmentResult,
    top_n_per_category: int = 10,
    figure_settings: dict | None = None,
) -> dict[str, bytes]:
    figure_settings = {
        **_default_figure_settings(top_n_per_category),
        **(figure_settings or {}),
    }

    files: dict[str, bytes] = {
        "publication_enrichment/summary.json": json.dumps(
            result.summary, indent=2
        ).encode("utf-8"),
        "publication_enrichment/METHODS.txt": publication_methods_text(
            result.summary
        ).encode("utf-8"),
        "publication_enrichment/figure_settings.json": json.dumps(
            figure_settings,
            indent=2,
        ).encode("utf-8"),
    }
    for key, df in result.tables.items():
        files[f"publication_enrichment/tables/{key}.csv"] = df.to_csv(
            index=False
        ).encode("utf-8")

    fig = publication_enrichment_figure(
        result.tables,
        plot_type=figure_settings.get("plot_type", "Mirrored"),
        categories=figure_settings.get("categories"),
        top_n_by_category=figure_settings.get("top_n_by_category"),
        sort_by=figure_settings.get("sort_by", "Adjusted P-value"),
        display_cutoff=figure_settings.get("display_cutoff"),
        selected_term_ids=figure_settings.get("selected_term_ids"),
        wrap_width=int(figure_settings.get("wrap_width", 48)),
    )
    for fmt in ["png", "pdf", "svg"]:
        files[
            f"publication_enrichment/figures/publication_enrichment.{fmt}"
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


def _filter_by_cutoff(df: pd.DataFrame, cutoff: float) -> pd.DataFrame:
    if df.empty or "p.adjust" not in df.columns:
        return df.copy()
    out = df.copy()
    out["p.adjust"] = pd.to_numeric(out["p.adjust"], errors="coerce")
    out = out[
        out["p.adjust"].notna()
        & (out["p.adjust"] < float(cutoff))
    ]
    if "Count" in out.columns:
        out["Count"] = pd.to_numeric(out["Count"], errors="coerce")
        out = out.sort_values(
            ["p.adjust", "Count"],
            ascending=[True, False],
        )
    else:
        out = out.sort_values("p.adjust")
    return out.reset_index(drop=True)


def _view_tables(
    tables: dict[str, pd.DataFrame],
    display_cutoff: float,
) -> dict[str, pd.DataFrame]:
    out = dict(tables)
    for _, (all_key, raw_key, reduced_key) in TABLE_CATEGORY_KEYS.items():
        all_df = tables.get(all_key, pd.DataFrame())
        if not all_df.empty:
            out[raw_key] = _filter_by_cutoff(all_df, display_cutoff)
        out[reduced_key] = _filter_by_cutoff(
            tables.get(reduced_key, pd.DataFrame()),
            display_cutoff,
        )
    return out


def _render_database_counts(
    tables: dict[str, pd.DataFrame],
) -> None:
    rows = []
    for category, (all_key, raw_key, reduced_key) in TABLE_CATEGORY_KEYS.items():
        tested = len(tables.get(all_key, pd.DataFrame()))
        significant = len(tables.get(raw_key, pd.DataFrame()))
        reduced = len(tables.get(reduced_key, pd.DataFrame()))
        if tested or significant or reduced:
            rows.append(
                {
                    "Category": category,
                    "Tested": tested,
                    "Significant": significant,
                    "Non-redundant": reduced,
                }
            )
    if rows:
        st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)


def render_publication_enrichment(
    mapping: pd.DataFrame,
    settings: dict,
) -> None:
    st.markdown("### Publication enrichment")
    st.caption(
        "Publication-oriented R/Bioconductor over-representation analysis with explicit "
        "background, multiple-testing control, redundancy reduction and reproducibility metadata."
    )

    foreground = _submitted_symbols(settings)
    foreground_source = "submitted target symbols"
    if not foreground:
        foreground = _mapped_symbols(mapping)
        foreground_source = "STRING-mapped preferred gene symbols (fallback)"

    if len(foreground) < 2:
        st.warning("At least two target symbols are required.")
        return

    taxon_id = int(settings["taxon_id"])
    executor = publication_execution_mode()
    if executor == "unavailable":
        st.warning(
            "Publication Enrichment compute is not configured and Rscript is not available locally. "
            "Quick STRING enrichment remains fully available."
        )
        return

    ready, missing = publication_environment_status(taxon_id)
    if not ready:
        st.error(
            "Publication Enrichment compute is not ready. "
            "Quick STRING enrichment remains available."
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
            "The background defines which genes had an opportunity to enter the foreground. "
            "Changing it changes the statistical question."
        ),
    )

    custom_background: list[str] = []
    if background_label == "Custom study background":
        left, right = st.columns([1.2, 1])
        with left:
            bg_text = st.text_area(
                "Paste custom background genes",
                placeholder="All genes/proteins that could realistically have entered the foreground...",
                height=150,
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
            "Recommended general reference for prediction/database-derived target lists: "
            "all Entrez IDs represented in the selected organism annotation package."
        )
    else:
        st.warning(
            "Package/default background is retained mainly for reproducibility of older analyses. "
            "For a new study, prefer an explicit scientifically justified background."
        )

    with st.expander("Basic Settings", expanded=True):
        selected_databases = st.multiselect(
            "Databases / ontologies",
            list(DATABASE_OPTIONS),
            default=list(DATABASE_OPTIONS),
            key="publication_databases",
        )

        p_adjust_label = st.selectbox(
            "Multiple-testing correction",
            list(P_ADJUST_OPTIONS),
            index=0,
            key="publication_p_adjust_method",
        )

        cutoff_mode = st.selectbox(
            "Adjusted P-value threshold",
            ["0.001", "0.01", "0.05 — publication default", "0.10", "Custom"],
            index=2,
            key="publication_cutoff_mode",
        )
        if cutoff_mode == "Custom":
            analysis_cutoff = float(
                st.number_input(
                    "Custom adjusted P-value threshold",
                    min_value=0.0001,
                    max_value=1.0,
                    value=0.05,
                    step=0.001,
                    format="%.4f",
                    key="publication_custom_cutoff",
                )
            )
        else:
            analysis_cutoff = float(cutoff_mode.split()[0])

    selected_codes = [DATABASE_OPTIONS[x] for x in selected_databases]
    selected_ontologies = [x for x in selected_codes if x in {"BP", "CC", "MF"}]
    include_reactome = "REACTOME" in selected_codes

    with st.expander("Advanced Settings", expanded=False):
        s1, s2 = st.columns(2)
        min_gs_size = int(
            s1.number_input(
                "Minimum genes per term",
                min_value=1,
                max_value=5000,
                value=10,
                step=1,
                key="publication_min_gs_size",
            )
        )
        max_gs_size = int(
            s2.number_input(
                "Maximum genes per term",
                min_value=1,
                max_value=10000,
                value=500,
                step=10,
                key="publication_max_gs_size",
            )
        )

        reduce_go = False
        go_similarity_cutoff = 0.70
        bp_preselect_n: int | None = 100
        if selected_ontologies:
            reduce_go = st.checkbox(
                "Reduce redundant GO terms",
                value=True,
                key="publication_reduce_go",
            )
            if reduce_go:
                go_similarity_cutoff = float(
                    st.slider(
                        "GO Wang semantic-similarity cutoff",
                        min_value=0.50,
                        max_value=0.95,
                        value=0.70,
                        step=0.05,
                        key="publication_go_cutoff",
                    )
                )
                if "BP" in selected_ontologies:
                    bp_mode = st.selectbox(
                        "GO-BP terms entering semantic reduction",
                        [
                            "All significant terms",
                            "50",
                            "100 — validated reference setting",
                            "200",
                            "Custom",
                        ],
                        index=2,
                        key="publication_bp_preselect_mode",
                    )
                    if bp_mode == "All significant terms":
                        bp_preselect_n = None
                    elif bp_mode == "Custom":
                        bp_preselect_n = int(
                            st.number_input(
                                "Custom GO-BP preselection maximum",
                                min_value=1,
                                max_value=5000,
                                value=100,
                                step=10,
                                key="publication_bp_preselect_custom",
                            )
                        )
                    else:
                        bp_preselect_n = int(bp_mode.split()[0])

        reduce_reactome = False
        reactome_jaccard_cutoff = 0.70
        if include_reactome:
            reduce_reactome = st.checkbox(
                "Reduce redundant Reactome pathways",
                value=True,
                key="publication_reduce_reactome",
            )
            if reduce_reactome:
                reactome_jaccard_cutoff = float(
                    st.slider(
                        "Reactome Jaccard similarity cutoff",
                        min_value=0.50,
                        max_value=0.95,
                        value=0.70,
                        step=0.05,
                        key="publication_reactome_cutoff",
                    )
                )

    try:
        analysis_settings = normalize_publication_analysis_settings(
            {
                "ontologies": selected_ontologies,
                "include_reactome": include_reactome,
                "p_adjust_method": P_ADJUST_OPTIONS[p_adjust_label],
                "p_adjust_cutoff": analysis_cutoff,
                "min_gs_size": min_gs_size,
                "max_gs_size": max_gs_size,
                "reduce_go": reduce_go,
                "go_similarity_cutoff": go_similarity_cutoff,
                "go_bp_preselect_n": bp_preselect_n,
                "reduce_reactome": reduce_reactome,
                "reactome_jaccard_cutoff": reactome_jaccard_cutoff,
            }
        )
    except PublicationSettingsError as exc:
        st.error(str(exc))
        return

    bg_mode = {
        "Custom study background": "custom",
        "All annotated genes for organism": "annotated",
        "Package/default background": "default",
    }[background_label]

    st.markdown("#### Analysis preview")
    p1, p2, p3, p4 = st.columns(4)
    p1.metric("Foreground symbols", len(foreground))
    p2.metric(
        "Databases",
        len(selected_ontologies) + int(include_reactome),
    )
    p3.metric("Adjusted P cutoff", f"{analysis_settings['p_adjust_cutoff']:.3g}")
    p4.metric(
        "Background",
        "Custom" if bg_mode == "custom" else bg_mode.title(),
    )
    st.caption(
        "AnnotationDbi symbol→Entrez mapping and foreground/background containment QC "
        "are reported immediately after the run."
    )

    signature_now = {
        "taxon_id": taxon_id,
        "foreground": foreground,
        "background_label": background_label,
        "background_mode": bg_mode,
        "custom_background": custom_background,
        "source_type": source_type,
        "analysis_settings": analysis_settings,
    }
    stored_signature = st.session_state.get("publication_result_signature", {})
    has_result = st.session_state.get("publication_result") is not None
    settings_changed = has_result and stored_signature != signature_now

    if settings_changed:
        st.info(
            "Publication-enrichment settings have changed. The previous result is preserved "
            "but is hidden from interpretation/export until you update the analysis or restore "
            "the matching settings."
        )

    run_label = (
        "Update publication enrichment"
        if has_result
        else "Run publication enrichment"
    )
    if st.button(
        run_label,
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
                    analysis_settings=analysis_settings,
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
                st.session_state["publication_result_signature"] = signature_now
                st.session_state.pop("publication_figure_settings", None)
                st.success("Publication enrichment completed.")

    result = st.session_state.get("publication_result")
    signature = st.session_state.get("publication_result_signature", {})
    if result is None:
        st.session_state["publication_result_current"] = False
        return

    if signature != signature_now:
        st.session_state["publication_result_current"] = False
        st.warning(
            "The controls above differ from the stored result. Results below are from the "
            "previous run. Use **Update publication enrichment** before interpreting or exporting them."
        )
        return

    st.session_state["publication_result_current"] = True
    summary = result.summary
    tables = result.tables
    run_settings = summary.get("analysis_settings", analysis_settings)
    run_cutoff = float(run_settings.get("p_adjust_cutoff", 0.05))
    result_key = str(
        summary.get("execution", {}).get("request_id")
        or summary.get("execution", {}).get("run_id")
        or abs(hash(json.dumps(signature, sort_keys=True, default=str)))
    )

    st.divider()
    st.markdown("### Results workspace")

    c1, c2, c3, c4 = st.columns(4)
    submitted = int(summary.get("foreground_symbols_submitted", 0) or 0)
    mapped = int(summary.get("foreground_entrez_mapped", 0) or 0)
    c1.metric("Submitted", submitted)
    c2.metric("Mapped Entrez", mapped)
    c3.metric(
        "Mapping rate",
        f"{(100.0 * mapped / submitted):.1f}%" if submitted else "—",
    )
    c4.metric(
        "Background genes",
        summary.get("background_entrez_count") or "package default",
    )

    if is_egcg_risi_reference(foreground) and bg_mode == "default":
        observed = summary.get("counts", {})
        expected = EGCG_RISI_ENRICHMENT_DEFAULT_EXPECTED
        if analysis_settings == DEFAULT_PUBLICATION_ANALYSIS_SETTINGS:
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
                st.error("EGCG/RISI reference parity FAILED. " + detail)
            else:
                st.success(
                    "EGCG/RISI reference parity PASSED: "
                    "GO-BP 1432→49 · GO-CC 23→17 · GO-MF 54→26 · Reactome 276→91."
                )

    display_cutoff = float(
        st.number_input(
            "Display adjusted P-value cutoff",
            min_value=0.0001,
            max_value=run_cutoff,
            value=run_cutoff,
            step=min(0.001, run_cutoff / 10),
            format="%.4f",
            key=f"publication_display_cutoff_{result_key}",
            help=(
                "This can be made stricter without rerunning R because the complete tested "
                "tables are saved. To use a looser threshold than the analysis cutoff, update the analysis."
            ),
        )
    )
    view_tables = _view_tables(tables, display_cutoff)

    tab_overview, tab_qc, tab_raw, tab_reduced, tab_redundancy, tab_figure, tab_methods = st.tabs(
        [
            "Overview",
            "Mapping & QC",
            "Significant terms",
            "Non-redundant",
            "Redundancy explorer",
            "Figure editor",
            "Methods & export",
        ]
    )

    with tab_overview:
        st.info(
            f"Background: {summary.get('background_description', 'unknown')} · "
            f"Analysis threshold: adjusted P < {run_cutoff:g} · "
            f"Current display: adjusted P < {display_cutoff:g}"
        )
        _render_database_counts(view_tables)

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
        else:
            st.success("All submitted foreground symbols mapped to Entrez IDs.")

        if "mapping_background" in tables:
            st.markdown("#### Background mapping")
            st.dataframe(
                tables["mapping_background"],
                use_container_width=True,
                hide_index=True,
            )
        if "background_universe" in tables:
            st.caption(
                f"Explicit background universe contains {len(tables['background_universe'])} Entrez IDs."
            )

    with tab_raw:
        for category, (_, raw_key, _) in TABLE_CATEGORY_KEYS.items():
            df = view_tables.get(raw_key, pd.DataFrame())
            if raw_key not in tables and df.empty:
                continue
            with st.expander(
                f"{category} — {len(df)} significant terms",
                expanded=category == "GO-BP",
            ):
                st.dataframe(df, use_container_width=True, hide_index=True)

    with tab_reduced:
        for category, (_, _, reduced_key) in TABLE_CATEGORY_KEYS.items():
            df = view_tables.get(reduced_key, pd.DataFrame())
            if reduced_key not in tables and df.empty:
                continue
            with st.expander(
                f"{category} — {len(df)} non-redundant terms",
                expanded=True,
            ):
                st.dataframe(df, use_container_width=True, hide_index=True)

    with tab_redundancy:
        st.caption(
            "Reduction changes presentation, not the underlying tested enrichment table. "
            "Raw significant and non-redundant tables are both retained."
        )
        rows = []
        for category, (_, raw_key, reduced_key) in TABLE_CATEGORY_KEYS.items():
            raw_df = view_tables.get(raw_key, pd.DataFrame())
            reduced_df = view_tables.get(reduced_key, pd.DataFrame())
            if raw_key in tables or reduced_key in tables:
                rows.append(
                    {
                        "Category": category,
                        "Significant": len(raw_df),
                        "Representatives": len(reduced_df),
                        "Removed as redundant": max(0, len(raw_df) - len(reduced_df)),
                    }
                )
        if rows:
            st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)

        reactome = view_tables.get("reactome_raw", pd.DataFrame())
        if not reactome.empty and "RedundancyCluster" in reactome.columns:
            st.markdown("#### Reactome redundancy clusters")
            cluster_summary = (
                reactome.groupby("RedundancyCluster", dropna=False)
                .agg(
                    Terms=("ID", "count"),
                    Best_adjusted_P=("p.adjust", "min"),
                )
                .reset_index()
                .sort_values(["Terms", "Best_adjusted_P"], ascending=[False, True])
            )
            st.dataframe(cluster_summary, hide_index=True, use_container_width=True)

        st.info(
            f"GO reduction: {'ON' if run_settings.get('reduce_go', True) else 'OFF'}"
            f" · Wang cutoff {run_settings.get('go_similarity_cutoff', 0.70):.2f}. "
            f"Reactome reduction: {'ON' if run_settings.get('reduce_reactome', True) else 'OFF'}"
            f" · Jaccard cutoff {run_settings.get('reactome_jaccard_cutoff', 0.70):.2f}."
        )

    with tab_figure:
        available_categories = [
            category
            for category, (_, _, reduced_key) in TABLE_CATEGORY_KEYS.items()
            if not view_tables.get(reduced_key, pd.DataFrame()).empty
        ]
        if not available_categories:
            st.info("No non-redundant enrichment terms are available for the current display cutoff.")
        else:
            plot_type = st.radio(
                "Figure type",
                ["Mirrored", "Dot plot", "Horizontal bar"],
                horizontal=True,
                key=f"publication_plot_type_{result_key}",
            )
            figure_categories = st.multiselect(
                "Categories shown",
                available_categories,
                default=available_categories,
                key=f"publication_figure_categories_{result_key}",
            )
            sort_by = st.selectbox(
                "Rank terms by",
                ["Adjusted P-value", "Gene count", "Gene ratio"],
                key=f"publication_figure_sort_{result_key}",
            )
            wrap_width = int(
                st.slider(
                    "Term-label wrap width",
                    25,
                    80,
                    48,
                    key=f"publication_wrap_width_{result_key}",
                )
            )

            top_n_by_category: dict[str, int] = {}
            if figure_categories:
                cols = st.columns(min(4, len(figure_categories)))
                for idx, category in enumerate(figure_categories):
                    top_n_by_category[category] = int(
                        cols[idx % len(cols)].number_input(
                            f"{category} Top N",
                            min_value=1,
                            max_value=50,
                            value=10,
                            step=1,
                            key=f"publication_topn_{result_key}_{category}",
                        )
                    )

            manual_selection = st.checkbox(
                "Select figure terms manually",
                value=False,
                key=f"publication_manual_terms_{result_key}",
                help=(
                    "Manual figure selection changes only presentation. The complete statistical "
                    "tables remain in the export package."
                ),
            )
            selected_term_ids: dict[str, list[str]] = {}
            if manual_selection:
                for category in figure_categories:
                    reduced_key = TABLE_CATEGORY_KEYS[category][2]
                    df = view_tables.get(reduced_key, pd.DataFrame())
                    if df.empty or "ID" not in df.columns:
                        continue
                    option_ids = df["ID"].astype(str).tolist()
                    descriptions = {
                        str(row["ID"]): str(row.get("Description", row["ID"]))
                        for _, row in df.iterrows()
                    }
                    chosen = st.multiselect(
                        f"{category} terms",
                        option_ids,
                        default=[],
                        format_func=lambda term_id, desc=descriptions: desc.get(term_id, term_id),
                        key=f"publication_manual_{result_key}_{category}",
                    )
                    selected_term_ids[category] = chosen

            figure_settings = {
                "plot_type": plot_type,
                "categories": figure_categories,
                "top_n_by_category": top_n_by_category,
                "sort_by": sort_by,
                "display_cutoff": display_cutoff,
                "selected_term_ids": selected_term_ids,
                "wrap_width": wrap_width,
                "manual_selection": manual_selection,
            }
            st.session_state["publication_figure_settings"] = figure_settings

            fig = publication_enrichment_figure(
                view_tables,
                plot_type=plot_type,
                categories=figure_categories,
                top_n_by_category=top_n_by_category,
                sort_by=sort_by,
                display_cutoff=display_cutoff,
                selected_term_ids=selected_term_ids,
                wrap_width=wrap_width,
            )
            st.pyplot(fig, use_container_width=True)
            d1, d2, d3 = st.columns(3)
            d1.download_button(
                "PNG (600 dpi)",
                figure_bytes(fig, "png"),
                "publication_enrichment.png",
                "image/png",
                use_container_width=True,
                key="publication_plot_png",
            )
            d2.download_button(
                "PDF",
                figure_bytes(fig, "pdf"),
                "publication_enrichment.pdf",
                "application/pdf",
                use_container_width=True,
                key="publication_plot_pdf",
            )
            d3.download_button(
                "SVG",
                figure_bytes(fig, "svg"),
                "publication_enrichment.svg",
                "image/svg+xml",
                use_container_width=True,
                key="publication_plot_svg",
            )
            if manual_selection:
                st.warning(
                    "The figure contains manually selected terms. This is recorded in "
                    "figure_settings.json; the complete statistical tables are unchanged."
                )
            plt.close(fig)

    with tab_methods:
        methods = publication_methods_text(summary)
        st.code(methods, language=None)

        figure_settings = st.session_state.get(
            "publication_figure_settings",
            _default_figure_settings(),
        )
        memory = BytesIO()
        with ZipFile(memory, "w", ZIP_DEFLATED) as zf:
            for name, payload in publication_export_files(
                result,
                figure_settings=figure_settings,
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
