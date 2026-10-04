from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd


class PublicationEnrichmentError(RuntimeError):
    pass


SOURCE_TYPES = [
    "Experimental omics",
    "Target-prediction/database intersection",
    "Candidate/panel study",
    "Other / I know my background",
]

BACKGROUND_MODES = {
    "all_annotated": "All annotated human genes (org.Hs.eg.db)",
    "custom": "Custom study-specific background",
    "package_default": "Package/default background (legacy reproduction)",
}


def background_recommendation(source_type: str) -> dict[str, str]:
    if source_type == "Experimental omics":
        return {
            "mode": "custom",
            "title": "Use the genes that could actually have been selected",
            "reason": (
                "For transcriptomics/proteomics, the preferred universe is the measured, "
                "detectable or QC-passed genes/proteins—not every annotated human gene."
            ),
        }
    if source_type == "Target-prediction/database intersection":
        return {
            "mode": "all_annotated",
            "title": "Use an explicit annotated-human reference unless a candidate universe is reconstructable",
            "reason": (
                "Database/prediction intersections usually do not have an experimental detection "
                "universe. An explicit all-annotated human background is a transparent general "
                "reference; a reconstructed candidate universe is preferable when defensible."
            ),
        }
    if source_type == "Candidate/panel study":
        return {
            "mode": "custom",
            "title": "Use the assayed or eligible panel",
            "reason": (
                "Only genes/proteins that could have entered the candidate set should define "
                "the enrichment universe."
            ),
        }
    return {
        "mode": "custom",
        "title": "Provide the scientifically eligible universe",
        "reason": (
            "Use a custom background whenever you know which genes/proteins could have entered "
            "the foreground. Use all annotated human genes only when that is the intended reference."
        ),
    }


def _script_path() -> Path:
    return Path(__file__).resolve().parents[1] / "publication_enrichment" / "run_enrichment.R"


def backend_status() -> dict:
    rscript = shutil.which("Rscript")
    script = _script_path()

    status = {
        "available": False,
        "rscript": rscript,
        "script": str(script),
        "missing_packages": [],
        "message": "",
    }

    if not rscript:
        status["message"] = "Rscript is not available in this environment."
        return status
    if not script.exists():
        status["message"] = f"Publication enrichment script not found: {script}"
        return status

    packages = [
        "clusterProfiler",
        "org.Hs.eg.db",
        "ReactomePA",
        "AnnotationDbi",
        "jsonlite",
    ]
    expr = (
        "p <- c(" + ",".join(repr(p) for p in packages) + ");"
        "m <- p[!vapply(p, requireNamespace, logical(1), quietly=TRUE)];"
        "cat(paste(m, collapse='|'))"
    )
    try:
        check = subprocess.run(
            [rscript, "-e", expr],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
    except Exception as exc:
        status["message"] = f"Could not check R packages: {exc}"
        return status

    if check.returncode != 0:
        status["message"] = (check.stderr or check.stdout or "R package check failed").strip()
        return status

    missing = [x for x in check.stdout.strip().split("|") if x]
    status["missing_packages"] = missing
    if missing:
        status["message"] = "Missing R package(s): " + ", ".join(missing)
        return status

    status["available"] = True
    status["message"] = "R/Bioconductor publication enrichment backend is ready."
    return status


def _read_csv_if_present(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def run_publication_enrichment(
    targets: list[str],
    background_mode: str,
    custom_background: list[str] | None = None,
    *,
    fdr_cutoff: float = 0.05,
    go_similarity_cutoff: float = 0.70,
    reactome_similarity_cutoff: float = 0.70,
    analysis_profile: str = "publication",
    timeout_seconds: int = 900,
) -> dict:
    if background_mode not in BACKGROUND_MODES:
        raise PublicationEnrichmentError(f"Unsupported background mode: {background_mode}")
    if analysis_profile not in {"publication", "legacy_notebook"}:
        raise PublicationEnrichmentError(f"Unsupported analysis profile: {analysis_profile}")

    cleaned_targets = list(dict.fromkeys(str(x).strip() for x in targets if str(x).strip()))
    if len(cleaned_targets) < 2:
        raise PublicationEnrichmentError("At least two target gene symbols are required.")

    custom_background = custom_background or []
    cleaned_background = list(
        dict.fromkeys(str(x).strip() for x in custom_background if str(x).strip())
    )
    if background_mode == "custom" and len(cleaned_background) < 2:
        raise PublicationEnrichmentError(
            "Custom background mode requires at least two eligible background gene symbols."
        )

    status = backend_status()
    if not status["available"]:
        raise PublicationEnrichmentError(status["message"])

    rscript = status["rscript"]
    script = _script_path()

    with tempfile.TemporaryDirectory(prefix="publication_enrichment_") as tmp:
        work = Path(tmp)
        foreground_file = work / "foreground.csv"
        background_file = work / "background.csv"
        outdir = work / "results"

        pd.DataFrame({"GeneSymbol": cleaned_targets}).to_csv(foreground_file, index=False)

        background_arg = "NONE"
        if background_mode == "custom":
            pd.DataFrame({"GeneSymbol": cleaned_background}).to_csv(background_file, index=False)
            background_arg = str(background_file)

        cmd = [
            str(rscript),
            str(script),
            str(foreground_file),
            background_mode,
            background_arg,
            str(outdir),
            str(float(fdr_cutoff)),
            str(float(go_similarity_cutoff)),
            str(float(reactome_similarity_cutoff)),
            analysis_profile,
        ]

        run = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
        )
        if run.returncode != 0:
            details = (run.stderr or run.stdout or "Unknown R error").strip()
            raise PublicationEnrichmentError(details)

        summary_path = outdir / "summary.json"
        if not summary_path.exists():
            raise PublicationEnrichmentError(
                "R backend completed without producing summary.json."
            )

        summary = json.loads(summary_path.read_text(encoding="utf-8"))

        table_files = {
            "mapping": "symbol_to_entrez_mapping.csv",
            "unmapped": "unmapped_symbols.csv",
            "background": "background_universe_entrez.csv",
            "background_mapping": "background_symbol_to_entrez_mapping.csv",
            "go_bp_significant": "go_bp_significant_ranked.csv",
            "go_cc_significant": "go_cc_significant_ranked.csv",
            "go_mf_significant": "go_mf_significant_ranked.csv",
            "reactome_significant": "reactome_significant_ranked.csv",
            "go_bp_nonredundant": "go_bp_nonredundant.csv",
            "go_cc_nonredundant": "go_cc_nonredundant.csv",
            "go_mf_nonredundant": "go_mf_nonredundant.csv",
            "reactome_clustered": "reactome_significant_clustered.csv",
            "reactome_nonredundant": "reactome_nonredundant.csv",
            "software_versions": "software_versions.csv",
        }
        tables = {
            key: _read_csv_if_present(outdir / filename)
            for key, filename in table_files.items()
        }

        artifacts = {}
        for path in outdir.iterdir():
            if path.is_file():
                artifacts[path.name] = path.read_bytes()

        return {
            "summary": summary,
            "tables": tables,
            "artifacts": artifacts,
            "stdout": run.stdout,
            "stderr": run.stderr,
        }
