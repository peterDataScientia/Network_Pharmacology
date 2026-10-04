from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


class PublicationEnrichmentError(RuntimeError):
    pass


@dataclass
class PublicationEnrichmentResult:
    outdir: Path
    summary: dict
    tables: dict[str, pd.DataFrame]


ORGANISM_CODE = {
    9606: "human",
    10090: "mouse",
    10116: "rat",
}


def rscript_available() -> bool:
    return shutil.which("Rscript") is not None


def _write_lines(path: Path, values: list[str]) -> None:
    path.write_text("\n".join(str(v).strip() for v in values if str(v).strip()) + "\n")


def run_publication_enrichment(
    targets: list[str],
    taxon_id: int,
    background_mode: str,
    custom_background: list[str] | None = None,
    workdir: str | Path | None = None,
) -> PublicationEnrichmentResult:
    """Run the reproducible R/Bioconductor ORA workflow via Rscript.

    background_mode:
      - default: package/default enrichment universe
      - annotated: all Entrez IDs represented in the selected organism OrgDb
      - custom: user-supplied gene symbols
    """
    if taxon_id not in ORGANISM_CODE:
        raise PublicationEnrichmentError(
            "Publication enrichment currently supports human, mouse and rat."
        )

    targets = list(dict.fromkeys(str(x).strip() for x in targets if str(x).strip()))
    if len(targets) < 2:
        raise PublicationEnrichmentError("At least two foreground targets are required.")

    background_mode = background_mode.lower().strip()
    if background_mode not in {"default", "annotated", "custom"}:
        raise PublicationEnrichmentError(
            "Background mode must be default, annotated or custom."
        )

    if background_mode == "custom":
        custom_background = list(
            dict.fromkeys(
                str(x).strip() for x in (custom_background or []) if str(x).strip()
            )
        )
        if not custom_background:
            raise PublicationEnrichmentError(
                "Custom background mode requires a non-empty background gene list."
            )

    if not rscript_available():
        raise PublicationEnrichmentError(
            "Rscript is not available in this environment. "
            "The quick STRING enrichment remains available."
        )

    repo_root = Path(__file__).resolve().parents[1]
    r_script = repo_root / "r" / "enrichment_publication.R"
    if not r_script.exists():
        raise PublicationEnrichmentError(
            f"Publication enrichment R script not found: {r_script}"
        )

    if workdir is None:
        outdir = Path(tempfile.mkdtemp(prefix="publication_enrichment_"))
    else:
        outdir = Path(workdir).resolve()
        outdir.mkdir(parents=True, exist_ok=True)

    foreground_file = outdir / "foreground_symbols.txt"
    _write_lines(foreground_file, targets)

    if background_mode == "custom":
        background_file = outdir / "custom_background_symbols.txt"
        _write_lines(background_file, custom_background or [])
        background_arg = str(background_file)
    else:
        background_arg = "NONE"

    cmd = [
        "Rscript",
        str(r_script),
        str(foreground_file),
        str(outdir),
        ORGANISM_CODE[taxon_id],
        background_mode,
        background_arg,
    ]

    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=900,
        check=False,
    )
    if proc.returncode != 0:
        details = (proc.stderr or proc.stdout or "Unknown R error").strip()
        raise PublicationEnrichmentError(
            "Publication enrichment failed in R. " + details[-4000:]
        )

    summary_path = outdir / "summary.json"
    if not summary_path.exists():
        raise PublicationEnrichmentError(
            "R completed without producing summary.json."
        )

    summary = json.loads(summary_path.read_text())

    table_files = {
        "mapping_foreground": "mapping_foreground.csv",
        "unmapped_foreground": "unmapped_foreground.csv",
        "go_bp_raw": "go_bp_raw_significant.csv",
        "go_bp_reduced": "go_bp_reduced.csv",
        "go_cc_raw": "go_cc_raw_significant.csv",
        "go_cc_reduced": "go_cc_reduced.csv",
        "go_mf_raw": "go_mf_raw_significant.csv",
        "go_mf_reduced": "go_mf_reduced.csv",
        "reactome_raw": "reactome_raw_significant.csv",
        "reactome_reduced": "reactome_reduced.csv",
    }
    if (outdir / "background_universe.csv").exists():
        table_files["background_universe"] = "background_universe.csv"
    if (outdir / "mapping_background.csv").exists():
        table_files["mapping_background"] = "mapping_background.csv"

    tables: dict[str, pd.DataFrame] = {}
    for key, filename in table_files.items():
        path = outdir / filename
        if path.exists():
            try:
                tables[key] = pd.read_csv(path)
            except pd.errors.EmptyDataError:
                tables[key] = pd.DataFrame()

    return PublicationEnrichmentResult(
        outdir=outdir,
        summary=summary,
        tables=tables,
    )
