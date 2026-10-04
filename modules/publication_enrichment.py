from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
import urllib.error
import urllib.request
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

ORGANISM_DB = {
    9606: "org.Hs.eg.db",
    10090: "org.Mm.eg.db",
    10116: "org.Rn.eg.db",
}

BASE_R_PACKAGES = [
    "clusterProfiler",
    "ReactomePA",
    "AnnotationDbi",
    "GOSemSim",
    "jsonlite",
]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _r_env() -> dict[str, str]:
    env = os.environ.copy()
    user_lib = (
        env.get("PUBLICATION_R_LIB")
        or env.get("R_LIBS_USER")
        or str(_repo_root() / ".r-library")
    )
    Path(user_lib).mkdir(parents=True, exist_ok=True)
    env["R_LIBS_USER"] = user_lib
    return env


def rscript_available() -> bool:
    return shutil.which("Rscript") is not None


def publication_backend_url() -> str | None:
    value = os.environ.get("PUBLICATION_BACKEND_URL", "").strip()
    return value.rstrip("/") if value else None


def _remote_health(base_url: str) -> tuple[bool, list[str]]:
    try:
        with urllib.request.urlopen(base_url + "/health", timeout=45) as response:
            payload = json.load(response)
    except Exception as exc:
        return False, [f"backend unavailable: {exc}"]
    if payload.get("status") != "ok":
        return False, ["backend health check failed"]
    return True, []


def publication_environment_status(taxon_id: int) -> tuple[bool, list[str]]:
    if taxon_id not in ORGANISM_CODE:
        return False, ["unsupported organism"]

    backend = publication_backend_url()
    if backend:
        return _remote_health(backend)

    if not rscript_available():
        return False, ["Rscript"]

    packages = BASE_R_PACKAGES + [ORGANISM_DB[taxon_id]]
    r_expr = (
        "pkgs <- c(" + ",".join(json.dumps(p) for p in packages) + "); "
        "missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; "
        "cat(paste(missing, collapse='\\n'))"
    )
    try:
        proc = subprocess.run(
            ["Rscript", "-e", r_expr],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
            env=_r_env(),
        )
    except subprocess.TimeoutExpired:
        return False, ["R environment check timed out"]
    if proc.returncode != 0:
        return False, ["R environment check failed"]
    missing = [x.strip() for x in proc.stdout.splitlines() if x.strip()]
    return len(missing) == 0, missing


def ensure_publication_environment(taxon_id: int) -> None:
    ready, missing = publication_environment_status(taxon_id)
    if ready:
        return

    if not rscript_available():
        raise PublicationEnrichmentError(
            "Rscript is not available on this host. Quick STRING enrichment remains available."
        )

    installer = _repo_root() / "scripts" / "install_publication_enrichment.R"
    if not installer.exists():
        raise PublicationEnrichmentError(
            f"Bioconductor installer not found: {installer}"
        )

    organism = ORGANISM_CODE.get(taxon_id)
    if organism is None:
        raise PublicationEnrichmentError(
            "Publication enrichment currently supports human, mouse and rat."
        )

    try:
        proc = subprocess.run(
            ["Rscript", str(installer), organism],
            capture_output=True,
            text=True,
            timeout=1200,
            check=False,
            env=_r_env(),
        )
    except subprocess.TimeoutExpired as exc:
        partial = ""
        if exc.stdout:
            partial += exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else str(exc.stdout)
        if exc.stderr:
            partial += "\n" + (
                exc.stderr.decode(errors="replace")
                if isinstance(exc.stderr, bytes)
                else str(exc.stderr)
            )
        tail = partial.strip()[-3000:]
        detail = f" Last installer output: {tail}" if tail else ""
        raise PublicationEnrichmentError(
            "R/Bioconductor installation exceeded 20 minutes on this host and was stopped. "
            "Do not retry installation from the analysis button; use a prebuilt R environment instead."
            + detail
        ) from None
    if proc.returncode != 0:
        details = (proc.stderr or proc.stdout or "Unknown R installation error").strip()
        raise PublicationEnrichmentError(
            "Could not prepare the R/Bioconductor environment. " + details[-5000:]
        )

    ready, still_missing = publication_environment_status(taxon_id)
    if not ready:
        raise PublicationEnrichmentError(
            "R/Bioconductor setup completed but required packages are still missing: "
            + ", ".join(still_missing or missing)
        )


def _run_remote_publication_enrichment(
    base_url: str,
    targets: list[str],
    taxon_id: int,
    background_mode: str,
    custom_background: list[str] | None,
) -> PublicationEnrichmentResult:
    payload = json.dumps(
        {
            "targets": targets,
            "taxon_id": taxon_id,
            "background_mode": background_mode,
            "custom_background": custom_background,
        }
    ).encode("utf-8")
    headers = {"Content-Type": "application/json"}
    backend_key = os.environ.get("PUBLICATION_BACKEND_KEY", "").strip()
    if backend_key:
        headers["X-API-Key"] = backend_key

    request = urllib.request.Request(
        base_url + "/enrich",
        data=payload,
        headers=headers,
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=660) as response:
            data = json.load(response)
    except urllib.error.HTTPError as exc:
        try:
            detail = json.loads(exc.read().decode("utf-8")).get("detail", str(exc))
        except Exception:
            detail = str(exc)
        raise PublicationEnrichmentError(
            "Publication enrichment backend failed: " + str(detail)
        ) from None
    except Exception as exc:
        raise PublicationEnrichmentError(
            "Could not reach the publication enrichment backend: " + str(exc)
        ) from None

    summary = data.get("summary")
    raw_tables = data.get("tables", {})
    if not isinstance(summary, dict):
        raise PublicationEnrichmentError("Backend response did not contain a valid summary.")

    tables: dict[str, pd.DataFrame] = {}
    for key, records in raw_tables.items():
        tables[key] = pd.DataFrame(records or [])

    return PublicationEnrichmentResult(
        outdir=Path("."),
        summary=summary,
        tables=tables,
    )


def _write_lines(path: Path, values: list[str]) -> None:
    path.write_text("\n".join(str(v).strip() for v in values if str(v).strip()) + "\n")


def run_publication_enrichment(
    targets: list[str],
    taxon_id: int,
    background_mode: str,
    custom_background: list[str] | None = None,
    workdir: str | Path | None = None,
    auto_install: bool = True,
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

    backend = publication_backend_url()
    if backend:
        return _run_remote_publication_enrichment(
            backend,
            targets,
            taxon_id,
            background_mode,
            custom_background,
        )

    if auto_install:
        ensure_publication_environment(taxon_id)
    else:
        ready, missing = publication_environment_status(taxon_id)
        if not ready:
            raise PublicationEnrichmentError(
                "Publication enrichment environment is not ready: "
                + ", ".join(missing)
            )

    r_script = _repo_root() / "r" / "enrichment_publication.R"
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

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=900,
            check=False,
            env=_r_env(),
        )
    except subprocess.TimeoutExpired:
        raise PublicationEnrichmentError(
            "Publication enrichment exceeded 15 minutes and was stopped."
        ) from None
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
