from __future__ import annotations

import json
import os
import subprocess
import tempfile
import threading
from pathlib import Path
from typing import Literal

import pandas as pd
from fastapi import Depends, FastAPI, Header, HTTPException
from pydantic import BaseModel, Field


APP_ROOT = Path("/app")
R_SCRIPT = APP_ROOT / "r" / "enrichment_publication.R"
API_TOKEN = os.getenv("BACKEND_API_TOKEN", "").strip()
ANALYSIS_TIMEOUT = int(os.getenv("ANALYSIS_TIMEOUT_SECONDS", "1200"))
MAX_TARGETS = int(os.getenv("MAX_TARGETS", "500"))
MAX_BACKGROUND = int(os.getenv("MAX_BACKGROUND", "100000"))

ORGANISM_CODE = {
    "human": "human",
    "mouse": "mouse",
    "rat": "rat",
}

TABLE_FILES = {
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

_analysis_lock = threading.Lock()


class EnrichmentRequest(BaseModel):
    targets: list[str] = Field(min_length=2)
    organism: Literal["human", "mouse", "rat"] = "human"
    background_mode: Literal["default", "annotated", "custom"] = "annotated"
    custom_background: list[str] | None = None


def _auth(authorization: str | None = Header(default=None)) -> None:
    if not API_TOKEN:
        return
    expected = f"Bearer {API_TOKEN}"
    if authorization != expected:
        raise HTTPException(status_code=401, detail="Unauthorized")


def _dedupe(values: list[str] | None) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for value in values or []:
        value = str(value).strip()
        if not value:
            continue
        key = value.upper()
        if key not in seen:
            seen.add(key)
            out.append(value)
    return out


def _write_lines(path: Path, values: list[str]) -> None:
    path.write_text("\n".join(values) + "\n", encoding="utf-8")


def _package_versions() -> dict:
    expr = r'''
pkgs <- c(
  "clusterProfiler", "ReactomePA", "AnnotationDbi", "GOSemSim",
  "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db"
)
for (p in pkgs) {
  cat(p, "=", as.character(utils::packageVersion(p)), "\n", sep="")
}
cat("R=", R.version.string, "\n", sep="")
'''
    proc = subprocess.run(
        ["Rscript", "-e", expr],
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    if proc.returncode != 0:
        return {"ready": False, "error": (proc.stderr or proc.stdout)[-2000:]}

    versions: dict[str, str | bool] = {"ready": True}
    for line in proc.stdout.splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            versions[key.strip()] = value.strip()
    return versions


VERSIONS = _package_versions()

app = FastAPI(
    title="Network Pharmacology Publication Enrichment Backend",
    version="1.0.0",
)


@app.get("/health")
def health() -> dict:
    return {
        "status": "ok" if VERSIONS.get("ready") else "not_ready",
        "engine": "R/Bioconductor",
        "versions": VERSIONS,
    }


@app.post("/v1/enrichment", dependencies=[Depends(_auth)])
def enrichment(request: EnrichmentRequest) -> dict:
    if not VERSIONS.get("ready"):
        raise HTTPException(status_code=503, detail="R/Bioconductor engine is not ready.")

    targets = _dedupe(request.targets)
    if len(targets) < 2:
        raise HTTPException(status_code=422, detail="At least two unique target genes are required.")
    if len(targets) > MAX_TARGETS:
        raise HTTPException(
            status_code=422,
            detail=f"Maximum foreground size is {MAX_TARGETS} genes.",
        )

    custom_background = _dedupe(request.custom_background)
    if request.background_mode == "custom":
        if not custom_background:
            raise HTTPException(
                status_code=422,
                detail="Custom background mode requires a non-empty background.",
            )
        if len(custom_background) > MAX_BACKGROUND:
            raise HTTPException(
                status_code=422,
                detail=f"Maximum custom background size is {MAX_BACKGROUND} genes.",
            )

    if not _analysis_lock.acquire(blocking=False):
        raise HTTPException(
            status_code=429,
            detail="The publication enrichment engine is currently processing another analysis.",
        )

    try:
        with tempfile.TemporaryDirectory(prefix="publication_backend_") as tmp:
            outdir = Path(tmp)
            foreground_file = outdir / "foreground_symbols.txt"
            _write_lines(foreground_file, targets)

            if request.background_mode == "custom":
                background_file = outdir / "custom_background_symbols.txt"
                _write_lines(background_file, custom_background)
                background_arg = str(background_file)
            else:
                background_arg = "NONE"

            cmd = [
                "Rscript",
                str(R_SCRIPT),
                str(foreground_file),
                str(outdir),
                ORGANISM_CODE[request.organism],
                request.background_mode,
                background_arg,
            ]
            try:
                proc = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=ANALYSIS_TIMEOUT,
                    check=False,
                )
            except subprocess.TimeoutExpired:
                raise HTTPException(
                    status_code=504,
                    detail=f"Publication enrichment exceeded {ANALYSIS_TIMEOUT} seconds.",
                ) from None

            if proc.returncode != 0:
                detail = (proc.stderr or proc.stdout or "Unknown R error").strip()[-5000:]
                raise HTTPException(
                    status_code=500,
                    detail="R enrichment failed: " + detail,
                )

            summary_path = outdir / "summary.json"
            if not summary_path.exists():
                raise HTTPException(
                    status_code=500,
                    detail="R completed without producing summary.json.",
                )

            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            tables: dict[str, list[dict]] = {}

            for key, filename in TABLE_FILES.items():
                path = outdir / filename
                if not path.exists():
                    continue
                try:
                    df = pd.read_csv(path)
                except pd.errors.EmptyDataError:
                    df = pd.DataFrame()
                tables[key] = json.loads(df.to_json(orient="records"))

            for optional_key, filename in {
                "background_universe": "background_universe.csv",
                "mapping_background": "mapping_background.csv",
            }.items():
                path = outdir / filename
                if path.exists():
                    try:
                        df = pd.read_csv(path)
                    except pd.errors.EmptyDataError:
                        df = pd.DataFrame()
                    tables[optional_key] = json.loads(df.to_json(orient="records"))

            return {
                "status": "ok",
                "summary": summary,
                "tables": tables,
            }
    finally:
        _analysis_lock.release()
