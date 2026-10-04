from __future__ import annotations

import csv
import json
import subprocess
import tempfile
from pathlib import Path
from typing import Literal

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

APP_ROOT = Path("/app")
R_SCRIPT = APP_ROOT / "r" / "enrichment_publication.R"

ORGANISM_CODE = {
    9606: "human",
    10090: "mouse",
    10116: "rat",
}

class EnrichmentRequest(BaseModel):
    targets: list[str] = Field(min_length=2, max_length=1000)
    taxon_id: int
    background_mode: Literal["default", "annotated", "custom"]
    custom_background: list[str] | None = None

app = FastAPI(
    title="Network Pharmacology Publication Enrichment API",
    version="1.0.0",
)

def _r_version_report() -> dict:
    expr = r'''
cat(jsonlite::toJSON(list(
  R=R.version.string,
  clusterProfiler=as.character(packageVersion("clusterProfiler")),
  ReactomePA=as.character(packageVersion("ReactomePA")),
  AnnotationDbi=as.character(packageVersion("AnnotationDbi")),
  GOSemSim=as.character(packageVersion("GOSemSim")),
  org_Hs=as.character(packageVersion("org.Hs.eg.db")),
  org_Mm=as.character(packageVersion("org.Mm.eg.db")),
  org_Rn=as.character(packageVersion("org.Rn.eg.db"))
), auto_unbox=TRUE))
'''
    proc = subprocess.run(
        ["Rscript", "-e", expr],
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError((proc.stderr or proc.stdout).strip())
    return json.loads(proc.stdout)

@app.get("/health")
def health() -> dict:
    try:
        versions = _r_version_report()
    except Exception as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    return {"status": "ok", "versions": versions}

@app.post("/enrich")
def enrich(payload: EnrichmentRequest) -> dict:
    organism = ORGANISM_CODE.get(payload.taxon_id)
    if organism is None:
        raise HTTPException(status_code=400, detail="Supported taxon IDs: 9606, 10090, 10116")

    targets = list(dict.fromkeys(x.strip() for x in payload.targets if x.strip()))
    if len(targets) < 2:
        raise HTTPException(status_code=400, detail="At least two non-empty target symbols are required.")

    custom_background = None
    if payload.background_mode == "custom":
        custom_background = list(
            dict.fromkeys(x.strip() for x in (payload.custom_background or []) if x.strip())
        )
        if not custom_background:
            raise HTTPException(status_code=400, detail="Custom background mode requires background genes.")

    with tempfile.TemporaryDirectory(prefix="publication_api_") as tmp:
        outdir = Path(tmp)
        fg = outdir / "foreground_symbols.txt"
        fg.write_text("\n".join(targets) + "\n")

        if payload.background_mode == "custom":
            bg = outdir / "custom_background_symbols.txt"
            bg.write_text("\n".join(custom_background or []) + "\n")
            bg_arg = str(bg)
        else:
            bg_arg = "NONE"

        cmd = [
            "Rscript",
            str(R_SCRIPT),
            str(fg),
            str(outdir),
            organism,
            payload.background_mode,
            bg_arg,
        ]
        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            raise HTTPException(status_code=504, detail="Enrichment exceeded 10 minutes.") from exc

        if proc.returncode != 0:
            detail = (proc.stderr or proc.stdout or "Unknown R error").strip()[-4000:]
            raise HTTPException(status_code=500, detail=detail)

        summary = json.loads((outdir / "summary.json").read_text())
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
            "background_universe": "background_universe.csv",
            "mapping_background": "mapping_background.csv",
        }
        tables = {}
        for key, filename in table_files.items():
            path = outdir / filename
            if not path.exists():
                continue
            with path.open(newline="") as fh:
                tables[key] = list(csv.DictReader(fh))

        return {"summary": summary, "tables": tables}
