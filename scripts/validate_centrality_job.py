from __future__ import annotations

import csv
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

EXPECTED_IGRAPH = "2.3.4"
EXPECTED_R = "4.6.1"


def fail(message: str) -> None:
    raise SystemExit(message)


def _csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main() -> None:
    request_path = Path(sys.argv[1] if len(sys.argv) > 1 else "centrality_input/request.json")
    outdir = Path(sys.argv[2] if len(sys.argv) > 2 else "centrality_output")

    request = json.loads(request_path.read_text(encoding="utf-8"))
    summary_path = outdir / "centrality_summary.json"
    if not summary_path.exists():
        fail("centrality_summary.json is missing")

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "ok":
        fail("Centrality summary did not report status=ok")

    engine = summary.get("engine", {})
    if engine.get("language") != "R" or engine.get("package") != "igraph":
        fail(f"Unexpected centrality engine: {engine}")

    if EXPECTED_R not in str(engine.get("R", "")):
        fail(f"Unexpected R version: {engine.get('R')!r}")
    if str(engine.get("igraph", "")) != EXPECTED_IGRAPH:
        fail(
            f"igraph version mismatch: observed {engine.get('igraph')}, "
            f"expected {EXPECTED_IGRAPH}"
        )

    expected_nodes = int(request.get("expected_nodes", -1))
    expected_edges = int(request.get("expected_edges", -1))
    expected_top_n = min(int(request["top_n"]), max(expected_nodes, 0))

    if expected_nodes >= 0 and int(summary.get("nodes", -1)) != expected_nodes:
        fail(
            f"Node-count mismatch: observed {summary.get('nodes')}, "
            f"expected {expected_nodes}"
        )
    if expected_edges >= 0 and int(summary.get("edges", -1)) != expected_edges:
        fail(
            f"Edge-count mismatch: observed {summary.get('edges')}, "
            f"expected {expected_edges}"
        )
    if expected_nodes >= 0 and int(summary.get("top_n_effective", -1)) != expected_top_n:
        fail(
            f"Effective Top-N mismatch: observed {summary.get('top_n_effective')}, "
            f"expected {expected_top_n}"
        )

    required = [
        outdir / "centrality_all_genes.csv",
        outdir / "centrality_consensus_rankings.csv",
        outdir / "consensus_hubs_4of4.csv",
    ]
    for path in required:
        if not path.exists():
            fail(f"Missing centrality output: {path.name}")

    if expected_nodes > 0:
        centrality = _csv_rows(outdir / "centrality_all_genes.csv")
        ranked = _csv_rows(outdir / "centrality_consensus_rankings.csv")
        hubs = _csv_rows(outdir / "consensus_hubs_4of4.csv")

        if len(centrality) != expected_nodes:
            fail(
                f"Centrality table has {len(centrality)} rows, expected {expected_nodes}"
            )
        if len(ranked) != expected_nodes:
            fail(
                f"Ranked table has {len(ranked)} rows, expected {expected_nodes}"
            )
        if len(hubs) != int(summary.get("hub_count", -1)):
            fail(
                f"Hub-count mismatch: table has {len(hubs)}, "
                f"summary has {summary.get('hub_count')}"
            )

    metadata = {
        "request_id": request.get("request_id"),
        "executor": "github-actions",
        "centrality_engine": "R/igraph",
        "workflow_run_id": os.environ.get("GITHUB_RUN_ID"),
        "workflow_run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
        "git_sha": os.environ.get("GITHUB_SHA"),
        "validated_at_utc": datetime.now(timezone.utc).isoformat(),
        "engine": engine,
        "string_version": request.get("string_version"),
        "required_score": request.get("required_score"),
        "network_type": request.get("network_type"),
        "network_flavor": request.get("network_flavor"),
        "active_sources": request.get("active_sources", []),
        "add_nodes": request.get("add_nodes", 0),
        "edge_hash": request.get("observed_edge_hash"),
    }
    (outdir / "centrality_job_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )
    print("R_IGRAPH_CENTRALITY_JOB_OK")


if __name__ == "__main__":
    main()
