from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from modules.publication_settings import (
    is_reference_default_settings,
    normalize_publication_analysis_settings,
)
from modules.reference_validation import (
    EGCG_RISI_ENRICHMENT_DEFAULT_EXPECTED,
    PUBLICATION_VERSION_EXPECTED,
    is_egcg_risi_reference,
)


def fail(message: str) -> None:
    raise SystemExit(message)


def main() -> None:
    request_path = Path(sys.argv[1] if len(sys.argv) > 1 else "job_input/request.json")
    summary_path = Path(sys.argv[2] if len(sys.argv) > 2 else "job_output/summary.json")

    request = json.loads(request_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))

    if summary.get("status") != "ok":
        fail("summary.json did not report status=ok")

    expected_mode = request["background_mode"]
    if summary.get("background_mode") != expected_mode:
        fail(
            f"Background mode mismatch: observed {summary.get('background_mode')!r}, "
            f"expected {expected_mode!r}"
        )

    expected_settings = normalize_publication_analysis_settings(
        request.get("analysis_settings")
    )
    observed_settings = summary.get("analysis_settings", {})
    observed_normalized = normalize_publication_analysis_settings({
        **observed_settings,
        "go_bp_preselect_n": (
            None
            if observed_settings.get("go_bp_preselect_n") == "all"
            else observed_settings.get("go_bp_preselect_n")
        ),
    })
    if observed_normalized != expected_settings:
        fail(
            "Analysis settings mismatch between request and R output: "
            + json.dumps(
                {"observed": observed_normalized, "expected": expected_settings},
                sort_keys=True,
            )
        )

    expected_submitted = len(request["targets"])
    if summary.get("foreground_symbols_submitted") != expected_submitted:
        fail(
            "Foreground count mismatch: "
            f"observed {summary.get('foreground_symbols_submitted')}, "
            f"expected {expected_submitted}"
        )

    versions = summary.get("versions", {})
    r_string = str(versions.get("R", ""))
    if PUBLICATION_VERSION_EXPECTED["R"] not in r_string:
        fail(f"Unexpected R version: {r_string!r}")

    for package in ("clusterProfiler", "ReactomePA", "AnnotationDbi", "GOSemSim"):
        observed = str(versions.get(package, ""))
        expected = PUBLICATION_VERSION_EXPECTED[package]
        if observed != expected:
            fail(f"{package} version mismatch: observed {observed}, expected {expected}")

    expected_orgdb = PUBLICATION_VERSION_EXPECTED["organism_db"][request["taxon_id"]]
    observed_orgdb = str(versions.get("organism_db", ""))
    if observed_orgdb != expected_orgdb:
        fail(
            f"Organism DB version mismatch: observed {observed_orgdb}, "
            f"expected {expected_orgdb}"
        )

    reference_parity = None
    if (
        request["taxon_id"] == 9606
        and expected_mode == "default"
        and is_egcg_risi_reference(request["targets"])
        and is_reference_default_settings(expected_settings)
    ):
        observed_counts = summary.get("counts", {})
        mismatches = {
            key: {"observed": observed_counts.get(key), "expected": expected}
            for key, expected in EGCG_RISI_ENRICHMENT_DEFAULT_EXPECTED.items()
            if observed_counts.get(key) != expected
        }
        if mismatches:
            fail("EGCG/RISI reference parity failed: " + json.dumps(mismatches, sort_keys=True))
        reference_parity = "passed"
        print("EGCG_RISI_REFERENCE_PARITY_OK")

    metadata = {
        "request_id": request["request_id"],
        "executor": "github-actions",
        "workflow_run_id": os.environ.get("GITHUB_RUN_ID"),
        "workflow_run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
        "git_sha": os.environ.get("GITHUB_SHA"),
        "validated_at_utc": datetime.now(timezone.utc).isoformat(),
        "reference_parity": reference_parity,
        "analysis_settings": expected_settings,
        "versions": versions,
    }
    out = summary_path.parent / "job_metadata.json"
    out.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print("PUBLICATION_JOB_VALIDATION_OK")


if __name__ == "__main__":
    main()
