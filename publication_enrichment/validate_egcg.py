from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from modules.publication_enrichment import (  # noqa: E402
    PublicationEnrichmentError,
    run_publication_enrichment,
)


EXPECTED_LEGACY_COUNTS = {
    "significant_counts": {
        "go_bp": 1432,
        "go_cc": 23,
        "go_mf": 54,
        "reactome": 276,
    },
    "nonredundant_counts": {
        "go_bp": 49,
        "go_cc": 17,
        "go_mf": 26,
        "reactome": 91,
    },
}


def compare_counts(summary: dict) -> list[str]:
    mismatches = []
    for section, expected_values in EXPECTED_LEGACY_COUNTS.items():
        observed = summary.get(section, {})
        for key, expected in expected_values.items():
            actual = observed.get(key)
            if actual != expected:
                mismatches.append(
                    f"{section}.{key}: expected {expected}, observed {actual}"
                )
    return mismatches


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Regression-test the publication enrichment backend against the "
            "executed EGCG notebook counts."
        )
    )
    parser.add_argument(
        "--allow-annotation-drift",
        action="store_true",
        help=(
            "Report count mismatches without failing. Use only when package/database "
            "versions differ from the notebook environment."
        ),
    )
    args = parser.parse_args()

    target_file = (
        ROOT
        / "publication_enrichment"
        / "example"
        / "egcg_risi_shared_targets.csv"
    )
    targets = pd.read_csv(target_file)["GeneSymbol"].dropna().astype(str).tolist()

    try:
        result = run_publication_enrichment(
            targets,
            background_mode="package_default",
            fdr_cutoff=0.05,
            go_similarity_cutoff=0.70,
            reactome_similarity_cutoff=0.70,
            analysis_profile="legacy_notebook",
            timeout_seconds=1800,
        )
    except PublicationEnrichmentError as exc:
        print(f"VALIDATION ERROR: {exc}", file=sys.stderr)
        return 2

    summary = result["summary"]
    print(json.dumps(summary, indent=2))

    mismatches = compare_counts(summary)
    if mismatches:
        print("\nCOUNT MISMATCHES:")
        for mismatch in mismatches:
            print(f"  - {mismatch}")
        print(
            "\nSoftware versions:",
            json.dumps(summary.get("software_versions", {}), indent=2),
        )
        if not args.allow_annotation_drift:
            print(
                "\nFAIL: counts differ from the executed reference notebook. "
                "Check annotation/package versions before accepting the backend."
            )
            return 1
        print("\nWARNING: annotation drift allowed; structural validation only.")
    else:
        print(
            "\nPASS: significant and non-redundant counts exactly reproduce "
            "the executed EGCG reference notebook."
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
