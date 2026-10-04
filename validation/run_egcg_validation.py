from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from modules.publication_enrichment import run_publication_enrichment
TARGETS = [
    x.strip()
    for x in (ROOT / "validation" / "egcg_risi_shared_targets.txt").read_text().splitlines()
    if x.strip()
]

REFERENCE_COUNTS = {
    "go_bp_raw": 1432,
    "go_bp_reduced": 49,
    "go_cc_raw": 23,
    "go_cc_reduced": 17,
    "go_mf_raw": 54,
    "go_mf_reduced": 26,
    "reactome_raw": 276,
    "reactome_reduced": 91,
}


def main() -> None:
    base = ROOT / "validation" / "outputs"
    base.mkdir(parents=True, exist_ok=True)

    for mode in ("default", "annotated"):
        result = run_publication_enrichment(
            targets=TARGETS,
            taxon_id=9606,
            background_mode=mode,
            workdir=base / mode,
        )
        print("\n===", mode.upper(), "===")
        print(json.dumps(result.summary, indent=2))

        if mode == "default":
            observed = result.summary.get("counts", {})
            print("\nHistorical notebook count comparison:")
            for key, expected in REFERENCE_COUNTS.items():
                actual = observed.get(key)
                marker = "MATCH" if actual == expected else "DIFF"
                print(f"{key:20s} expected={expected:4d} observed={actual!s:>4}  {marker}")

    print(
        "\nValidation complete. Review term-level outputs and package versions before merging."
    )


if __name__ == "__main__":
    main()
