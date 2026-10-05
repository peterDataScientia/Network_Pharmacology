from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from modules.reference_validation import (
    EGCG_RISI_NETWORK_EXPECTED,
    EGCG_RISI_TOP10_EXPECTED,
)

outdir = Path(sys.argv[1] if len(sys.argv) > 1 else "centrality_output")

summary = json.loads((outdir / "centrality_summary.json").read_text(encoding="utf-8"))
centrality = pd.read_csv(outdir / "centrality_all_genes.csv")
ranked = pd.read_csv(outdir / "centrality_consensus_rankings.csv")
hubs = pd.read_csv(outdir / "consensus_hubs_4of4.csv")

assert summary["status"] == "ok"
assert summary["engine"]["language"] == "R"
assert summary["engine"]["package"] == "igraph"
assert summary["engine"]["igraph"] == "2.3.4"
assert "4.6.1" in summary["engine"]["R"]

assert summary["nodes"] == EGCG_RISI_NETWORK_EXPECTED["connected"]
assert summary["edges"] == EGCG_RISI_NETWORK_EXPECTED["edges"]
assert summary["top_n_effective"] == EGCG_RISI_NETWORK_EXPECTED["top_n"]
assert set(hubs["Gene"]) == set(EGCG_RISI_NETWORK_EXPECTED["hubs"])

for metric, expected in EGCG_RISI_TOP10_EXPECTED.items():
    observed = (
        centrality.sort_values(
            [metric, "Gene"],
            ascending=[False, True],
            na_position="last",
        )
        .head(10)["Gene"]
        .tolist()
    )
    assert observed == list(expected), f"{metric}: {observed} != {list(expected)}"

assert len(ranked) == EGCG_RISI_NETWORK_EXPECTED["connected"]

print("R_IGRAPH_NETWORK_REFERENCE_PARITY_OK")
print("nodes=30 edges=90 hubs=AKT1,BCL2L1,CASP3,STAT3,TP53")
