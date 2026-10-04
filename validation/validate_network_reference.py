from pathlib import Path

import pandas as pd

from modules.hub_consensus_v2 import (
    PRIMARY_METRICS,
    build_graph,
    centrality_table,
    consensus_hub_analysis,
)
from modules.reference_validation import (
    EGCG_RISI_NETWORK_EXPECTED,
    EGCG_RISI_TOP10_EXPECTED,
)

ROOT = Path(__file__).resolve().parents[1]
network = pd.read_csv(
    ROOT / "validation" / "egcg_risi_string_v12_network.tsv",
    sep="\t",
)
network["score"] = 1.0

graph = build_graph(network)
centrality = centrality_table(graph)
ranked, hubs, effective_top_n = consensus_hub_analysis(
    centrality,
    top_n=EGCG_RISI_NETWORK_EXPECTED["top_n"],
)

assert graph.number_of_nodes() == EGCG_RISI_NETWORK_EXPECTED["connected"]
assert graph.number_of_edges() == EGCG_RISI_NETWORK_EXPECTED["edges"]
assert effective_top_n == 10
assert set(hubs["Gene"]) == set(EGCG_RISI_NETWORK_EXPECTED["hubs"])

for metric in PRIMARY_METRICS:
    observed = (
        centrality.sort_values(
            [metric, "Gene"],
            ascending=[False, True],
            na_position="last",
        )
        .head(10)["Gene"]
        .tolist()
    )
    expected = list(EGCG_RISI_TOP10_EXPECTED[metric])
    assert observed == expected, f"{metric}: {observed} != {expected}"

print("NETWORK_REFERENCE_PARITY_OK")
print("nodes=30 edges=90 hubs=AKT1,BCL2L1,CASP3,STAT3,TP53")
