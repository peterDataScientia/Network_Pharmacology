from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from modules.hub_consensus_v2 import build_display_graph, build_graph
from modules.hub_plots_v2 import network_figure


def main() -> None:
    mapping = pd.DataFrame(
        {
            "queryItem": ["TP53", "MDM2", "AKT1"],
            "preferredName": ["TP53", "MDM2", "AKT1"],
        }
    )
    network = pd.DataFrame(
        {
            "preferredName_A": ["TP53"],
            "preferredName_B": ["MDM2"],
            "score": [0.95],
        }
    )

    topology = build_graph(network)
    assert set(topology.nodes) == {"TP53", "MDM2"}

    shown = build_display_graph(
        network,
        mapping,
        hide_disconnected_nodes=False,
    )
    assert set(shown.nodes) == {"TP53", "MDM2", "AKT1"}

    hidden = build_display_graph(
        network,
        mapping,
        hide_disconnected_nodes=True,
    )
    assert set(hidden.nodes) == {"TP53", "MDM2"}

    query_labels = build_display_graph(
        network,
        pd.DataFrame(
            {
                "queryItem": ["P53_INPUT", "MDM2_INPUT", "AKT1_INPUT"],
                "preferredName": ["TP53", "MDM2", "AKT1"],
            }
        ),
        hide_disconnected_nodes=False,
        use_query_labels=True,
    )
    assert "AKT1_INPUT" in query_labels.nodes

    regulatory = pd.DataFrame(
        {
            "source_preferred_name": ["TP53", "MDM2"],
            "target_preferred_name": ["MDM2", "TP53"],
            "combined_score": [0.91, 0.88],
            "sign": ["pos", "neg"],
        }
    )
    directed = build_graph(regulatory)
    assert directed.is_directed()
    assert directed["TP53"]["MDM2"]["sign"] == "pos"
    assert directed["MDM2"]["TP53"]["sign"] == "neg"

    for centered in (False, True):
        for signs in (False, True):
            fig = network_figure(
                directed,
                layout="circular",
                show_labels=True,
                label_font_size=12,
                colorblind_friendly=True,
                center_node_labels=centered,
                show_regulatory_signs=signs,
            )
            assert fig is not None
            plt.close(fig)

    print("STRING_DISPLAY_SETTINGS_OK")


if __name__ == "__main__":
    main()
