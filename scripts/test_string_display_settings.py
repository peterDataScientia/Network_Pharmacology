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
from modules.string_settings_ui import string_settings_signature


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

    options = {
        "species": 9606,
        "string_version": "12.0",
        "network_type": "functional",
        "network_flavor": "evidence",
        "active_sources": [
            "textmining",
            "experiments",
            "databases",
            "coexpression",
            "neighborhood",
            "fusion",
            "cooccurrence",
        ],
        "required_score": 900,
        "first_shell": 0,
        "second_shell": 0,
        "layout": "force_directed",
        "colorblind_friendly": True,
        "bubble_3d": True,
        "block_structure_pics": False,
        "center_node_labels": False,
        "show_query_node_labels": False,
        "hide_disconnected_nodes": False,
        "hide_node_labels": False,
        "label_font_size": 12,
        "typed_physical_edges": True,
        "typed_regulatory_edges": True,
        "show_regulatory_signs": True,
    }
    baseline = string_settings_signature(options)
    advanced_changes = {
        "layout": "circular",
        "colorblind_friendly": False,
        "bubble_3d": False,
        "block_structure_pics": True,
        "center_node_labels": True,
        "show_query_node_labels": True,
        "hide_disconnected_nodes": True,
        "hide_node_labels": True,
        "label_font_size": 18,
        "typed_physical_edges": False,
        "typed_regulatory_edges": False,
        "show_regulatory_signs": False,
    }
    for key, changed_value in advanced_changes.items():
        changed = dict(options)
        changed[key] = changed_value
        assert string_settings_signature(changed) != baseline, key

    print("STRING_DISPLAY_SETTINGS_OK")


if __name__ == "__main__":
    main()
