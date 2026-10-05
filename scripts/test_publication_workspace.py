from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from modules.publication_plots import publication_enrichment_figure
from modules.publication_settings import (
    DEFAULT_PUBLICATION_ANALYSIS_SETTINGS,
    PublicationSettingsError,
    is_reference_default_settings,
    normalize_publication_analysis_settings,
)
from modules.publication_ui import _filter_by_cutoff


def main() -> None:
    defaults = normalize_publication_analysis_settings(None)
    assert defaults == DEFAULT_PUBLICATION_ANALYSIS_SETTINGS
    assert is_reference_default_settings(defaults)

    custom = normalize_publication_analysis_settings(
        {
            "ontologies": ["BP", "MF"],
            "include_reactome": False,
            "p_adjust_method": "bonferroni",
            "p_adjust_cutoff": 0.01,
            "min_gs_size": 5,
            "max_gs_size": 300,
            "reduce_go": True,
            "go_similarity_cutoff": 0.80,
            "go_bp_preselect_n": None,
            "reduce_reactome": False,
            "reactome_jaccard_cutoff": 0.75,
        }
    )
    assert custom["ontologies"] == ["BP", "MF"]
    assert custom["include_reactome"] is False
    assert custom["p_adjust_method"] == "bonferroni"
    assert custom["go_bp_preselect_n"] is None

    try:
        normalize_publication_analysis_settings(
            {"ontologies": [], "include_reactome": False}
        )
    except PublicationSettingsError:
        pass
    else:
        raise AssertionError("Expected an error when no enrichment source is selected.")

    tested = pd.DataFrame(
        {
            "ID": ["A", "B", "C"],
            "Description": ["Alpha", "Beta", "Gamma"],
            "p.adjust": [0.001, 0.02, 0.08],
            "Count": [8, 5, 4],
            "GeneRatio": ["8/100", "5/100", "4/100"],
        }
    )
    strict = _filter_by_cutoff(tested, 0.01)
    assert strict["ID"].tolist() == ["A"]

    tables = {
        "go_bp_reduced": tested.iloc[:2].copy(),
        "reactome_reduced": pd.DataFrame(
            {
                "ID": ["R1"],
                "Description": ["Reactome one"],
                "p.adjust": [0.005],
                "Count": [6],
                "GeneRatio": ["6/100"],
            }
        ),
    }
    for plot_type in ("Mirrored", "Dot plot", "Horizontal bar"):
        fig = publication_enrichment_figure(
            tables,
            plot_type=plot_type,
            categories=["GO-BP", "Reactome"],
            top_n_by_category={"GO-BP": 2, "Reactome": 1},
            sort_by="Adjusted P-value",
            display_cutoff=0.05,
            wrap_width=40,
        )
        assert fig is not None
        fig.clear()

    print("PUBLICATION_WORKSPACE_TESTS_OK")


if __name__ == "__main__":
    main()
