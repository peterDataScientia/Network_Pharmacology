from __future__ import annotations

from math import prod
from typing import Iterable

STRING_SCORE_PRIOR = 0.041

SOURCE_LABELS = {
    "textmining": "Textmining",
    "experiments": "Experiments",
    "databases": "Databases",
    "coexpression": "Co-expression",
    "neighborhood": "Neighborhood",
    "fusion": "Gene Fusion",
    "cooccurrence": "Co-occurrence",
}

FUNCTIONAL_SOURCES = tuple(SOURCE_LABELS)
PHYSICAL_SOURCES = ("textmining", "experiments", "databases")
REGULATORY_SOURCES = ("textmining", "experiments", "databases")

STANDARD_FIELDS = {
    "textmining": "tscore",
    "experiments": "escore",
    "databases": "dscore",
    "coexpression": "ascore",
    "neighborhood": "nscore",
    "fusion": "fscore",
    "cooccurrence": "pscore",
}

TYPED_FUNCTIONAL_FIELDS = {
    "textmining": "functional_textmining_score",
    "experiments": "functional_experimental_score",
    "databases": "functional_database_score",
    "coexpression": "functional_coexpression_score",
    "neighborhood": "functional_neighborhood_on_chromosome_score",
    "fusion": "functional_gene_fusion_score",
    "cooccurrence": "functional_phylogenetic_cooccurrence_score",
}

REGULATORY_FIELDS = {
    "textmining": "textmining_score",
    "experiments": "experimental_score",
    "databases": "database_score",
}


def applicable_sources(network_type: str) -> tuple[str, ...]:
    if network_type == "regulatory":
        return REGULATORY_SOURCES
    if network_type == "physical":
        return PHYSICAL_SOURCES
    return FUNCTIONAL_SOURCES


def _float_score(value) -> float:
    try:
        score = float(value)
    except (TypeError, ValueError):
        return 0.0
    if score < 0:
        return 0.0
    if score > 1:
        score /= 1000.0
    return min(score, 1.0)


def combine_channel_scores(scores: Iterable[float], prior: float = STRING_SCORE_PRIOR) -> float:
    """Combine STRING evidence-channel probabilities using STRING's documented rule."""
    adjusted = []
    for raw in scores:
        score = _float_score(raw)
        if score <= prior:
            continue
        adjusted.append(max(0.0, min(1.0, (score - prior) / (1.0 - prior))))

    if not adjusted:
        return 0.0

    combined_adjusted = 1.0 - prod(1.0 - value for value in adjusted)
    return combined_adjusted + prior * (1.0 - combined_adjusted)


def _score_fields(network_type: str, network_flavor: str) -> dict[str, str]:
    if network_type == "regulatory":
        return REGULATORY_FIELDS
    if network_type == "functional" and network_flavor == "typed":
        return TYPED_FUNCTIONAL_FIELDS
    return STANDARD_FIELDS


def _official_score(row: dict, network_type: str, network_flavor: str) -> float:
    if network_type == "regulatory":
        return _float_score(row.get("combined_score"))
    if network_type == "functional" and network_flavor == "typed":
        return _float_score(row.get("functional_combined_score"))
    return _float_score(row.get("score"))


def _normalize_endpoints(row: dict, network_type: str) -> dict:
    out = dict(row)
    if network_type == "regulatory":
        out.setdefault("stringId_A", row.get("source_string_id", ""))
        out.setdefault("stringId_B", row.get("target_string_id", ""))
        out.setdefault("preferredName_A", row.get("source_preferred_name", ""))
        out.setdefault("preferredName_B", row.get("target_preferred_name", ""))
    return out


def filter_and_normalize_rows(
    rows: list[dict],
    *,
    active_sources: Iterable[str] | None,
    required_score: int,
    network_type: str,
    network_flavor: str = "evidence",
) -> list[dict]:
    """Apply user-selected STRING evidence channels and normalize endpoint names.

    STRING's public network API exposes channel scores but not the website's
    direct-vs-transferred evidence split. This function therefore filters only
    documented evidence channels and leaves evidence transfer unchanged.
    """
    applicable = applicable_sources(network_type)
    selected = tuple(
        source
        for source in (active_sources or applicable)
        if source in applicable
    )
    if not selected:
        return []

    all_applicable_selected = set(selected) == set(applicable)
    fields = _score_fields(network_type, network_flavor)
    threshold = max(0.0, min(float(required_score) / 1000.0, 1.0))
    out_rows: list[dict] = []

    for original in rows:
        row = _normalize_endpoints(original, network_type)
        official = _official_score(row, network_type, network_flavor)
        if all_applicable_selected:
            active_score = official
        else:
            active_score = combine_channel_scores(
                row.get(fields[source], 0.0)
                for source in selected
                if source in fields
            )

        if active_score + 1e-12 < threshold:
            continue

        row["analysis_score"] = active_score
        row["active_source_count"] = len(selected)
        if network_type == "regulatory":
            row.setdefault("score", official)
        elif network_type == "functional" and network_flavor == "typed":
            row.setdefault("score", official)
        out_rows.append(row)

    return out_rows
