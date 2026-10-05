from __future__ import annotations

DEFAULT_PUBLICATION_ANALYSIS_SETTINGS = {
    "ontologies": ["BP", "CC", "MF"],
    "include_reactome": True,
    "p_adjust_method": "BH",
    "p_adjust_cutoff": 0.05,
    "min_gs_size": 10,
    "max_gs_size": 500,
    "reduce_go": True,
    "go_similarity_cutoff": 0.70,
    "go_bp_preselect_n": 100,
    "reduce_reactome": True,
    "reactome_jaccard_cutoff": 0.70,
}

ALLOWED_ONTOLOGIES = ("BP", "CC", "MF")
ALLOWED_P_ADJUST_METHODS = ("BH", "bonferroni", "holm", "BY")


class PublicationSettingsError(ValueError):
    pass


def normalize_publication_analysis_settings(settings: dict | None = None) -> dict:
    out = dict(DEFAULT_PUBLICATION_ANALYSIS_SETTINGS)
    if settings:
        out.update(settings)

    ontologies = []
    for value in out.get("ontologies") or []:
        code = str(value).strip().upper()
        if code in ALLOWED_ONTOLOGIES and code not in ontologies:
            ontologies.append(code)
    out["ontologies"] = ontologies
    out["include_reactome"] = bool(out.get("include_reactome", True))

    method = str(out.get("p_adjust_method", "BH")).strip()
    aliases = {
        "bh": "BH",
        "benjamini-hochberg": "BH",
        "benjamini hochberg": "BH",
        "bonferroni": "bonferroni",
        "holm": "holm",
        "by": "BY",
        "benjamini-yekutieli": "BY",
        "benjamini yekutieli": "BY",
    }
    method = aliases.get(method.lower(), method)
    if method not in ALLOWED_P_ADJUST_METHODS:
        raise PublicationSettingsError(
            "p_adjust_method must be one of: "
            + ", ".join(ALLOWED_P_ADJUST_METHODS)
        )
    out["p_adjust_method"] = method

    try:
        cutoff = float(out.get("p_adjust_cutoff", 0.05))
    except (TypeError, ValueError) as exc:
        raise PublicationSettingsError("p_adjust_cutoff must be numeric.") from exc
    if not 0 < cutoff <= 1:
        raise PublicationSettingsError("p_adjust_cutoff must be > 0 and <= 1.")
    out["p_adjust_cutoff"] = cutoff

    try:
        min_size = int(out.get("min_gs_size", 10))
        max_size = int(out.get("max_gs_size", 500))
    except (TypeError, ValueError) as exc:
        raise PublicationSettingsError("Gene-set size limits must be integers.") from exc
    if min_size < 1:
        raise PublicationSettingsError("min_gs_size must be at least 1.")
    if max_size < min_size:
        raise PublicationSettingsError("max_gs_size must be >= min_gs_size.")
    out["min_gs_size"] = min_size
    out["max_gs_size"] = max_size

    out["reduce_go"] = bool(out.get("reduce_go", True))
    try:
        go_cutoff = float(out.get("go_similarity_cutoff", 0.70))
    except (TypeError, ValueError) as exc:
        raise PublicationSettingsError("go_similarity_cutoff must be numeric.") from exc
    if not 0 < go_cutoff <= 1:
        raise PublicationSettingsError("go_similarity_cutoff must be > 0 and <= 1.")
    out["go_similarity_cutoff"] = go_cutoff

    bp_n = out.get("go_bp_preselect_n", 100)
    if bp_n in (None, "", 0, "all", "All"):
        bp_n = None
    else:
        try:
            bp_n = int(bp_n)
        except (TypeError, ValueError) as exc:
            raise PublicationSettingsError(
                "go_bp_preselect_n must be a positive integer or null."
            ) from exc
        if bp_n < 1:
            raise PublicationSettingsError(
                "go_bp_preselect_n must be a positive integer or null."
            )
    out["go_bp_preselect_n"] = bp_n

    out["reduce_reactome"] = bool(out.get("reduce_reactome", True))
    try:
        reactome_cutoff = float(out.get("reactome_jaccard_cutoff", 0.70))
    except (TypeError, ValueError) as exc:
        raise PublicationSettingsError(
            "reactome_jaccard_cutoff must be numeric."
        ) from exc
    if not 0 < reactome_cutoff <= 1:
        raise PublicationSettingsError(
            "reactome_jaccard_cutoff must be > 0 and <= 1."
        )
    out["reactome_jaccard_cutoff"] = reactome_cutoff

    if not out["ontologies"] and not out["include_reactome"]:
        raise PublicationSettingsError(
            "Select at least one enrichment database/ontology."
        )

    return out


def is_reference_default_settings(settings: dict | None) -> bool:
    return (
        normalize_publication_analysis_settings(settings)
        == DEFAULT_PUBLICATION_ANALYSIS_SETTINGS
    )
