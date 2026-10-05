from __future__ import annotations

import time
from typing import Iterable

import pandas as pd
import requests

from modules.string_filters import (
    applicable_sources,
    filter_and_normalize_rows,
)

# Explicit STRING version pinning is required for reproducible network topology.
# The EGCG/RISI manuscript reference workflow used STRING v12.0.
STRING_API_BASES = {
    "12.0": "https://version-12.string-db.org/api",
    "12.5": "https://version-12-5.string-db.org/api",
}
CALLER_IDENTITY = "Network_Pharmacology_Streamlit_App"


class StringAPIError(RuntimeError):
    pass


def _api_base(string_version: str) -> str:
    try:
        return STRING_API_BASES[string_version]
    except KeyError as exc:
        raise StringAPIError(
            f"Unsupported STRING version {string_version!r}. "
            f"Choose one of: {', '.join(STRING_API_BASES)}"
        ) from exc


def _post(
    method: str,
    output_format: str,
    payload: dict,
    timeout: int = 90,
    string_version: str = "12.0",
) -> requests.Response:
    url = f"{_api_base(string_version)}/{output_format}/{method}"
    try:
        response = requests.post(url, data=payload, timeout=timeout)
    except requests.RequestException as exc:
        raise StringAPIError(f"Could not contact STRING: {exc}") from exc

    if response.status_code != 200:
        msg = response.text.strip().replace("\n", " ")[:300]
        raise StringAPIError(f"STRING returned HTTP {response.status_code}: {msg}")

    return response


def _post_json(
    method: str,
    payload: dict,
    timeout: int = 60,
    string_version: str = "12.0",
):
    response = _post(
        method,
        "json",
        payload,
        timeout=timeout,
        string_version=string_version,
    )
    try:
        return response.json()
    except ValueError as exc:
        raise StringAPIError("STRING returned an unreadable JSON response.") from exc


def _post_binary(
    output_format: str,
    method: str,
    payload: dict,
    timeout: int = 90,
    string_version: str = "12.0",
) -> bytes:
    return _post(
        method,
        output_format,
        payload,
        timeout=timeout,
        string_version=string_version,
    ).content


def _post_text(
    output_format: str,
    method: str,
    payload: dict,
    timeout: int = 90,
    string_version: str = "12.0",
) -> str:
    return _post(
        method,
        output_format,
        payload,
        timeout=timeout,
        string_version=string_version,
    ).text


def _identifiers(values: Iterable[str]) -> str:
    return "\r".join(str(v).strip() for v in values if str(v).strip())


def map_identifiers(
    targets: list[str],
    species: int,
    string_version: str = "12.0",
) -> pd.DataFrame:
    payload = {
        "identifiers": _identifiers(targets),
        "species": species,
        "echo_query": 1,
        "caller_identity": CALLER_IDENTITY,
    }
    rows = _post_json(
        "get_string_ids",
        payload,
        string_version=string_version,
    )
    if not rows:
        return pd.DataFrame(columns=["queryItem", "stringId", "preferredName", "annotation"])
    return pd.DataFrame(rows)


def get_network(
    string_ids: list[str],
    species: int,
    required_score: int,
    network_type: str,
    *,
    network_flavor: str = "evidence",
    active_sources: list[str] | None = None,
    add_nodes: int = 0,
    typed_physical_edges: bool = True,
    typed_regulatory_edges: bool = True,
    show_query_node_labels: bool = False,
    string_version: str = "12.0",
) -> pd.DataFrame:
    payload = {
        "identifiers": _identifiers(string_ids),
        "species": species,
        "required_score": int(required_score),
        "network_type": network_type,
        "add_nodes": max(0, int(add_nodes)),
        "show_query_node_labels": int(bool(show_query_node_labels)),
        "caller_identity": CALLER_IDENTITY,
    }
    if network_type == "functional" and network_flavor == "typed":
        payload.update(
            network_flavor="typed",
            typed_physical_edges=int(bool(typed_physical_edges)),
            typed_regulatory_edges=int(bool(typed_regulatory_edges)),
        )

    rows = _post_json(
        "network",
        payload,
        string_version=string_version,
    )
    selected_sources = active_sources or list(applicable_sources(network_type))
    filtered = filter_and_normalize_rows(
        rows if isinstance(rows, list) else [],
        active_sources=selected_sources,
        required_score=required_score,
        network_type=network_type,
        network_flavor=network_flavor,
    )
    return pd.DataFrame(filtered)


def get_enrichment(
    string_ids: list[str],
    species: int,
    string_version: str = "12.0",
) -> pd.DataFrame:
    payload = {
        "identifiers": _identifiers(string_ids),
        "species": species,
        "caller_identity": CALLER_IDENTITY,
    }
    rows = _post_json(
        "enrichment",
        payload,
        string_version=string_version,
    )
    return pd.DataFrame(rows)


def get_network_media(
    string_ids: list[str],
    species: int,
    required_score: int,
    network_type: str,
    *,
    network_flavor: str = "evidence",
    first_shell: int = 0,
    second_shell: int = 0,
    active_sources: list[str] | None = None,
    typed_physical_edges: bool = True,
    typed_regulatory_edges: bool = True,
    show_regulatory_signs: bool = True,
    bubble_3d: bool = True,
    block_structure_pics: bool = False,
    center_node_labels: bool = False,
    show_query_node_labels: bool = False,
    hide_disconnected_nodes: bool = False,
    hide_node_labels: bool = False,
    label_font_size: int = 12,
    string_version: str = "12.0",
) -> dict:
    """Retrieve STRING's own high-resolution PNG, SVG and stable network link.

    The public STRING image/link API exposes network type, flavor, neighborhood
    and display controls. It does not expose the website's evidence-channel or
    direct-vs-transferred evidence toggles, so source-filtered analyses are
    accompanied by a provenance warning rather than a misleading native image claim.
    """
    image_payload = {
        "identifiers": _identifiers(string_ids),
        "species": species,
        "required_score": int(required_score),
        "network_type": network_type,
        "network_flavor": network_flavor,
        "add_color_nodes": max(0, int(first_shell)),
        "add_white_nodes": max(0, int(second_shell)),
        "typed_physical_edges": int(bool(typed_physical_edges)),
        "typed_regulatory_edges": int(bool(typed_regulatory_edges)),
        "show_regulatory_signs": int(bool(show_regulatory_signs)),
        "hide_node_labels": int(bool(hide_node_labels)),
        "hide_disconnected_nodes": int(bool(hide_disconnected_nodes)),
        "show_query_node_labels": int(bool(show_query_node_labels)),
        "block_structure_pics_in_bubbles": int(bool(block_structure_pics)),
        "flat_node_design": int(not bool(bubble_3d)),
        "center_node_labels": int(bool(center_node_labels)),
        "custom_label_font_size": max(5, min(50, int(label_font_size))),
        "caller_identity": CALLER_IDENTITY,
    }

    link_payload = {
        key: value
        for key, value in image_payload.items()
        if key
        not in {
            "flat_node_design",
            "center_node_labels",
            "custom_label_font_size",
        }
    }

    media = {
        "highres_png": None,
        "svg": None,
        "link": None,
        "errors": [],
        "warnings": [],
    }

    selected = set(active_sources or applicable_sources(network_type))
    if selected != set(applicable_sources(network_type)):
        media["warnings"].append(
            "The app analysis is filtered to the selected evidence channels. "
            "STRING's public image/link API does not expose that channel filter, "
            "so the native STRING figure/link may contain additional all-source edges. "
            "Use the app-generated network as the authoritative filtered topology."
        )

    try:
        media["highres_png"] = _post_binary(
            "highres_image",
            "network",
            image_payload,
            string_version=string_version,
        )
    except StringAPIError as exc:
        media["errors"].append(f"High-resolution PNG: {exc}")

    time.sleep(1.0)

    try:
        svg_text = _post_text(
            "svg",
            "network",
            image_payload,
            string_version=string_version,
        )
        media["svg"] = svg_text.encode("utf-8")
    except StringAPIError as exc:
        media["errors"].append(f"SVG: {exc}")

    time.sleep(1.0)

    try:
        link = _post_text(
            "tsv-no-header",
            "get_link",
            link_payload,
            string_version=string_version,
        ).strip()
        if link:
            fields = link.split("\t")
            media["link"] = next(
                (field.strip() for field in fields if field.strip().startswith("http")),
                link,
            )
    except StringAPIError as exc:
        media["errors"].append(f"STRING link: {exc}")

    return media


def run_string_workflow(
    targets: list[str],
    species: int,
    required_score: int,
    network_type: str,
    *,
    network_flavor: str = "evidence",
    active_sources: list[str] | None = None,
    first_shell: int = 0,
    second_shell: int = 0,
    typed_physical_edges: bool = True,
    typed_regulatory_edges: bool = True,
    show_regulatory_signs: bool = True,
    bubble_3d: bool = True,
    block_structure_pics: bool = False,
    center_node_labels: bool = False,
    show_query_node_labels: bool = False,
    hide_disconnected_nodes: bool = False,
    hide_node_labels: bool = False,
    label_font_size: int = 12,
    string_version: str = "12.0",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """Run the complete STRING workflow with courteous sequential API calls."""
    mapping = map_identifiers(
        targets,
        species,
        string_version=string_version,
    )
    if mapping.empty:
        raise StringAPIError("None of the submitted identifiers could be mapped by STRING.")

    ids = mapping["stringId"].dropna().astype(str).drop_duplicates().tolist()
    selected_sources = active_sources or list(applicable_sources(network_type))

    time.sleep(1.0)
    network = get_network(
        ids,
        species,
        required_score,
        network_type,
        network_flavor=network_flavor,
        active_sources=selected_sources,
        add_nodes=max(0, int(first_shell) + int(second_shell)),
        typed_physical_edges=typed_physical_edges,
        typed_regulatory_edges=typed_regulatory_edges,
        show_query_node_labels=show_query_node_labels,
        string_version=string_version,
    )

    time.sleep(1.0)
    enrichment = (
        get_enrichment(ids, species, string_version=string_version)
        if len(ids) >= 2
        else pd.DataFrame()
    )

    time.sleep(1.0)
    native_media = get_network_media(
        ids,
        species,
        required_score,
        network_type,
        network_flavor=network_flavor,
        first_shell=first_shell,
        second_shell=second_shell,
        active_sources=selected_sources,
        typed_physical_edges=typed_physical_edges,
        typed_regulatory_edges=typed_regulatory_edges,
        show_regulatory_signs=show_regulatory_signs,
        bubble_3d=bubble_3d,
        block_structure_pics=block_structure_pics,
        center_node_labels=center_node_labels,
        show_query_node_labels=show_query_node_labels,
        hide_disconnected_nodes=hide_disconnected_nodes,
        hide_node_labels=hide_node_labels,
        label_font_size=label_font_size,
        string_version=string_version,
    )

    return mapping, network, enrichment, native_media
