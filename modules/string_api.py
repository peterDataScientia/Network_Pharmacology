from __future__ import annotations

import time
from typing import Iterable

import pandas as pd
import requests

# Pin STRING for reproducible production behavior.
STRING_API_BASE = "https://version-12-5.string-db.org/api"
CALLER_IDENTITY = "Network_Pharmacology_Streamlit_App"


class StringAPIError(RuntimeError):
    pass


def _post(method: str, output_format: str, payload: dict, timeout: int = 90) -> requests.Response:
    url = f"{STRING_API_BASE}/{output_format}/{method}"
    try:
        response = requests.post(url, data=payload, timeout=timeout)
    except requests.RequestException as exc:
        raise StringAPIError(f"Could not contact STRING: {exc}") from exc

    if response.status_code != 200:
        msg = response.text.strip().replace("\n", " ")[:300]
        raise StringAPIError(f"STRING returned HTTP {response.status_code}: {msg}")

    return response


def _post_json(method: str, payload: dict, timeout: int = 60):
    response = _post(method, "json", payload, timeout=timeout)
    try:
        return response.json()
    except ValueError as exc:
        raise StringAPIError("STRING returned an unreadable JSON response.") from exc


def _post_binary(output_format: str, method: str, payload: dict, timeout: int = 90) -> bytes:
    return _post(method, output_format, payload, timeout=timeout).content


def _post_text(output_format: str, method: str, payload: dict, timeout: int = 90) -> str:
    return _post(method, output_format, payload, timeout=timeout).text


def _identifiers(values: Iterable[str]) -> str:
    return "\r".join(str(v).strip() for v in values if str(v).strip())


def map_identifiers(targets: list[str], species: int) -> pd.DataFrame:
    payload = {
        "identifiers": _identifiers(targets),
        "species": species,
        "echo_query": 1,
        "caller_identity": CALLER_IDENTITY,
    }
    rows = _post_json("get_string_ids", payload)
    if not rows:
        return pd.DataFrame(columns=["queryItem", "stringId", "preferredName", "annotation"])
    return pd.DataFrame(rows)


def get_network(
    string_ids: list[str],
    species: int,
    required_score: int,
    network_type: str,
) -> pd.DataFrame:
    payload = {
        "identifiers": _identifiers(string_ids),
        "species": species,
        "required_score": int(required_score),
        "network_type": network_type,
        "add_nodes": 0,
        "caller_identity": CALLER_IDENTITY,
    }
    rows = _post_json("network", payload)
    return pd.DataFrame(rows)


def get_enrichment(string_ids: list[str], species: int) -> pd.DataFrame:
    payload = {
        "identifiers": _identifiers(string_ids),
        "species": species,
        "caller_identity": CALLER_IDENTITY,
    }
    rows = _post_json("enrichment", payload)
    return pd.DataFrame(rows)


def get_network_media(
    string_ids: list[str],
    species: int,
    required_score: int,
    network_type: str,
    network_flavor: str = "evidence",
) -> dict:
    """Retrieve STRING's own high-resolution PNG, SVG and stable network link.

    Media failures are recorded independently so a temporary image-rendering
    problem does not discard otherwise valid network/enrichment results.
    """
    image_payload = {
        "identifiers": _identifiers(string_ids),
        "species": species,
        "required_score": int(required_score),
        "network_type": network_type,
        "network_flavor": network_flavor,
        "add_color_nodes": 0,
        "add_white_nodes": 0,
        "hide_node_labels": 0,
        "hide_disconnected_nodes": 0,
        "block_structure_pics_in_bubbles": 0,
        "flat_node_design": 0,
        "center_node_labels": 0,
        "custom_label_font_size": 12,
        "caller_identity": CALLER_IDENTITY,
    }

    link_payload = {
        "identifiers": _identifiers(string_ids),
        "species": species,
        "required_score": int(required_score),
        "network_type": network_type,
        "network_flavor": network_flavor,
        "add_color_nodes": 0,
        "add_white_nodes": 0,
        "hide_node_labels": 0,
        "hide_disconnected_nodes": 0,
        "caller_identity": CALLER_IDENTITY,
    }

    media = {
        "highres_png": None,
        "svg": None,
        "link": None,
        "errors": [],
    }

    try:
        media["highres_png"] = _post_binary("highres_image", "network", image_payload)
    except StringAPIError as exc:
        media["errors"].append(f"High-resolution PNG: {exc}")

    time.sleep(1.0)

    try:
        svg_text = _post_text("svg", "network", image_payload)
        media["svg"] = svg_text.encode("utf-8")
    except StringAPIError as exc:
        media["errors"].append(f"SVG: {exc}")

    time.sleep(1.0)

    try:
        link = _post_text("tsv-no-header", "get_link", link_payload).strip()
        if link:
            # The endpoint normally returns only the stable URL. If tabular
            # output ever includes extra columns, keep the URL-like field.
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
    network_flavor: str = "evidence",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """Run the complete STRING workflow with courteous sequential API calls."""
    mapping = map_identifiers(targets, species)
    if mapping.empty:
        raise StringAPIError("None of the submitted identifiers could be mapped by STRING.")

    ids = mapping["stringId"].dropna().astype(str).drop_duplicates().tolist()

    time.sleep(1.0)
    network = get_network(ids, species, required_score, network_type)

    time.sleep(1.0)
    enrichment = get_enrichment(ids, species) if len(ids) >= 2 else pd.DataFrame()

    time.sleep(1.0)
    native_media = get_network_media(
        ids,
        species,
        required_score,
        network_type,
        network_flavor=network_flavor,
    )

    return mapping, network, enrichment, native_media
