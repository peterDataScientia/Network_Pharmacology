from __future__ import annotations

import time
from typing import Iterable

import pandas as pd
import requests

STRING_API_BASE = "https://string-db.org/api"
CALLER_IDENTITY = "Network_Pharmacology_Streamlit_App"


class StringAPIError(RuntimeError):
    pass


def _post_json(method: str, payload: dict, timeout: int = 60):
    url = f"{STRING_API_BASE}/json/{method}"
    try:
        response = requests.post(url, data=payload, timeout=timeout)
    except requests.RequestException as exc:
        raise StringAPIError(f"Could not contact STRING: {exc}") from exc
    if response.status_code != 200:
        msg = response.text.strip().replace("\n", " ")[:300]
        raise StringAPIError(f"STRING returned HTTP {response.status_code}: {msg}")
    try:
        return response.json()
    except ValueError as exc:
        raise StringAPIError("STRING returned an unreadable response.") from exc


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


def get_network(string_ids: list[str], species: int, required_score: int, network_type: str) -> pd.DataFrame:
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


def run_string_workflow(
    targets: list[str], species: int, required_score: int, network_type: str
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run three courteous sequential STRING calls.

    STRING asks API users not to flood the service; a one-second pause between
    calls is intentionally retained here.
    """
    mapping = map_identifiers(targets, species)
    if mapping.empty:
        raise StringAPIError("None of the submitted identifiers could be mapped by STRING.")

    ids = mapping["stringId"].dropna().astype(str).drop_duplicates().tolist()
    time.sleep(1.0)
    network = get_network(ids, species, required_score, network_type)
    time.sleep(1.0)
    enrichment = get_enrichment(ids, species) if len(ids) >= 2 else pd.DataFrame()
    return mapping, network, enrichment
