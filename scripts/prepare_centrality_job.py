from __future__ import annotations

import base64
import csv
import gzip
import hashlib
import io
import json
import os
import sys
import urllib.parse
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from modules.string_filters import applicable_sources, filter_and_normalize_rows

STRING_API_BASES = {
    "12.0": "https://version-12.string-db.org/api",
    "12.5": "https://version-12-5.string-db.org/api",
}
CALLER_IDENTITY = "Network_Pharmacology_Streamlit_App"


def _canonical_edges(rows: list[dict]) -> list[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for row in rows:
        a = str(row.get("preferredName_A", "")).strip()
        b = str(row.get("preferredName_B", "")).strip()
        if not a or not b or a == b:
            continue
        pairs.add(tuple(sorted((a, b))))
    return sorted(pairs)


def _edge_hash(pairs: list[tuple[str, str]]) -> str:
    payload = "".join(f"{a}\t{b}\n" for a, b in pairs).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _decode_payload() -> dict:
    encoded = os.environ.get("PAYLOAD_B64", "").strip()
    if not encoded:
        raise SystemExit("PAYLOAD_B64 is empty.")
    try:
        raw = gzip.decompress(base64.urlsafe_b64decode(encoded.encode("ascii")))
        return json.loads(raw.decode("utf-8"))
    except Exception as exc:
        raise SystemExit(f"Could not decode centrality payload: {exc}") from exc


def _fetch_network(payload: dict) -> list[dict]:
    version = str(payload["string_version"])
    base = STRING_API_BASES.get(version)
    if base is None:
        raise SystemExit(f"Unsupported STRING version: {version}")

    identifiers = [
        str(x).strip()
        for x in payload.get("string_ids", [])
        if str(x).strip()
    ]
    if not identifiers:
        raise SystemExit("No STRING identifiers were supplied.")

    network_type = str(payload["network_type"])
    network_flavor = str(payload.get("network_flavor", "evidence"))
    form_values = {
        "identifiers": "\r".join(identifiers),
        "species": int(payload["taxon_id"]),
        "required_score": int(payload["required_score"]),
        "network_type": network_type,
        "add_nodes": max(0, int(payload.get("add_nodes", 0))),
        "caller_identity": CALLER_IDENTITY,
    }
    if network_type == "functional" and network_flavor == "typed":
        form_values["network_flavor"] = "typed"

    form = urllib.parse.urlencode(form_values).encode("utf-8")

    request = urllib.request.Request(
        f"{base}/json/network",
        data=form,
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            rows = json.load(response)
    except Exception as exc:
        raise SystemExit(f"STRING network request failed: {exc}") from exc

    if not isinstance(rows, list):
        raise SystemExit("STRING returned an invalid network response.")

    active_sources = payload.get("active_sources") or list(
        applicable_sources(network_type)
    )
    return filter_and_normalize_rows(
        rows,
        active_sources=active_sources,
        required_score=int(payload["required_score"]),
        network_type=network_type,
        network_flavor=network_flavor,
    )


def main() -> None:
    payload = _decode_payload()

    request_id = str(payload.get("request_id", "")).strip()
    if not request_id:
        raise SystemExit("request_id is required.")

    top_n = int(payload.get("top_n", 0))
    if top_n < 1:
        raise SystemExit("top_n must be positive.")

    rows = _fetch_network(payload)
    pairs = _canonical_edges(rows)
    observed_hash = _edge_hash(pairs)

    expected_hash = str(payload.get("expected_edge_hash", "")).strip()
    expected_edges = int(payload.get("expected_edges", -1))
    expected_nodes = int(payload.get("expected_nodes", -1))

    if expected_hash and observed_hash != expected_hash:
        raise SystemExit(
            "STRING edge-set mismatch between Streamlit and GitHub runner: "
            f"observed {observed_hash}, expected {expected_hash}"
        )
    if expected_edges >= 0 and len(pairs) != expected_edges:
        raise SystemExit(
            f"STRING edge-count mismatch: observed {len(pairs)}, expected {expected_edges}"
        )

    nodes = sorted({node for pair in pairs for node in pair})
    if expected_nodes >= 0 and len(nodes) != expected_nodes:
        raise SystemExit(
            f"STRING connected-node mismatch: observed {len(nodes)}, expected {expected_nodes}"
        )

    input_dir = Path("centrality_input")
    output_dir = Path("centrality_output")
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    normalized = dict(payload)
    normalized["observed_edge_hash"] = observed_hash
    normalized["observed_edges"] = len(pairs)
    normalized["observed_nodes"] = len(nodes)
    (input_dir / "request.json").write_text(
        json.dumps(normalized, indent=2),
        encoding="utf-8",
    )

    with (input_dir / "network.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["preferredName_A", "preferredName_B"])
        writer.writerows(pairs)

    github_output = os.environ.get("GITHUB_OUTPUT")
    if github_output:
        with open(github_output, "a", encoding="utf-8") as fh:
            fh.write(f"request_id={request_id}\n")
            fh.write(f"top_n={top_n}\n")

    print(
        f"Prepared centrality job {request_id}: "
        f"{len(nodes)} connected nodes, {len(pairs)} edges, top_n={top_n}"
    )


if __name__ == "__main__":
    main()
