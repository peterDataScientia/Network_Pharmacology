from __future__ import annotations

import base64
import gzip
import hashlib
import io
import json
import time
import zipfile
from typing import Callable

import pandas as pd
import requests

from modules.publication_github import (
    GitHubPublicationError,
    _api_url,
    _headers,
    _json_request,
    github_actions_configured,
    github_ref,
)

CENTRALITY_WORKFLOW = "network-centrality-job.yml"
MAX_DISPATCH_PAYLOAD_CHARS = 55000


class CentralityGitHubError(RuntimeError):
    pass


def _emit(callback: Callable[[str], None] | None, message: str) -> None:
    if callback is None:
        return
    try:
        callback(message)
    except Exception:
        pass


def _canonical_edge_pairs(network: pd.DataFrame) -> list[tuple[str, str]]:
    if network.empty:
        return []
    pairs: set[tuple[str, str]] = set()
    for _, row in network.iterrows():
        a = str(row.get("preferredName_A", "")).strip()
        b = str(row.get("preferredName_B", "")).strip()
        if not a or not b or a == b:
            continue
        pairs.add(tuple(sorted((a, b))))
    return sorted(pairs)


def network_edge_signature(network: pd.DataFrame) -> tuple[str, int, int]:
    pairs = _canonical_edge_pairs(network)
    payload = "".join(f"{a}\t{b}\n" for a, b in pairs).encode("utf-8")
    digest = hashlib.sha256(payload).hexdigest()
    nodes = {node for pair in pairs for node in pair}
    return digest, len(pairs), len(nodes)


def centrality_runner_status() -> tuple[bool, list[str]]:
    if not github_actions_configured():
        return False, ["PUBLICATION_GITHUB_TOKEN"]
    try:
        workflow = _json_request(
            "GET",
            f"/actions/workflows/{CENTRALITY_WORKFLOW}",
            timeout=30,
        )
    except Exception as exc:
        return False, [str(exc)]

    if not isinstance(workflow, dict) or workflow.get("state") != "active":
        return False, ["R/igraph centrality workflow is not active."]
    return True, []


def _current_ref_sha() -> str:
    try:
        data = _json_request("GET", f"/commits/{github_ref()}", timeout=30)
    except GitHubPublicationError as exc:
        raise CentralityGitHubError(str(exc)) from None
    if not isinstance(data, dict) or not data.get("sha"):
        raise CentralityGitHubError("Could not determine the current GitHub commit.")
    return str(data["sha"])


def _encode_payload(
    *,
    string_ids: list[str],
    taxon_id: int,
    required_score: int,
    network_type: str,
    network_flavor: str,
    active_sources: list[str],
    add_nodes: int,
    string_version: str,
    top_n: int,
    expected_edge_hash: str,
    expected_edges: int,
    expected_nodes: int,
    ref_sha: str,
) -> tuple[str, str]:
    core = {
        "string_ids": list(dict.fromkeys(str(x).strip() for x in string_ids if str(x).strip())),
        "taxon_id": int(taxon_id),
        "required_score": int(required_score),
        "network_type": str(network_type),
        "network_flavor": str(network_flavor),
        "active_sources": list(dict.fromkeys(str(x) for x in active_sources)),
        "add_nodes": max(0, int(add_nodes)),
        "string_version": str(string_version),
        "top_n": int(top_n),
        "expected_edge_hash": str(expected_edge_hash),
        "expected_edges": int(expected_edges),
        "expected_nodes": int(expected_nodes),
    }
    canonical = json.dumps(
        core,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")

    request_id = hashlib.sha256(
        ref_sha.encode("ascii") + b"\0" + canonical
    ).hexdigest()[:20]

    payload = dict(core)
    payload["request_id"] = request_id
    packed = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    encoded = base64.urlsafe_b64encode(
        gzip.compress(packed, compresslevel=9)
    ).decode("ascii")

    if len(encoded) > MAX_DISPATCH_PAYLOAD_CHARS:
        raise CentralityGitHubError(
            "The centrality request is too large for GitHub Actions dispatch."
        )
    return request_id, encoded


def _matching_runs(request_id: str, ref_sha: str) -> list[dict]:
    try:
        data = _json_request(
            "GET",
            (
                f"/actions/workflows/{CENTRALITY_WORKFLOW}/runs"
                f"?event=workflow_dispatch&branch={github_ref()}&per_page=100"
            ),
            timeout=45,
        )
    except GitHubPublicationError as exc:
        raise CentralityGitHubError(str(exc)) from None
    if not isinstance(data, dict):
        return []

    expected_title = f"R centrality - {request_id}"
    runs = [
        run
        for run in data.get("workflow_runs", [])
        if run.get("display_title") == expected_title
        and run.get("head_sha") == ref_sha
    ]
    return sorted(runs, key=lambda run: int(run.get("id", 0)), reverse=True)


def _wait_for_new_run(
    request_id: str,
    ref_sha: str,
    existing_ids: set[int],
    callback: Callable[[str], None] | None,
) -> dict:
    deadline = time.monotonic() + 120
    while time.monotonic() < deadline:
        for run in _matching_runs(request_id, ref_sha):
            if int(run.get("id", 0)) not in existing_ids:
                return run
        _emit(callback, "GitHub accepted the centrality request; waiting for a runner…")
        time.sleep(4)
    raise CentralityGitHubError(
        "GitHub accepted the centrality workflow, but the run did not appear within 2 minutes."
    )


def _wait_for_completion(
    run: dict,
    callback: Callable[[str], None] | None,
    timeout_seconds: int = 900,
) -> dict:
    run_id = int(run["id"])
    deadline = time.monotonic() + timeout_seconds
    last_status = None

    while time.monotonic() < deadline:
        try:
            current = _json_request("GET", f"/actions/runs/{run_id}", timeout=30)
        except GitHubPublicationError as exc:
            raise CentralityGitHubError(str(exc)) from None
        if not isinstance(current, dict):
            raise CentralityGitHubError("GitHub returned an invalid centrality run response.")

        status = str(current.get("status", "unknown"))
        if status != last_status:
            if status == "queued":
                _emit(callback, "R/igraph centrality job is queued…")
            elif status == "in_progress":
                _emit(callback, "GitHub runner started; R/igraph centrality is running…")
            else:
                _emit(callback, f"Centrality workflow status: {status}")
            last_status = status

        if status == "completed":
            conclusion = str(current.get("conclusion", ""))
            if conclusion != "success":
                raise CentralityGitHubError(
                    f"R/igraph centrality workflow finished with {conclusion!r}. "
                    f"Workflow run: {current.get('html_url', '')}"
                )
            _emit(callback, "R/igraph centrality completed; retrieving results…")
            return current

        time.sleep(5)

    raise CentralityGitHubError("R/igraph centrality workflow exceeded 15 minutes.")


def _download_artifact(
    run_id: int,
    request_id: str,
    callback: Callable[[str], None] | None,
) -> bytes:
    artifact_name = f"network-centrality-{request_id}"
    deadline = time.monotonic() + 90

    while time.monotonic() < deadline:
        try:
            data = _json_request(
                "GET",
                f"/actions/runs/{run_id}/artifacts?per_page=100",
                timeout=30,
            )
        except GitHubPublicationError as exc:
            raise CentralityGitHubError(str(exc)) from None
        if isinstance(data, dict):
            for artifact in data.get("artifacts", []):
                if artifact.get("name") == artifact_name and not artifact.get("expired", False):
                    artifact_id = int(artifact["id"])
                    try:
                        response = requests.get(
                            _api_url(f"/actions/artifacts/{artifact_id}/zip"),
                            headers=_headers(),
                            timeout=90,
                            allow_redirects=True,
                        )
                    except requests.RequestException as exc:
                        raise CentralityGitHubError(
                            f"Could not download R/igraph centrality artifact: {exc}"
                        ) from exc
                    if not response.ok:
                        raise CentralityGitHubError(
                            f"Centrality artifact download failed with HTTP {response.status_code}."
                        )
                    _emit(callback, "R/igraph centrality artifact downloaded.")
                    return response.content
        time.sleep(3)

    raise CentralityGitHubError(
        "The centrality workflow completed, but its result artifact was unavailable."
    )


def _member_bytes(zf: zipfile.ZipFile, filename: str) -> bytes | None:
    direct = [name for name in zf.namelist() if name == filename]
    candidates = direct or [
        name for name in zf.namelist() if name.endswith("/" + filename)
    ]
    if not candidates:
        return None
    return zf.read(candidates[0])


def _parse_artifact(data: bytes) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict, dict]:
    try:
        zf = zipfile.ZipFile(io.BytesIO(data))
    except zipfile.BadZipFile as exc:
        raise CentralityGitHubError("GitHub returned a damaged centrality artifact.") from exc

    def read_csv(filename: str) -> pd.DataFrame:
        payload = _member_bytes(zf, filename)
        if payload is None:
            raise CentralityGitHubError(f"Centrality artifact is missing {filename}.")
        try:
            return pd.read_csv(io.BytesIO(payload))
        except pd.errors.EmptyDataError:
            return pd.DataFrame()

    summary_bytes = _member_bytes(zf, "centrality_summary.json")
    metadata_bytes = _member_bytes(zf, "centrality_job_metadata.json")
    if summary_bytes is None or metadata_bytes is None:
        raise CentralityGitHubError("Centrality artifact is missing provenance metadata.")

    summary = json.loads(summary_bytes.decode("utf-8"))
    metadata = json.loads(metadata_bytes.decode("utf-8"))

    return (
        read_csv("centrality_all_genes.csv"),
        read_csv("centrality_consensus_rankings.csv"),
        read_csv("consensus_hubs_4of4.csv"),
        summary,
        metadata,
    )


def run_r_igraph_centrality(
    *,
    mapping: pd.DataFrame,
    network: pd.DataFrame,
    taxon_id: int,
    required_score: int,
    network_type: str,
    string_version: str,
    top_n: int,
    network_flavor: str = "evidence",
    active_sources: list[str] | None = None,
    add_nodes: int = 0,
    status_callback: Callable[[str], None] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, int, dict]:
    if network_type == "regulatory":
        raise CentralityGitHubError(
            "The validated 4/4 consensus hub workflow is undirected and is not "
            "applied to STRING regulatory networks."
        )

    if network.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), 0, {
            "executor": "github-actions",
            "centrality_engine": "R/igraph",
            "nodes": 0,
            "edges": 0,
        }

    ready, missing = centrality_runner_status()
    if not ready:
        raise CentralityGitHubError("; ".join(missing))

    if "stringId" not in mapping.columns:
        raise CentralityGitHubError("STRING mapping does not contain stringId values.")

    string_ids = (
        mapping["stringId"].dropna().astype(str).drop_duplicates().tolist()
    )
    edge_hash, expected_edges, expected_nodes = network_edge_signature(network)
    ref_sha = _current_ref_sha()

    request_id, payload_b64 = _encode_payload(
        string_ids=string_ids,
        taxon_id=taxon_id,
        required_score=required_score,
        network_type=network_type,
        network_flavor=network_flavor,
        active_sources=active_sources or [],
        add_nodes=add_nodes,
        string_version=string_version,
        top_n=top_n,
        expected_edge_hash=edge_hash,
        expected_edges=expected_edges,
        expected_nodes=expected_nodes,
        ref_sha=ref_sha,
    )

    matches = _matching_runs(request_id, ref_sha)
    for run in matches:
        if run.get("status") == "completed" and run.get("conclusion") == "success":
            try:
                _emit(status_callback, "Reusing identical validated R/igraph centrality results…")
                artifact = _download_artifact(
                    int(run["id"]),
                    request_id,
                    status_callback,
                )
                centrality, ranked, hubs, summary, metadata = _parse_artifact(artifact)
                return (
                    centrality,
                    ranked,
                    hubs,
                    int(summary.get("top_n_effective", 0)),
                    {
                        **metadata,
                        "request_id": request_id,
                        "run_id": int(run["id"]),
                        "run_url": run.get("html_url"),
                        "reused": True,
                        "ref_sha": ref_sha,
                        "summary": summary,
                    },
                )
            except CentralityGitHubError:
                break

    active = next(
        (
            run
            for run in matches
            if run.get("status") in {"queued", "in_progress"}
        ),
        None,
    )

    if active is None:
        existing_ids = {int(run.get("id", 0)) for run in matches}
        _emit(status_callback, "Submitting network centrality to R/igraph on GitHub Actions…")
        try:
            _json_request(
                "POST",
                f"/actions/workflows/{CENTRALITY_WORKFLOW}/dispatches",
                payload={
                    "ref": github_ref(),
                    "inputs": {
                        "request_id": request_id,
                        "payload_b64": payload_b64,
                    },
                },
                timeout=45,
            )
        except GitHubPublicationError as exc:
            raise CentralityGitHubError(str(exc)) from None
        active = _wait_for_new_run(
            request_id,
            ref_sha,
            existing_ids,
            status_callback,
        )
    else:
        _emit(status_callback, "An identical R/igraph centrality job is already running; joining it…")

    completed = _wait_for_completion(active, status_callback)
    run_id = int(completed["id"])
    artifact = _download_artifact(run_id, request_id, status_callback)
    centrality, ranked, hubs, summary, metadata = _parse_artifact(artifact)

    return (
        centrality,
        ranked,
        hubs,
        int(summary.get("top_n_effective", 0)),
        {
            **metadata,
            "request_id": request_id,
            "run_id": run_id,
            "run_url": completed.get("html_url"),
            "reused": False,
            "ref_sha": ref_sha,
            "summary": summary,
        },
    )
