from __future__ import annotations

import base64
import gzip
import hashlib
import io
import json
import os
import time
import zipfile
from typing import Callable

import pandas as pd
import requests

from modules.publication_settings import normalize_publication_analysis_settings

DEFAULT_REPOSITORY = "peterDataScientia/Network_Pharmacology"
DEFAULT_WORKFLOW = "publication-enrichment-job.yml"
DEFAULT_REF = "main"
MAX_DISPATCH_PAYLOAD_CHARS = 55000

TABLE_FILES = {
    "mapping_foreground": "mapping_foreground.csv",
    "go_bp_all": "go_bp_all_tested.csv",
    "go_cc_all": "go_cc_all_tested.csv",
    "go_mf_all": "go_mf_all_tested.csv",
    "reactome_all": "reactome_all_tested.csv",
    "unmapped_foreground": "unmapped_foreground.csv",
    "go_bp_raw": "go_bp_raw_significant.csv",
    "go_bp_reduced": "go_bp_reduced.csv",
    "go_cc_raw": "go_cc_raw_significant.csv",
    "go_cc_reduced": "go_cc_reduced.csv",
    "go_mf_raw": "go_mf_raw_significant.csv",
    "go_mf_reduced": "go_mf_reduced.csv",
    "reactome_raw": "reactome_raw_significant.csv",
    "reactome_reduced": "reactome_reduced.csv",
    "background_universe": "background_universe.csv",
    "mapping_background": "mapping_background.csv",
}


class GitHubPublicationError(RuntimeError):
    pass


def _config_value(name: str, default: str = "") -> str:
    value = os.environ.get(name)
    if value is not None and str(value).strip():
        return str(value).strip()

    # Streamlit root-level secrets are normally exposed as environment
    # variables, but this fallback makes the integration robust to deployment
    # configuration differences.
    try:
        import streamlit as st

        value = st.secrets.get(name, default)
        if value is not None and str(value).strip():
            return str(value).strip()
    except Exception:
        pass

    return default


def github_token() -> str:
    return _config_value("PUBLICATION_GITHUB_TOKEN")


def github_repository() -> str:
    return _config_value("PUBLICATION_GITHUB_REPOSITORY", DEFAULT_REPOSITORY)


def github_workflow() -> str:
    return _config_value("PUBLICATION_GITHUB_WORKFLOW", DEFAULT_WORKFLOW)


def github_ref() -> str:
    return _config_value("PUBLICATION_GITHUB_REF", DEFAULT_REF)


def github_actions_configured() -> bool:
    return bool(github_token())


def _headers() -> dict[str, str]:
    token = github_token()
    headers = {
        "Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
        "User-Agent": "Network-Pharmacology-Publication-Enrichment",
    }
    if token:
        headers["Authorization"] = f"Bearer {token}"
    return headers


def _api_url(path: str) -> str:
    return f"https://api.github.com/repos/{github_repository()}{path}"


def _json_request(
    method: str,
    path: str,
    *,
    payload: dict | None = None,
    timeout: int = 45,
) -> dict | list | None:
    try:
        response = requests.request(
            method,
            _api_url(path),
            headers=_headers(),
            json=payload,
            timeout=timeout,
        )
    except requests.RequestException as exc:
        raise GitHubPublicationError(f"Could not reach GitHub Actions: {exc}") from exc

    if response.status_code == 204:
        return None

    if not response.ok:
        detail = response.text.strip().replace("\n", " ")[:1000]
        raise GitHubPublicationError(
            f"GitHub API returned HTTP {response.status_code}: {detail}"
        )

    try:
        return response.json()
    except ValueError as exc:
        raise GitHubPublicationError("GitHub returned an unreadable JSON response.") from exc


def github_actions_status() -> tuple[bool, list[str]]:
    if not github_actions_configured():
        return False, ["PUBLICATION_GITHUB_TOKEN"]

    try:
        workflow = _json_request(
            "GET",
            f"/actions/workflows/{github_workflow()}",
            timeout=30,
        )
    except GitHubPublicationError as exc:
        return False, [str(exc)]

    if not isinstance(workflow, dict) or workflow.get("state") != "active":
        return False, ["Publication enrichment workflow is not active."]

    return True, []


def _emit(callback: Callable[[str], None] | None, message: str) -> None:
    if callback is None:
        return
    try:
        callback(message)
    except Exception:
        pass


def _current_ref_sha() -> str:
    data = _json_request("GET", f"/commits/{github_ref()}", timeout=30)
    if not isinstance(data, dict) or not data.get("sha"):
        raise GitHubPublicationError("Could not determine the current GitHub commit.")
    return str(data["sha"])


def _encode_payload(
    targets: list[str],
    taxon_id: int,
    background_mode: str,
    custom_background: list[str] | None,
    analysis_settings: dict | None,
    ref_sha: str,
) -> tuple[str, str]:
    normalized_settings = normalize_publication_analysis_settings(analysis_settings)
    core = {
        "targets": targets,
        "taxon_id": int(taxon_id),
        "background_mode": background_mode,
        "custom_background": custom_background if background_mode == "custom" else None,
        "analysis_settings": normalized_settings,
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
    encoded = base64.urlsafe_b64encode(gzip.compress(packed, compresslevel=9)).decode("ascii")

    if len(encoded) > MAX_DISPATCH_PAYLOAD_CHARS:
        raise GitHubPublicationError(
            "This request is too large for a GitHub Actions dispatch payload. "
            "Use a smaller custom background or the annotated/package-default background."
        )

    return request_id, encoded


def _matching_runs(request_id: str, ref_sha: str) -> list[dict]:
    data = _json_request(
        "GET",
        (
            f"/actions/workflows/{github_workflow()}/runs"
            f"?event=workflow_dispatch&branch={github_ref()}&per_page=100"
        ),
        timeout=45,
    )
    if not isinstance(data, dict):
        return []

    expected_title = f"Publication enrichment - {request_id}"
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
        matches = _matching_runs(request_id, ref_sha)
        for run in matches:
            if int(run.get("id", 0)) not in existing_ids:
                return run
        _emit(callback, "GitHub accepted the request; waiting for a runner to be assigned…")
        time.sleep(4)
    raise GitHubPublicationError(
        "GitHub accepted the workflow dispatch, but the new run did not appear within 2 minutes."
    )


def _wait_for_completion(
    run: dict,
    callback: Callable[[str], None] | None,
    timeout_seconds: int,
) -> dict:
    run_id = int(run["id"])
    deadline = time.monotonic() + timeout_seconds
    last_status = None

    while time.monotonic() < deadline:
        current = _json_request("GET", f"/actions/runs/{run_id}", timeout=30)
        if not isinstance(current, dict):
            raise GitHubPublicationError("GitHub returned an invalid workflow-run response.")

        status = str(current.get("status", "unknown"))
        if status != last_status:
            if status == "queued":
                _emit(callback, "GitHub job is queued…")
            elif status == "in_progress":
                _emit(callback, "GitHub runner started; R/Bioconductor analysis is running…")
            else:
                _emit(callback, f"GitHub workflow status: {status}")
            last_status = status

        if status == "completed":
            conclusion = str(current.get("conclusion", ""))
            if conclusion != "success":
                url = current.get("html_url", "")
                raise GitHubPublicationError(
                    f"GitHub publication workflow finished with {conclusion!r}. "
                    f"Workflow run: {url}"
                )
            _emit(callback, "GitHub analysis completed; retrieving result artifact…")
            return current

        time.sleep(5)

    raise GitHubPublicationError(
        f"GitHub publication workflow exceeded {timeout_seconds // 60} minutes."
    )


def _download_artifact(
    run_id: int,
    request_id: str,
    callback: Callable[[str], None] | None,
) -> bytes:
    artifact_name = f"publication-enrichment-{request_id}"
    deadline = time.monotonic() + 90

    while time.monotonic() < deadline:
        data = _json_request(
            "GET",
            f"/actions/runs/{run_id}/artifacts?per_page=100",
            timeout=30,
        )
        if isinstance(data, dict):
            for artifact in data.get("artifacts", []):
                if (
                    artifact.get("name") == artifact_name
                    and not artifact.get("expired", False)
                ):
                    artifact_id = int(artifact["id"])
                    try:
                        response = requests.get(
                            _api_url(f"/actions/artifacts/{artifact_id}/zip"),
                            headers=_headers(),
                            timeout=90,
                            allow_redirects=True,
                        )
                    except requests.RequestException as exc:
                        raise GitHubPublicationError(
                            f"Could not download GitHub result artifact: {exc}"
                        ) from exc

                    if not response.ok:
                        raise GitHubPublicationError(
                            "GitHub artifact download failed with "
                            f"HTTP {response.status_code}."
                        )
                    _emit(callback, "Publication result artifact downloaded.")
                    return response.content

        time.sleep(3)

    raise GitHubPublicationError(
        "The workflow completed, but its publication result artifact was not available."
    )


def _member_bytes(zf: zipfile.ZipFile, filename: str) -> bytes | None:
    direct = [name for name in zf.namelist() if name == filename]
    candidates = direct or [
        name for name in zf.namelist() if name.endswith("/" + filename)
    ]
    if not candidates:
        return None
    return zf.read(candidates[0])


def _parse_artifact(data: bytes) -> tuple[dict, dict[str, pd.DataFrame]]:
    try:
        zf = zipfile.ZipFile(io.BytesIO(data))
    except zipfile.BadZipFile as exc:
        raise GitHubPublicationError("GitHub returned a damaged result artifact.") from exc

    summary_bytes = _member_bytes(zf, "summary.json")
    if summary_bytes is None:
        raise GitHubPublicationError("Result artifact does not contain summary.json.")

    try:
        summary = json.loads(summary_bytes.decode("utf-8"))
    except Exception as exc:
        raise GitHubPublicationError("Result summary.json could not be read.") from exc

    tables: dict[str, pd.DataFrame] = {}
    for key, filename in TABLE_FILES.items():
        payload = _member_bytes(zf, filename)
        if payload is None:
            continue
        try:
            tables[key] = pd.read_csv(io.BytesIO(payload))
        except pd.errors.EmptyDataError:
            tables[key] = pd.DataFrame()

    return summary, tables


def run_github_publication_enrichment(
    targets: list[str],
    taxon_id: int,
    background_mode: str,
    custom_background: list[str] | None = None,
    analysis_settings: dict | None = None,
    *,
    status_callback: Callable[[str], None] | None = None,
    timeout_seconds: int = 1800,
) -> tuple[dict, dict[str, pd.DataFrame], dict]:
    if not github_actions_configured():
        raise GitHubPublicationError("PUBLICATION_GITHUB_TOKEN is not configured.")

    _emit(status_callback, "Checking GitHub Actions publication runner…")
    ready, missing = github_actions_status()
    if not ready:
        raise GitHubPublicationError("; ".join(missing))

    ref_sha = _current_ref_sha()
    request_id, payload_b64 = _encode_payload(
        targets,
        taxon_id,
        background_mode,
        custom_background,
        analysis_settings,
        ref_sha,
    )

    matches = _matching_runs(request_id, ref_sha)
    for run in matches:
        if run.get("status") == "completed" and run.get("conclusion") == "success":
            try:
                _emit(
                    status_callback,
                    "Reusing an identical successful GitHub analysis from this code revision…",
                )
                artifact = _download_artifact(
                    int(run["id"]),
                    request_id,
                    status_callback,
                )
                summary, tables = _parse_artifact(artifact)
                return summary, tables, {
                    "executor": "github-actions",
                    "request_id": request_id,
                    "run_id": int(run["id"]),
                    "run_url": run.get("html_url"),
                    "reused": True,
                    "ref_sha": ref_sha,
                }
            except GitHubPublicationError:
                # The previous artifact may have expired; dispatch a fresh run.
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
        _emit(status_callback, "Submitting publication enrichment to GitHub Actions…")
        _json_request(
            "POST",
            f"/actions/workflows/{github_workflow()}/dispatches",
            payload={
                "ref": github_ref(),
                "inputs": {
                    "request_id": request_id,
                    "payload_b64": payload_b64,
                },
            },
            timeout=45,
        )
        active = _wait_for_new_run(
            request_id,
            ref_sha,
            existing_ids,
            status_callback,
        )
    else:
        _emit(status_callback, "An identical GitHub analysis is already running; joining it…")

    completed = _wait_for_completion(
        active,
        status_callback,
        timeout_seconds,
    )
    run_id = int(completed["id"])
    artifact = _download_artifact(run_id, request_id, status_callback)
    summary, tables = _parse_artifact(artifact)

    return summary, tables, {
        "executor": "github-actions",
        "request_id": request_id,
        "run_id": run_id,
        "run_url": completed.get("html_url"),
        "reused": False,
        "ref_sha": ref_sha,
    }
