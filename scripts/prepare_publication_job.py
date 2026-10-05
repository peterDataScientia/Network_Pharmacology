from __future__ import annotations

import base64
import gzip
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from modules.publication_settings import (
    PublicationSettingsError,
    normalize_publication_analysis_settings,
)

TAXON_TO_ORGANISM = {
    9606: "human",
    10090: "mouse",
    10116: "rat",
}


def _dedupe(values):
    out = []
    seen = set()
    for value in values or []:
        text = str(value).strip()
        if not text:
            continue
        key = text.upper()
        if key in seen:
            continue
        seen.add(key)
        out.append(text)
    return out


def main() -> None:
    encoded = os.environ.get("PAYLOAD_B64", "").strip()
    if not encoded:
        raise SystemExit("PAYLOAD_B64 is empty.")

    try:
        raw = gzip.decompress(base64.urlsafe_b64decode(encoded.encode("ascii")))
        payload = json.loads(raw.decode("utf-8"))
    except Exception as exc:
        raise SystemExit(f"Could not decode job payload: {exc}") from exc

    targets = _dedupe(payload.get("targets"))
    if not 2 <= len(targets) <= 1000:
        raise SystemExit("Target count must be between 2 and 1000.")

    try:
        taxon_id = int(payload.get("taxon_id"))
    except Exception as exc:
        raise SystemExit("taxon_id must be an integer.") from exc

    organism = TAXON_TO_ORGANISM.get(taxon_id)
    if organism is None:
        raise SystemExit("Supported taxon IDs are 9606, 10090 and 10116.")

    background_mode = str(payload.get("background_mode", "")).strip().lower()
    if background_mode not in {"default", "annotated", "custom"}:
        raise SystemExit("background_mode must be default, annotated or custom.")

    custom_background = _dedupe(payload.get("custom_background"))
    if background_mode == "custom" and not custom_background:
        raise SystemExit("Custom background mode requires at least one background gene.")

    try:
        analysis_settings = normalize_publication_analysis_settings(
            payload.get("analysis_settings")
        )
    except PublicationSettingsError as exc:
        raise SystemExit(str(exc)) from exc

    request_id = str(payload.get("request_id", "")).strip()
    if not request_id:
        raise SystemExit("request_id is required.")

    root = Path("job_input")
    root.mkdir(parents=True, exist_ok=True)
    Path("job_output").mkdir(parents=True, exist_ok=True)

    (root / "foreground_symbols.txt").write_text(
        "\n".join(targets) + "\n",
        encoding="utf-8",
    )

    background_arg = "NONE"
    if background_mode == "custom":
        (root / "custom_background_symbols.txt").write_text(
            "\n".join(custom_background) + "\n",
            encoding="utf-8",
        )
        background_arg = "/job_input/custom_background_symbols.txt"

    settings_path = root / "analysis_settings.json"
    settings_path.write_text(
        json.dumps(analysis_settings, indent=2),
        encoding="utf-8",
    )

    normalized = {
        "request_id": request_id,
        "targets": targets,
        "taxon_id": taxon_id,
        "organism": organism,
        "background_mode": background_mode,
        "custom_background": custom_background if background_mode == "custom" else None,
        "analysis_settings": analysis_settings,
    }
    (root / "request.json").write_text(
        json.dumps(normalized, indent=2),
        encoding="utf-8",
    )

    github_output = os.environ.get("GITHUB_OUTPUT")
    if github_output:
        with open(github_output, "a", encoding="utf-8") as fh:
            fh.write(f"organism={organism}\n")
            fh.write(f"background_mode={background_mode}\n")
            fh.write(f"background_arg={background_arg}\n")
            fh.write(f"request_id={request_id}\n")
            fh.write("settings_path=job_input/analysis_settings.json\n")
            fh.write("settings_docker_path=/job_input/analysis_settings.json\n")

    print(
        f"Prepared publication job {request_id}: "
        f"{len(targets)} targets, {organism}, background={background_mode}"
    )


if __name__ == "__main__":
    main()
