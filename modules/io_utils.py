from __future__ import annotations

import io
import re
import zipfile
from datetime import datetime, timezone

import pandas as pd


def normalize_targets(text: str) -> list[str]:
    tokens = re.split(r"[\s,;|]+", text or "")
    out, seen = [], set()
    for token in tokens:
        token = token.strip()
        if token and token.upper() not in seen:
            seen.add(token.upper())
            out.append(token)
    return out


def targets_from_upload(uploaded_file) -> list[str]:
    if uploaded_file is None:
        return []
    raw = uploaded_file.getvalue()
    name = uploaded_file.name.lower()

    if name.endswith((".csv", ".tsv")):
        sep = "\t" if name.endswith(".tsv") else ","
        try:
            df = pd.read_csv(io.BytesIO(raw), sep=sep)
            preferred = [c for c in df.columns if str(c).lower() in {"gene", "genes", "symbol", "target", "targets", "gene_symbol"}]
            col = preferred[0] if preferred else df.columns[0]
            return normalize_targets("\n".join(df[col].dropna().astype(str)))
        except Exception:
            pass

    return normalize_targets(raw.decode("utf-8", errors="ignore"))


def dataframe_tsv(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False, sep="\t").encode("utf-8")


def build_results_zip(files: dict[str, bytes], metadata: str) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, content in files.items():
            zf.writestr(name, content)
        zf.writestr("METHODS_AND_SETTINGS.txt", metadata)
        zf.writestr(
            "generated_utc.txt",
            datetime.now(timezone.utc).isoformat(),
        )
    buf.seek(0)
    return buf.getvalue()
