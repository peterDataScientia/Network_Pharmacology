#!/usr/bin/env python3
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

FILES = [
    "mapping_foreground.csv",
    "unmapped_foreground.csv",
    "go_bp_raw_significant.csv",
    "go_bp_reduced.csv",
    "go_cc_raw_significant.csv",
    "go_cc_reduced.csv",
    "go_mf_raw_significant.csv",
    "go_mf_reduced.csv",
    "reactome_raw_significant.csv",
    "reactome_reduced.csv",
]

NUMERIC_COLUMNS = {
    "pvalue", "p.adjust", "qvalue", "Count", "BgRatio", "GeneRatio",
    "RedundancyCluster"
}

def maybe_number(s: str):
    s = s.strip()
    if s == "":
        return None
    try:
        return float(s)
    except ValueError:
        return None

def normalize_value(col: str, value: str):
    value = value.strip()

    # Gene membership order is scientifically irrelevant.
    if col == "geneID" and value:
        return "/".join(sorted(x for x in value.split("/") if x))

    # Ratios such as 5/123 are categorical strings, not decimal values.
    if col in {"GeneRatio", "BgRatio"}:
        return value

    n = maybe_number(value)
    if n is not None:
        return ("NUM", n)

    return value

def load_rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    return rows

def stable_key(row):
    for key in ("ID", "SYMBOL", "ENTREZID"):
        if key in row:
            return str(row.get(key, ""))
    return repr(sorted(row.items()))

def compare_rows(a_rows, b_rows, name):
    if len(a_rows) != len(b_rows):
        raise AssertionError(f"{name}: row count differs: {len(a_rows)} vs {len(b_rows)}")

    a_rows = sorted(a_rows, key=stable_key)
    b_rows = sorted(b_rows, key=stable_key)

    for i, (a, b) in enumerate(zip(a_rows, b_rows), 1):
        if set(a) != set(b):
            raise AssertionError(f"{name}: columns differ at row {i}")

        for col in a:
            av = normalize_value(col, a[col])
            bv = normalize_value(col, b[col])

            if isinstance(av, tuple) and av[0] == "NUM" and isinstance(bv, tuple) and bv[0] == "NUM":
                # Same calculation should be extremely close; tolerate only serialization-level noise.
                if not math.isclose(av[1], bv[1], rel_tol=1e-12, abs_tol=1e-15):
                    raise AssertionError(
                        f"{name}: numeric value differs for {stable_key(a)} column {col}: "
                        f"{a[col]!r} vs {b[col]!r}"
                    )
            elif av != bv:
                raise AssertionError(
                    f"{name}: value differs for {stable_key(a)} column {col}: "
                    f"{a[col]!r} vs {b[col]!r}"
                )

def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: compare_publication_outputs.py <baseline_dir> <candidate_dir>")

    base = Path(sys.argv[1])
    cand = Path(sys.argv[2])

    for name in FILES:
        a = base / name
        b = cand / name
        if not a.exists() or not b.exists():
            raise SystemExit(f"Missing comparison file: {name}")
        compare_rows(load_rows(a), load_rows(b), name)
        print(f"OK {name}")

    print("SEMANTIC_CACHE_SCIENTIFIC_PARITY_OK")

if __name__ == "__main__":
    main()
