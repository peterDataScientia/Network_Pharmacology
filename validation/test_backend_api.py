from __future__ import annotations

import os
import sys
import time

import requests

TARGETS = [
    "AKT1","BAD","BAX","BCL2L1","CASP3","CDKN1A","CDKN2A","CTNNB1",
    "DNMT1","EGF","EGFR","GSK3B","HIF1A","HMOX1","IL1B","IL6","MAPK1",
    "MAPK3","MAPK8","MLH1","MMP1","MMP2","MTOR","NFE2L2","PARP1","PTGS2",
    "RPSA","SOD2","STAT3","TERT","TNF","TP53",
]

EXPECTED = {
    "go_bp_raw": 1432,
    "go_bp_reduced": 49,
    "go_cc_raw": 23,
    "go_cc_reduced": 17,
    "go_mf_raw": 54,
    "go_mf_reduced": 26,
    "reactome_raw": 276,
    "reactome_reduced": 91,
}


def main() -> int:
    base = os.environ.get("BACKEND_TEST_URL", "http://127.0.0.1:10000").rstrip("/")
    token = os.environ.get("BACKEND_API_TOKEN", "test-token")
    headers = {"Authorization": f"Bearer {token}"}

    deadline = time.time() + 180
    health = None
    while time.time() < deadline:
        try:
            r = requests.get(f"{base}/health", timeout=10)
            if r.ok:
                health = r.json()
                if health.get("status") == "ok":
                    break
        except requests.RequestException:
            pass
        time.sleep(3)

    if not health or health.get("status") != "ok":
        print("Backend did not become ready:", health)
        return 2

    print("HEALTH:", health)
    versions = health.get("versions", {})
    required_versions = {
        "clusterProfiler": "4.20.0",
        "ReactomePA": "1.56.0",
        "org.Hs.eg.db": "3.23.1",
        "org.Mm.eg.db": "3.23.0",
        "org.Rn.eg.db": "3.23.0",
    }
    for package, expected in required_versions.items():
        actual = versions.get(package)
        if actual != expected:
            print(f"VERSION MISMATCH {package}: expected={expected}, actual={actual}")
            return 3

    payload = {
        "targets": TARGETS,
        "organism": "human",
        "background_mode": "default",
        "custom_background": [],
    }
    started = time.time()
    response = requests.post(
        f"{base}/v1/enrichment",
        headers=headers,
        json=payload,
        timeout=1250,
    )
    elapsed = time.time() - started
    print(f"API elapsed: {elapsed:.1f} s")
    if response.status_code != 200:
        print("API failure:", response.status_code, response.text[:5000])
        return 4

    data = response.json()
    summary = data.get("summary", {})
    counts = summary.get("counts", {})
    print("COUNTS:", counts)

    failed = False
    for key, expected in EXPECTED.items():
        actual = counts.get(key)
        marker = "MATCH" if actual == expected else "DIFF"
        print(f"{key:20s} expected={expected:4d} observed={str(actual):>4} {marker}")
        failed = failed or actual != expected

    tables = data.get("tables", {})
    if len(tables.get("go_bp_raw", [])) != EXPECTED["go_bp_raw"]:
        print("Returned GO-BP raw table length does not match summary.")
        failed = True
    if len(tables.get("reactome_reduced", [])) != EXPECTED["reactome_reduced"]:
        print("Returned Reactome reduced table length does not match summary.")
        failed = True

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
