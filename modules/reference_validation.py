from __future__ import annotations

EGCG_RISI_REFERENCE_TARGETS = frozenset({
    "AKT1", "BAD", "BAX", "BCL2L1", "CASP3", "CDKN1A", "CDKN2A", "CTNNB1",
    "DNMT1", "EGF", "EGFR", "GSK3B", "HIF1A", "HMOX1", "IL1B", "IL6",
    "MAPK1", "MAPK3", "MAPK8", "MLH1", "MMP1", "MMP2", "MTOR", "NFE2L2",
    "PARP1", "PTGS2", "RPSA", "SOD2", "STAT3", "TERT", "TNF", "TP53",
})

EGCG_RISI_NETWORK_EXPECTED = {
    "string_version": "12.0",
    "network_type": "functional",
    "required_score": 900,
    "top_n": 10,
    "mapped": 32,
    "connected": 30,
    "edges": 90,
    "hubs": frozenset({"AKT1", "BCL2L1", "CASP3", "STAT3", "TP53"}),
}

EGCG_RISI_TOP10_EXPECTED = {
    "Degree": (
        "TP53", "AKT1", "STAT3", "CASP3", "CTNNB1",
        "HIF1A", "IL6", "MAPK1", "MAPK3", "BCL2L1",
    ),
    "Betweenness": (
        "TP53", "STAT3", "AKT1", "MMP2", "NFE2L2",
        "CASP3", "BCL2L1", "EGFR", "TNF", "IL6",
    ),
    "Closeness": (
        "TP53", "AKT1", "STAT3", "CASP3", "HIF1A",
        "MAPK1", "BCL2L1", "MAPK3", "CDKN1A", "CTNNB1",
    ),
    "Eigenvector": (
        "TP53", "AKT1", "STAT3", "CASP3", "MAPK1",
        "HIF1A", "MAPK3", "CTNNB1", "BCL2L1", "EGFR",
    ),
}


PUBLICATION_VERSION_EXPECTED = {
    "R": "4.6.1",
    "clusterProfiler": "4.20.0",
    "ReactomePA": "1.56.0",
    "AnnotationDbi": "1.74.0",
    "GOSemSim": "2.38.3",
    "organism_db": {
        9606: "3.23.1",
        10090: "3.23.0",
        10116: "3.23.0",
    },
}

EGCG_RISI_ENRICHMENT_DEFAULT_EXPECTED = {
    "go_bp_raw": 1432,
    "go_bp_reduced": 49,
    "go_cc_raw": 23,
    "go_cc_reduced": 17,
    "go_mf_raw": 54,
    "go_mf_reduced": 26,
    "reactome_raw": 276,
    "reactome_reduced": 91,
}


def is_egcg_risi_reference(values) -> bool:
    genes = {
        str(value).strip().upper()
        for value in values
        if str(value).strip()
    }
    return genes == EGCG_RISI_REFERENCE_TARGETS
