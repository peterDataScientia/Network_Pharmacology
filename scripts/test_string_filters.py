from __future__ import annotations

from modules.string_filters import combine_channel_scores, filter_and_normalize_rows


def main() -> None:
    combined = combine_channel_scores([0.621, 0.585])
    assert 0.78 < combined < 0.82, combined

    standard = [{
        "stringId_A": "9606.A",
        "stringId_B": "9606.B",
        "preferredName_A": "A",
        "preferredName_B": "B",
        "score": 0.92,
        "nscore": 0.0,
        "fscore": 0.0,
        "pscore": 0.0,
        "ascore": 0.0,
        "escore": 0.30,
        "dscore": 0.0,
        "tscore": 0.0,
    }]
    all_sources = filter_and_normalize_rows(
        standard,
        active_sources=[
            "textmining",
            "experiments",
            "databases",
            "coexpression",
            "neighborhood",
            "fusion",
            "cooccurrence",
        ],
        required_score=900,
        network_type="functional",
        network_flavor="evidence",
    )
    assert len(all_sources) == 1
    assert abs(all_sources[0]["analysis_score"] - 0.92) < 1e-12

    experimental_only = filter_and_normalize_rows(
        standard,
        active_sources=["experiments"],
        required_score=400,
        network_type="functional",
        network_flavor="evidence",
    )
    assert experimental_only == []

    regulatory = [{
        "source_string_id": "9606.A",
        "target_string_id": "9606.B",
        "source_preferred_name": "A",
        "target_preferred_name": "B",
        "experimental_score": 0.80,
        "database_score": 0.0,
        "textmining_score": 0.0,
        "combined_score": 0.90,
        "sign": "pos",
    }]
    filtered_reg = filter_and_normalize_rows(
        regulatory,
        active_sources=["experiments"],
        required_score=700,
        network_type="regulatory",
        network_flavor="confidence",
    )
    assert len(filtered_reg) == 1
    assert filtered_reg[0]["preferredName_A"] == "A"
    assert filtered_reg[0]["preferredName_B"] == "B"
    assert abs(filtered_reg[0]["analysis_score"] - 0.80) < 1e-12

    typed = [{
        "stringId_A": "9606.A",
        "stringId_B": "9606.B",
        "preferredName_A": "A",
        "preferredName_B": "B",
        "functional_combined_score": 0.95,
        "functional_experimental_score": 0.75,
        "functional_database_score": 0.0,
        "functional_textmining_score": 0.0,
        "functional_coexpression_score": 0.0,
        "functional_neighborhood_on_chromosome_score": 0.0,
        "functional_gene_fusion_score": 0.0,
        "functional_phylogenetic_cooccurrence_score": 0.0,
    }]
    filtered_typed = filter_and_normalize_rows(
        typed,
        active_sources=["experiments"],
        required_score=700,
        network_type="functional",
        network_flavor="typed",
    )
    assert len(filtered_typed) == 1
    assert abs(filtered_typed[0]["analysis_score"] - 0.75) < 1e-12

    print("STRING_FILTER_TESTS_OK")


if __name__ == "__main__":
    main()
