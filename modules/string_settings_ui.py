from __future__ import annotations

from modules.string_filters import (
    FUNCTIONAL_SOURCES,
    PHYSICAL_SOURCES,
    REGULATORY_SOURCES,
    SOURCE_LABELS,
)


SCORE_PRESETS = {
    "Highest confidence (0.900)": 900,
    "High confidence (0.700)": 700,
    "Medium confidence (0.400)": 400,
    "Low confidence (0.150)": 150,
    "Custom value": None,
}

SHELL_PRESETS = {
    "Query proteins only / none": 0,
    "No more than 5 interactors": 5,
    "No more than 10 interactors": 10,
    "No more than 20 interactors": 20,
    "No more than 50 interactors": 50,
    "Custom value": None,
}


def _shell_control(st, label: str, key: str) -> int:
    option = st.selectbox(label, list(SHELL_PRESETS), key=f"{key}_preset")
    value = SHELL_PRESETS[option]
    if value is None:
        value = int(
            st.number_input(
                f"{label} custom maximum",
                min_value=0,
                max_value=500,
                value=25,
                step=1,
                key=f"{key}_custom",
            )
        )
    return int(value)


def render_string_settings(st) -> dict:
    st.header("STRING settings")

    species_label = st.selectbox(
        "Organism",
        [
            "Homo sapiens (Human)",
            "Mus musculus (Mouse)",
            "Rattus norvegicus (Rat)",
        ],
        index=0,
    )
    species = {
        "Homo sapiens (Human)": 9606,
        "Mus musculus (Mouse)": 10090,
        "Rattus norvegicus (Rat)": 10116,
    }[species_label]

    string_version_label = st.selectbox(
        "STRING database version",
        [
            "v12.0 — manuscript/reproducibility",
            "v12.5 — newer/current pinned release",
        ],
        index=0,
        help=(
            "Network topology can change between STRING releases. Keep the version "
            "pinned for reproducible analyses."
        ),
    )
    string_version = "12.0" if string_version_label.startswith("v12.0") else "12.5"

    with st.expander("Basic Settings", expanded=True):
        network_type_label = st.radio(
            "Network type",
            [
                "Full STRING network",
                "Physical subnetwork",
                "Regulatory subnetwork",
            ],
            index=0,
            help="Choose functional associations, physical relationships, or directed regulation.",
        )
        network_type = {
            "Full STRING network": "functional",
            "Physical subnetwork": "physical",
            "Regulatory subnetwork": "regulatory",
        }[network_type_label]

        flavor_options = ["Evidence", "Confidence"]
        if network_type == "functional":
            flavor_options.append("Typed")
        network_flavor_label = st.radio(
            "Meaning of network edges",
            flavor_options,
            index=0,
            help=(
                "Evidence uses source-specific edge styling; confidence emphasizes data support. "
                "Typed overlays physical and regulatory relationships on a functional network."
            ),
        )
        network_flavor = network_flavor_label.lower()

        if network_type == "functional":
            allowed_sources = FUNCTIONAL_SOURCES
        elif network_type == "physical":
            allowed_sources = PHYSICAL_SOURCES
        else:
            allowed_sources = REGULATORY_SOURCES

        st.markdown("**Active interaction sources**")
        active_sources = []
        for source in allowed_sources:
            if st.checkbox(
                SOURCE_LABELS[source],
                value=True,
                key=f"string_source_{source}_{network_type}",
            ):
                active_sources.append(source)

        unavailable = [
            key for key in FUNCTIONAL_SOURCES if key not in allowed_sources
        ]
        if unavailable:
            st.caption(
                "Not applicable to this network type: "
                + ", ".join(SOURCE_LABELS[key] for key in unavailable)
                + "."
            )

        st.checkbox(
            "Evidence transfer",
            value=True,
            disabled=True,
            help=(
                "STRING's public network API does not expose direct and transferred "
                "channel scores separately. This remains enabled so the app does not "
                "pretend to provide a filter it cannot reproduce exactly."
            ),
        )

        score_choice = st.selectbox(
            "Minimum required interaction score",
            list(SCORE_PRESETS),
            index=0,
        )
        required_score = SCORE_PRESETS[score_choice]
        if required_score is None:
            custom_score = st.number_input(
                "Custom confidence",
                min_value=0.0,
                max_value=0.999,
                value=0.900,
                step=0.001,
                format="%.3f",
            )
            required_score = int(round(float(custom_score) * 1000))

        st.markdown("**Max number of interactors to show**")
        first_shell = _shell_control(st, "1st shell", "string_first_shell")
        second_shell = _shell_control(st, "2nd shell", "string_second_shell")

    with st.expander("Advanced Settings", expanded=False):
        layout = st.radio(
            "Network layout",
            ["Force-directed", "Circular"],
            index=0,
            help=(
                "This controls the app-generated editable network. STRING's public image API "
                "does not expose its website layout switch."
            ),
        ).lower().replace("-", "_")

        st.markdown("**Network display options**")
        colorblind_friendly = st.checkbox("Colorblind-friendly local network", value=True)
        bubble_3d = st.checkbox("Enable 3D bubble design", value=True)
        block_structure_pics = st.checkbox(
            "Disable structure previews inside network bubbles",
            value=False,
        )
        center_node_labels = st.checkbox("Center protein names on nodes", value=False)
        show_query_node_labels = st.checkbox("Show your query protein names", value=False)
        hide_disconnected_nodes = st.checkbox("Hide disconnected nodes in the network", value=False)
        hide_node_labels = st.checkbox("Hide protein names", value=False)
        label_font_size = int(
            st.slider(
                "Protein name font size",
                min_value=5,
                max_value=50,
                value=12,
            )
        )

        typed_physical_edges = True
        typed_regulatory_edges = True
        show_regulatory_signs = True
        if network_flavor == "typed":
            st.markdown("**Typed overlay visibility**")
            typed_physical_edges = st.checkbox("Show physical overlay", value=True)
            typed_regulatory_edges = st.checkbox("Show regulatory overlay", value=True)
        if network_type == "regulatory" or network_flavor == "typed":
            show_regulatory_signs = st.checkbox(
                "Show positive/negative regulatory signs",
                value=True,
            )

    if not active_sources:
        st.error("Select at least one interaction source before running the analysis.")

    return {
        "species_label": species_label,
        "species": species,
        "string_version": string_version,
        "network_type": network_type,
        "network_type_label": network_type_label,
        "network_flavor": network_flavor,
        "active_sources": active_sources,
        "required_score": int(required_score),
        "first_shell": int(first_shell),
        "second_shell": int(second_shell),
        "add_nodes": int(first_shell + second_shell),
        "layout": layout,
        "colorblind_friendly": bool(colorblind_friendly),
        "bubble_3d": bool(bubble_3d),
        "block_structure_pics": bool(block_structure_pics),
        "center_node_labels": bool(center_node_labels),
        "show_query_node_labels": bool(show_query_node_labels),
        "hide_disconnected_nodes": bool(hide_disconnected_nodes),
        "hide_node_labels": bool(hide_node_labels),
        "label_font_size": int(label_font_size),
        "typed_physical_edges": bool(typed_physical_edges),
        "typed_regulatory_edges": bool(typed_regulatory_edges),
        "show_regulatory_signs": bool(show_regulatory_signs),
        "evidence_transfer": True,
    }
