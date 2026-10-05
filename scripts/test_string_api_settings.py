from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import modules.string_api as string_api


def main() -> None:
    json_calls = []
    binary_calls = []
    text_calls = []

    def fake_json(method, payload, timeout=60, string_version="12.0"):
        json_calls.append(
            {
                "method": method,
                "payload": dict(payload),
                "string_version": string_version,
            }
        )
        return []

    def fake_binary(
        output_format,
        method,
        payload,
        timeout=90,
        string_version="12.0",
    ):
        binary_calls.append(
            {
                "format": output_format,
                "method": method,
                "payload": dict(payload),
                "string_version": string_version,
            }
        )
        return b"PNG"

    def fake_text(
        output_format,
        method,
        payload,
        timeout=90,
        string_version="12.0",
    ):
        text_calls.append(
            {
                "format": output_format,
                "method": method,
                "payload": dict(payload),
                "string_version": string_version,
            }
        )
        if method == "get_link":
            return "https://version-12.string-db.org/cgi/network"
        return "<svg></svg>"

    old_json = string_api._post_json
    old_binary = string_api._post_binary
    old_text = string_api._post_text
    old_sleep = string_api.time.sleep

    string_api._post_json = fake_json
    string_api._post_binary = fake_binary
    string_api._post_text = fake_text
    string_api.time.sleep = lambda _: None

    try:
        string_api.get_network(
            ["9606.A", "9606.B"],
            9606,
            700,
            "functional",
            network_flavor="typed",
            active_sources=["experiments"],
            add_nodes=15,
            typed_physical_edges=False,
            typed_regulatory_edges=True,
            show_query_node_labels=True,
            string_version="12.0",
        )
        payload = json_calls[-1]["payload"]
        assert payload["required_score"] == 700
        assert payload["network_type"] == "functional"
        assert payload["network_flavor"] == "typed"
        assert payload["typed_physical_edges"] == 0
        assert payload["typed_regulatory_edges"] == 1
        assert payload["add_nodes"] == 15
        assert payload["show_query_node_labels"] == 1

        media = string_api.get_network_media(
            ["9606.A", "9606.B"],
            9606,
            900,
            "functional",
            network_flavor="typed",
            first_shell=5,
            second_shell=10,
            active_sources=[
                "textmining",
                "experiments",
                "databases",
                "coexpression",
                "neighborhood",
                "fusion",
                "cooccurrence",
            ],
            typed_physical_edges=False,
            typed_regulatory_edges=True,
            show_regulatory_signs=False,
            bubble_3d=False,
            block_structure_pics=True,
            center_node_labels=True,
            show_query_node_labels=True,
            hide_disconnected_nodes=True,
            hide_node_labels=True,
            label_font_size=21,
            string_version="12.0",
        )
        assert media["highres_png"] == b"PNG"
        assert media["svg"] == b"<svg></svg>"
        assert media["link"]

        image_payload = binary_calls[-1]["payload"]
        expected = {
            "required_score": 900,
            "network_type": "functional",
            "network_flavor": "typed",
            "add_color_nodes": 5,
            "add_white_nodes": 10,
            "typed_physical_edges": 0,
            "typed_regulatory_edges": 1,
            "show_regulatory_signs": 0,
            "hide_node_labels": 1,
            "hide_disconnected_nodes": 1,
            "show_query_node_labels": 1,
            "block_structure_pics_in_bubbles": 1,
            "flat_node_design": 1,
            "center_node_labels": 1,
            "custom_label_font_size": 21,
        }
        for key, value in expected.items():
            assert image_payload.get(key) == value, (key, image_payload.get(key), value)

        svg_payload = next(
            call["payload"]
            for call in text_calls
            if call["method"] == "network"
        )
        for key, value in expected.items():
            assert svg_payload.get(key) == value, (key, svg_payload.get(key), value)

        link_payload = next(
            call["payload"]
            for call in text_calls
            if call["method"] == "get_link"
        )
        for key in (
            "network_flavor",
            "add_color_nodes",
            "add_white_nodes",
            "typed_physical_edges",
            "typed_regulatory_edges",
            "show_regulatory_signs",
            "hide_node_labels",
            "hide_disconnected_nodes",
            "show_query_node_labels",
            "block_structure_pics_in_bubbles",
        ):
            assert key in link_payload, key

        # These parameters are documented for image/SVG but not the stable-link API.
        assert "flat_node_design" not in link_payload
        assert "center_node_labels" not in link_payload
        assert "custom_label_font_size" not in link_payload

    finally:
        string_api._post_json = old_json
        string_api._post_binary = old_binary
        string_api._post_text = old_text
        string_api.time.sleep = old_sleep

    print("STRING_API_SETTINGS_OK")


if __name__ == "__main__":
    main()
