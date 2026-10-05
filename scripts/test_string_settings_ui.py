from __future__ import annotations

from pathlib import Path

from streamlit.testing.v1 import AppTest


APP = Path(__file__).resolve().parents[1] / "app.py"


def by_label(items, label):
    return next(item for item in items if item.label == label)


def main() -> None:
    at = AppTest.from_file(APP).run(timeout=30)
    assert not list(at.exception), list(at.exception)

    network_type = by_label(at.sidebar.radio, "Network type")
    assert network_type.value == "Full STRING network"
    applied = at.session_state["applied_control_settings"]
    assert applied["string_options"]["network_type"] == "functional"

    by_label(at.sidebar.radio, "Meaning of network edges")
    by_label(at.sidebar.selectbox, "Minimum required interaction score")
    by_label(at.sidebar.selectbox, "1st shell")
    by_label(at.sidebar.selectbox, "2nd shell")
    by_label(at.sidebar.radio, "Network layout")
    by_label(at.sidebar.checkbox, "Textmining")
    by_label(at.sidebar.checkbox, "Experiments")
    by_label(at.sidebar.checkbox, "Databases")
    transfer = by_label(at.sidebar.checkbox, "Evidence transfer")
    assert transfer.disabled

    network_type.set_value("Regulatory subnetwork").run(timeout=30)
    assert not list(at.exception), list(at.exception)

    # Editing the fragment must not immediately change the settings used by
    # the main scientific page.
    applied = at.session_state["applied_control_settings"]
    assert applied["string_options"]["network_type"] == "functional"

    flavor = by_label(at.sidebar.radio, "Meaning of network edges")
    assert "Typed" not in flavor.options
    by_label(at.sidebar.checkbox, "Show positive/negative regulatory signs")

    apply_button = next(
        button for button in at.sidebar.button
        if button.label == "Apply settings"
    )
    assert not apply_button.disabled
    apply_button.click().run(timeout=30)
    assert not list(at.exception), list(at.exception)
    applied = at.session_state["applied_control_settings"]
    assert applied["string_options"]["network_type"] == "regulatory"

    network_type = by_label(at.sidebar.radio, "Network type")
    network_type.set_value("Full STRING network").run(timeout=30)
    flavor = by_label(at.sidebar.radio, "Meaning of network edges")
    flavor.set_value("Typed").run(timeout=30)
    assert not list(at.exception), list(at.exception)
    by_label(at.sidebar.checkbox, "Show physical overlay")
    by_label(at.sidebar.checkbox, "Show regulatory overlay")

    print("STRING_UI_SMOKE_OK")


if __name__ == "__main__":
    main()
