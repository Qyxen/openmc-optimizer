import xml.etree.ElementTree as ET

import pytest

from openmc_optimizer.material_adapter import apply_parameters


def _write_materials(path, material_xml):
    path.joinpath("materials.xml").write_text(
        f"<materials>{material_xml}</materials>", encoding="utf-8"
    )


def _material(path, material_id):
    root = ET.parse(path / "materials.xml").getroot()
    return root.find(f"material[@id='{material_id}']")


@pytest.mark.parametrize("basis", ["wo", "ao"])
def test_material_fraction_updates_component_and_rescales_remainder(tmp_path, basis):
    _write_materials(
        tmp_path,
        f'<material id="7"><density value="10" units="g/cm3"/>'
        f'<nuclide name="U235" {basis}="0.2"/>'
        f'<nuclide name="U238" {basis}="0.8"/></material>',
    )

    apply_parameters(
        {"u235": 0.3},
        tmp_path,
        {
            "bindings": {
                "u235": {
                    "type": "material_fraction",
                    "material_id": 7,
                    "component": "U235",
                    "component_type": "nuclide",
                    "basis": basis,
                }
            }
        },
    )

    material = _material(tmp_path, 7)
    fractions = {
        node.attrib["name"]: float(node.attrib[basis])
        for node in material
        if node.tag in ("nuclide", "element")
    }
    assert fractions == pytest.approx({"U235": 0.3, "U238": 0.7})


def test_material_fraction_adds_a_missing_component(tmp_path):
    _write_materials(
        tmp_path,
        '<material id="7"><density value="10" units="g/cm3"/>'
        '<nuclide name="U235" wo="0.2"/><nuclide name="U238" wo="0.8"/>'
        "</material>",
    )

    apply_parameters(
        {"np237": 0.1},
        tmp_path,
        {
            "bindings": {
                "np237": {
                    "type": "material_fraction",
                    "material_id": 7,
                    "component": "Np237",
                    "basis": "wo",
                }
            }
        },
    )

    material = _material(tmp_path, 7)
    fractions = {
        node.attrib["name"]: float(node.attrib["wo"])
        for node in material
        if node.tag in ("nuclide", "element")
    }
    assert fractions == pytest.approx({"U235": 0.18, "U238": 0.72, "Np237": 0.1})


def test_material_fraction_can_vary_within_a_normalization_group(tmp_path):
    _write_materials(
        tmp_path,
        '<material id="7"><density value="10" units="g/cm3"/>'
        '<nuclide name="U235" wo="0.1"/><nuclide name="U238" wo="0.7"/>'
        '<nuclide name="Pu239" wo="0.1"/><nuclide name="O16" wo="0.1"/>'
        "</material>",
    )

    apply_parameters(
        {"u235_enrichment": 0.04},
        tmp_path,
        {
            "bindings": {
                "u235_enrichment": {
                    "type": "material_fraction",
                    "material_id": 7,
                    "component": "U235",
                    "component_type": "nuclide",
                    "basis": "wo",
                    "normalization_group": ["U235", "U238"],
                }
            }
        },
    )

    material = _material(tmp_path, 7)
    fractions = {
        node.attrib["name"]: float(node.attrib["wo"])
        for node in material
        if node.tag in ("nuclide", "element")
    }
    assert fractions == pytest.approx(
        {"U235": 0.032, "U238": 0.768, "Pu239": 0.1, "O16": 0.1}
    )


def test_volume_fraction_mixes_composition_and_density(tmp_path):
    _write_materials(
        tmp_path,
        '<material id="1"><density value="2" units="g/cm3"/>'
        '<nuclide name="H1" wo="1.0"/></material>'
        '<material id="2"><density value="4" units="g/cm3"/>'
        '<nuclide name="O16" wo="1.0"/></material>',
    )

    apply_parameters(
        {"mix": 0.25},
        tmp_path,
        {
            "bindings": {
                "mix": {
                    "type": "volume_fraction_mix",
                    "target_material_id": 1,
                    "component_material_id": 2,
                }
            }
        },
    )

    target = _material(tmp_path, 1)
    assert float(target.find("density").attrib["value"]) == pytest.approx(2.5)
    fractions = {
        node.attrib["name"]: float(node.attrib["wo"])
        for node in target
        if node.tag in ("nuclide", "element")
    }
    assert fractions == pytest.approx({"H1": 0.6, "O16": 0.4})


def test_volume_fraction_requires_density_in_grams_per_cubic_centimeter(tmp_path):
    _write_materials(
        tmp_path,
        '<material id="1"><density value="2" units="g/cm3"/>'
        '<nuclide name="H1" wo="1.0"/></material>'
        '<material id="2"><density value="4" units="atom/b-cm"/>'
        '<nuclide name="O16" wo="1.0"/></material>',
    )

    with pytest.raises(ValueError, match="explicit g/cm3 density"):
        apply_parameters(
            {"mix": 0.25},
            tmp_path,
            {
                "bindings": {
                    "mix": {
                        "type": "volume_fraction_mix",
                        "target_material_id": 1,
                        "component_material_id": 2,
                    }
                }
            },
        )
