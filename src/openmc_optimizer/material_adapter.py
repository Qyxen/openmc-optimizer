"""Built-in adapter for changing OpenMC material compositions."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path
import re


def _element_atomic_weight(element: str) -> float:
    """Get an element's natural-abundance atomic weight from OpenMC."""
    try:
        import openmc
    except ImportError as error:
        raise RuntimeError(
            "OpenMC is required to convert atom fractions for element entries "
            "during volume-fraction mixing"
        ) from error

    try:
        return float(openmc.data.atomic_weight(element))
    except ValueError as error:
        raise ValueError(
            f"OpenMC has no natural atomic weight for element {element!r}; "
            "represent it with explicit nuclide entries instead"
        ) from error


def _species_weight_fractions(material: ET.Element) -> dict[tuple[str, str], float]:
    composition = {}
    total = 0.0
    for node in list(material):
        if node.tag not in ("nuclide", "element"):
            continue
        kind = node.tag
        name = node.attrib["name"]
        if "wo" in node.attrib:
            weight = float(node.attrib["wo"])
        elif "ao" in node.attrib:
            atom_fraction = float(node.attrib["ao"])
            if kind == "element":
                mass = _element_atomic_weight(name)
            else:
                match = re.match(r"^[A-Za-z]+(\d+)", name)
                if match is None:
                    raise ValueError(f"Invalid nuclide name: {name}")
                mass = float(match.group(1))
            weight = atom_fraction * mass
        else:
            raise ValueError(f"Composition entry {name} has neither ao nor wo")
        key = (kind, name)
        composition[key] = composition.get(key, 0.0) + weight
        total += weight
    if total <= 0.0:
        raise ValueError(f"Material {material.attrib.get('id')} has no composition")
    return {key: value / total for key, value in composition.items()}


def _material_density(material: ET.Element) -> float:
    density = material.find("density")
    if density is None or density.attrib.get("units") != "g/cm3":
        raise ValueError(
            f"Material {material.attrib.get('id')} needs an explicit g/cm3 density "
            "for volume-fraction mixing"
        )
    return float(density.attrib["value"])


def _apply_material_fraction(
    root: ET.Element,
    material_id: int,
    component: str,
    component_type: str,
    basis: str,
    fraction: float,
    normalization_group: tuple[str, ...] | None = None,
) -> None:
    """Set one material entry's fraction, scaling all other entries equally."""
    if basis not in ("ao", "wo"):
        raise ValueError("Material fraction basis must be 'ao' or 'wo'")
    if component_type not in ("nuclide", "element"):
        raise ValueError("Component type must be 'nuclide' or 'element'")
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("Material component fraction must be between 0 and 1")

    material = root.find(f"material[@id='{material_id}']")
    if material is None:
        raise ValueError(f"Material ID {material_id} was not found")
    entries = [node for node in list(material) if node.tag in ("nuclide", "element")]
    if not entries:
        raise ValueError(f"Material {material_id} has no nuclide or element entries")
    for node in entries:
        if basis not in node.attrib:
            raise ValueError(
                f"Material {material_id} must express every composition entry "
                f"using the requested {basis!r} basis"
            )

    selected = [
        node for node in entries
        if node.tag == component_type and node.attrib.get("name") == component
    ]

    if normalization_group is not None:
        group_names = set(normalization_group)
        if not group_names or component not in group_names:
            raise ValueError(
                "Normalization group must include the selected component"
            )
        group_entries = [
            node for node in entries
            if node.tag == component_type and node.attrib.get("name") in group_names
        ]
        present_names = {node.attrib["name"] for node in group_entries}
        missing_names = group_names - present_names
        if missing_names:
            raise ValueError(
                f"Material {material_id} is missing normalization-group entries: "
                f"{', '.join(sorted(missing_names))}"
            )
        selected_total = sum(float(node.attrib[basis]) for node in selected)
        group_total = sum(float(node.attrib[basis]) for node in group_entries)
        group_others = [node for node in group_entries if node not in selected]
        group_others_total = sum(float(node.attrib[basis]) for node in group_others)
        if group_total <= 0.0:
            raise ValueError(f"Normalization group in material {material_id} has zero fraction")
        if group_others_total <= 0.0 and fraction < 1.0:
            raise ValueError(
                f"Normalization group in material {material_id} has no other entries "
                "to scale to the remainder"
            )
        for node in selected:
            original = float(node.attrib[basis])
            share = original / selected_total if selected_total > 0.0 else 1.0 / len(selected)
            node.set(basis, f"{group_total * fraction * share:.12g}")
            node.attrib.pop("ao" if basis == "wo" else "wo", None)
        for node in group_others:
            original = float(node.attrib[basis])
            share = original / group_others_total if group_others_total > 0.0 else 0.0
            node.set(basis, f"{group_total * (1.0 - fraction) * share:.12g}")
            node.attrib.pop("ao" if basis == "wo" else "wo", None)
        return

    other = [node for node in entries if node not in selected]
    selected_total = sum(float(node.attrib[basis]) for node in selected)
    other_total = sum(float(node.attrib[basis]) for node in other)

    if other_total <= 0.0 and fraction < 1.0:
        raise ValueError(
            f"Material {material_id} has no other components to scale to the remainder"
        )
    if selected:
        for node in selected:
            original = float(node.attrib[basis])
            share = original / selected_total if selected_total > 0.0 else 1.0 / len(selected)
            node.set(basis, f"{fraction * share:.12g}")
            node.attrib.pop("ao" if basis == "wo" else "wo", None)
    else:
        selected = [ET.SubElement(material, component_type, {"name": component})]
        selected[0].set(basis, f"{fraction:.12g}")
    for node in other:
        original = float(node.attrib[basis])
        share = original / other_total if other_total > 0.0 else 0.0
        node.set(basis, f"{(1.0 - fraction) * share:.12g}")
        node.attrib.pop("ao" if basis == "wo" else "wo", None)


def _apply_volume_fraction_mix(
    root: ET.Element, target_id: int, component_id: int, volume_fraction: float
) -> None:
    if not 0.0 <= volume_fraction <= 1.0:
        raise ValueError("Volume fraction must be between 0 and 1")
    target = root.find(f"material[@id='{target_id}']")
    component = root.find(f"material[@id='{component_id}']")
    if target is None or component is None:
        raise ValueError(f"Could not find material IDs {target_id} and {component_id}")

    target_density = _material_density(target)
    component_density = _material_density(component)
    mixture_density = (1.0 - volume_fraction) * target_density + volume_fraction * component_density
    species_mass = {}
    for material, volume_share in ((target, 1.0 - volume_fraction), (component, volume_fraction)):
        density = _material_density(material)
        for key, weight_fraction in _species_weight_fractions(material).items():
            species_mass[key] = species_mass.get(key, 0.0) + volume_share * density * weight_fraction
    total_mass = sum(species_mass.values())

    for node in list(target):
        if node.tag in ("nuclide", "element"):
            target.remove(node)
    for (kind, name), mass in species_mass.items():
        ET.SubElement(target, kind, {"name": name, "wo": f"{mass / total_mass:.12g}"})
    density_node = target.find("density")
    density_node.set("value", f"{mixture_density:.12g}")


def apply_parameters(values: dict[str, float], model_dir: Path, options: dict) -> None:
    """Apply named composition and recipe parameters to materials.xml.

    The built-in operations are generic component fractions and physical
    volume mixing of two complete material recipes.
    """
    materials_path = model_dir / "materials.xml"
    tree = ET.parse(materials_path)
    root = tree.getroot()
    bindings = options.get("bindings", {})

    for parameter_name, value in values.items():
        if parameter_name not in bindings:
            raise KeyError(f"No adapter binding configured for {parameter_name}")
        binding = bindings[parameter_name]
        binding_type = binding.get("type", "material_fraction")
        if binding_type == "material_fraction":
            _apply_material_fraction(
                root,
                int(binding["material_id"]),
                str(binding["component"]),
                str(binding.get("component_type", "nuclide")),
                str(binding.get("basis", "ao")),
                float(value),
                (
                    None
                    if binding.get("normalization_group") is None
                    else tuple(binding["normalization_group"])
                ),
            )
            continue
        if binding_type == "volume_fraction_mix":
            _apply_volume_fraction_mix(
                root,
                int(binding["target_material_id"]),
                int(binding["component_material_id"]),
                float(value),
            )
            continue
        raise ValueError(f"Unsupported built-in material operation: {binding_type!r}")

    tree.write(materials_path, encoding="utf-8", xml_declaration=True)
