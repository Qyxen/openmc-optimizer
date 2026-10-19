import pytest

from openmc_optimizer import MaterialFraction, OpenMCOptimizer, VolumeFraction


def test_material_fraction_defaults_to_atom_fraction():
    variable = MaterialFraction(
        material_id=1,
        component="U235",
        bounds=(0.01, 0.05),
    )

    assert variable.basis == "ao"
    assert variable.binding()["basis"] == "ao"


def test_material_fraction_requires_bounds_within_zero_and_one(tmp_path):
    with pytest.raises(ValueError, match="MaterialFraction bounds"):
        OpenMCOptimizer(
            model_dir=tmp_path,
            variables=[
                MaterialFraction(
                    material_id=1,
                    component="U235",
                    bounds=(0.0, 1.1),
                )
            ],
        )


def test_volume_fraction_requires_bounds_within_zero_and_one(tmp_path):
    with pytest.raises(ValueError, match="VolumeFraction bounds"):
        OpenMCOptimizer(
            model_dir=tmp_path,
            variables=[
                VolumeFraction(
                    target_material_id=1,
                    component_material_id=2,
                    bounds=(-0.1, 0.5),
                )
            ],
        )


def test_optimizer_requires_unique_parameter_names(tmp_path):
    with pytest.raises(ValueError, match="names must be unique"):
        OpenMCOptimizer(
            model_dir=tmp_path,
            variables=[
                MaterialFraction(1, "U235", (0.01, 0.05), name="same"),
                MaterialFraction(2, "U235", (0.01, 0.05), name="same"),
            ],
        )
