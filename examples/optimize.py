"""Minimal Python entry point for optimizing an OpenMC model."""

from pathlib import Path

from openmc_optimizer import MaterialFraction, OpenMCOptimizer, VolumeFraction


def main() -> None:
    project_dir = Path(__file__).resolve().parent
    optimizer = OpenMCOptimizer(
        model_dir=project_dir / "model",  # Put geometry/materials/settings XML here.
        variables=[
            MaterialFraction(material_id=11, component="U235", bounds=(0.020, 0.050), name="A01_U235", normalization_group=("U235", "U238")),
            MaterialFraction(material_id=12, component="U235", bounds=(0.020, 0.050), name="A02_U235", normalization_group=("U235", "U238")),
            MaterialFraction(material_id=13, component="U235", bounds=(0.020, 0.050), name="B01_U235", normalization_group=("U235", "U238")),
            MaterialFraction(material_id=14, component="U235", bounds=(0.020, 0.050), name="B02_U235", normalization_group=("U235", "U238")),
            # Uncomment and replace IDs to optimize a Gd2O3 volume fraction too.
            # VolumeFraction(
            #     target_material_id=20,
            #     component_material_id=90,
            #     bounds=(0.0, 0.10),
            #     name="Gd2O3_loading",
            # ),
        ],
        target_ppf=1.35,
        target_keff=1.0,
        keff_tolerance=0.005,
        ppf_weight=0.7,
        keff_weight=0.3,
        ppf_scale=0.5,
        population_reduction_rate=0.20,
        work_dir=project_dir / "results",
        population_size=16,
        generations=25,
        elite_count=2,
        mutation_rate=0.20,
        mutation_scale=0.10,
        seed=14,
    )

    result = optimizer.run()
    print(f"Best PPF: {result.best_ppf:.5f}")
    print(f"Best keff: {result.best_keff:.5f} +/- {result.best_keff_std_dev:.5f}")
    print(f"Best weighted loss: {result.best_objective_loss:.5f}")
    print(f"Best values: {result.best_parameters}")
    print(f"Target reached: {result.target_reached}")
    print(f"Summary: {result.summary_file}")
    print(f"Generation history CSV: {result.history_csv}")
    print(f"Unique candidates CSV: {result.candidates_csv}")


if __name__ == "__main__":
    main()
