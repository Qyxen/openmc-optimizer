"""Small Python API for OpenMC model optimization."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias

from .config import Parameter, ProjectConfig
from .evaluator import OpenMCEvaluator
from .ga import optimize


@dataclass(frozen=True)
class MaterialFraction:
    """Optimize one nuclide or element fraction in any OpenMC material.

    The selected component is set to the candidate fraction and the other
    composition entries, or optional normalization group, are scaled
    proportionally to fill the remainder.
    Fractions use the ``wo`` (weight) or ``ao`` (atom) basis from materials.xml.
    """

    material_id: int
    component: str
    bounds: tuple[float, float]
    name: str | None = None
    component_type: str = "nuclide"
    basis: str = "wo"
    normalization_group: tuple[str, ...] | None = None

    @property
    def parameter_name(self) -> str:
        return self.name or f"{self.component}_{self.basis}_{self.material_id}"

    def binding(self) -> dict:
        binding = {
            "type": "material_fraction",
            "material_id": self.material_id,
            "component": self.component,
            "component_type": self.component_type,
            "basis": self.basis,
        }
        if self.normalization_group is not None:
            binding["normalization_group"] = self.normalization_group
        return binding


@dataclass(frozen=True)
class VolumeFraction:
    """Optimize volume fraction of one material mixed into another."""

    target_material_id: int
    component_material_id: int
    bounds: tuple[float, float]
    name: str | None = None

    @property
    def parameter_name(self) -> str:
        return self.name or f"volume_fraction_{self.target_material_id}_{self.component_material_id}"

    def binding(self) -> dict:
        return {
            "type": "volume_fraction_mix",
            "target_material_id": self.target_material_id,
            "component_material_id": self.component_material_id,
        }


DesignVariable: TypeAlias = MaterialFraction | VolumeFraction


@dataclass(frozen=True)
class OptimizationResult:
    best_parameters: dict[str, float]
    best_ppf: float
    best_keff: float
    best_keff_std_dev: float
    best_objective_loss: float | None
    target_ppf: float | None
    target_keff: float | None
    keff_tolerance: float
    target_reached: bool | None
    keff_target_reached: bool | None
    history: list[dict]
    summary_file: Path
    history_csv: Path
    candidates_csv: Path


class OpenMCOptimizer:
    """Minimize assembly PPF, optionally subject to a target k-effective.

    Parameters are provided as ``MaterialFraction`` or ``VolumeFraction``
    objects. The optimizer copies XML inputs for every
    candidate, edits the material recipe, adds an assembly-power tally, runs
    OpenMC, and minimizes a weighted PPF/keff loss. If ``target_keff`` is set,
    candidates must be within ``keff_tolerance`` before the loss ranks them.
    """

    def __init__(
        self,
        *,
        model_dir: str | Path,
        variables: list[DesignVariable],
        assembly_cell_ids: list[int] | tuple[int, ...] | None = None,
        target_ppf: float | None = None,
        target_keff: float | None = None,
        keff_tolerance: float = 0.005,
        ppf_weight: float = 0.7,
        keff_weight: float = 0.3,
        ppf_scale: float = 0.5,
        population_reduction_rate: float = 0.20,
        work_dir: str | Path = "results",
        population_size: int = 20,
        generations: int = 30,
        elite_count: int = 2,
        mutation_rate: float = 0.15,
        mutation_scale: float = 0.12,
        seed: int = 1,
        openmc_command: list[str] | tuple[str, ...] = ("openmc",),
        statepoint_file: str = "statepoint.100.h5",
    ):
        if not variables:
            raise ValueError("Provide at least one design variable")
        if assembly_cell_ids is not None and (
            not assembly_cell_ids or any(cell_id <= 0 for cell_id in assembly_cell_ids)
        ):
            raise ValueError("assembly_cell_ids must contain positive cell IDs")
        if len({variable.parameter_name for variable in variables}) != len(variables):
            raise ValueError("Design variable names must be unique")
        if population_size < 10 or not 1 <= elite_count < population_size:
            raise ValueError("Require initial population_size >= 10 and 1 <= elite_count < population_size")
        if generations < 1 or not 0.0 <= mutation_rate <= 1.0 or mutation_scale < 0.0:
            raise ValueError("Invalid genetic algorithm generation or mutation settings")
        if target_ppf is not None and target_ppf <= 0.0:
            raise ValueError("target_ppf must be positive")
        if target_keff is not None and target_keff <= 0.0:
            raise ValueError("target_keff must be positive")
        if keff_tolerance < 0.0:
            raise ValueError("keff_tolerance must be non-negative")
        if ppf_weight < 0.0 or keff_weight < 0.0 or ppf_weight + keff_weight <= 0.0:
            raise ValueError("Objective weights must be non-negative with a positive sum")
        if target_keff is None and ppf_weight == 0.0:
            raise ValueError("ppf_weight must be positive when target_keff is not set")
        if ppf_scale <= 0.0:
            raise ValueError("ppf_scale must be positive")
        if not 0.0 <= population_reduction_rate < 1.0:
            raise ValueError("population_reduction_rate must be in [0, 1)")

        parameters = []
        bindings = {}
        for variable in variables:
            lower, upper = map(float, variable.bounds)
            if lower >= upper:
                raise ValueError(f"Invalid bounds for {variable.parameter_name}")
            if isinstance(variable, MaterialFraction):
                if variable.material_id <= 0:
                    raise ValueError("MaterialFraction material_id must be positive")
                if variable.component_type not in ("nuclide", "element"):
                    raise ValueError("component_type must be 'nuclide' or 'element'")
                if variable.basis not in ("wo", "ao"):
                    raise ValueError("MaterialFraction basis must be 'wo' or 'ao'")
                if variable.normalization_group is not None and (
                    not variable.normalization_group
                    or variable.component not in variable.normalization_group
                    or len(set(variable.normalization_group)) != len(variable.normalization_group)
                ):
                    raise ValueError(
                        "normalization_group must contain unique component names, "
                        "including the selected component"
                    )
                if lower < 0.0 or upper > 1.0:
                    raise ValueError("MaterialFraction bounds must be within [0, 1]")
            elif isinstance(variable, VolumeFraction):
                if variable.target_material_id <= 0 or variable.component_material_id <= 0:
                    raise ValueError("VolumeFraction material IDs must be positive")
                if lower < 0.0 or upper > 1.0:
                    raise ValueError("VolumeFraction bounds must be within [0, 1]")
            parameters.append(Parameter(variable.parameter_name, lower, upper))
            bindings[variable.parameter_name] = variable.binding()

        self.config = ProjectConfig(
            source_file=Path.cwd() / "python-api",
            input_dir=Path(model_dir).expanduser().resolve(),
            work_dir=Path(work_dir).expanduser().resolve(),
            adapter=Path(__file__).with_name("material_adapter.py"),
            adapter_options={"bindings": bindings},
            openmc_command=tuple(str(item) for item in openmc_command),
            statepoint_file=statepoint_file,
            assembly_cell_ids=(
                None if assembly_cell_ids is None
                else tuple(int(item) for item in assembly_cell_ids)
            ),
            target_ppf=target_ppf,
            target_keff=target_keff,
            keff_tolerance=keff_tolerance,
            ppf_weight=ppf_weight,
            keff_weight=keff_weight,
            ppf_scale=ppf_scale,
            population_reduction_rate=population_reduction_rate,
            population_size=population_size,
            generations=generations,
            elite_count=elite_count,
            mutation_rate=mutation_rate,
            mutation_scale=mutation_scale,
            seed=seed,
            parameters=tuple(parameters),
        )

    def run(self) -> OptimizationResult:
        """Run the genetic algorithm and return its best candidate."""
        config = self.config
        if not config.input_dir.is_dir():
            raise FileNotFoundError(f"OpenMC model directory not found: {config.input_dir}")

        initial_data = [
            ("Model directory", str(config.input_dir)),
            ("Work directory", str(config.work_dir)),
            ("Target PPF", "Not set" if config.target_ppf is None else f"{config.target_ppf:g}"),
            ("Target keff", "Not set" if config.target_keff is None else f"{config.target_keff:g}"),
            ("keff tolerance", f"+/- {config.keff_tolerance:g}"),
            ("PPF objective weight", f"{config.ppf_weight:g}"),
            ("keff objective weight", f"{config.keff_weight:g}"),
            ("PPF normalization scale", f"{config.ppf_scale:g}"),
            ("Initial population size", str(config.population_size)),
            ("Minimum population size", "10"),
            ("Population reduction per generation", f"{config.population_reduction_rate:.0%}"),
            ("Generations", str(config.generations)),
            ("Elite count", str(config.elite_count)),
            ("Mutation rate", f"{config.mutation_rate:g}"),
            ("Mutation scale", f"{config.mutation_scale:g}"),
            ("Random seed", str(config.seed)),
            ("Statepoint file", config.statepoint_file),
        ]
        initial_data.extend(
            (f"{parameter.name} bounds", f"[{parameter.lower:g}, {parameter.upper:g}]")
            for parameter in config.parameters
        )
        key_width = max(len("Setting"), *(len(key) for key, _ in initial_data))
        value_width = max(len("Value"), *(len(value) for _, value in initial_data))
        border = f"+-{'-' * key_width}-+-{'-' * value_width}-+"
        print("\nOptimization input data:")
        print(border)
        print(f"| {'Setting'.ljust(key_width)} | {'Value'.ljust(value_width)} |")
        print(border)
        for key, value in initial_data:
            print(f"| {key.ljust(key_width)} | {value.ljust(value_width)} |")
        print(border)

        evaluator = OpenMCEvaluator(config)
        best_parameters, best_result, history = optimize(
            config.parameters,
            evaluator,
            population_size=config.population_size,
            generations=config.generations,
            elite_count=config.elite_count,
            mutation_rate=config.mutation_rate,
            mutation_scale=config.mutation_scale,
            seed=config.seed,
            target_ppf=config.target_ppf,
            target_keff=config.target_keff,
            keff_tolerance=config.keff_tolerance,
            ppf_weight=config.ppf_weight,
            keff_weight=config.keff_weight,
            ppf_scale=config.ppf_scale,
            population_reduction_rate=config.population_reduction_rate,
        )
        keff_target_reached = (
            None if config.target_keff is None
            else abs(best_result["keff"] - config.target_keff) <= config.keff_tolerance
        )
        target_reached = None
        if config.target_ppf is not None or config.target_keff is not None:
            ppf_reached = (
                config.target_ppf is None or best_result["ppf"] <= config.target_ppf
            )
            keff_reached = keff_target_reached is not False
            target_reached = ppf_reached and keff_reached
        summary = {
            "best_parameters": best_parameters,
            "best_result": best_result,
            "target_ppf": config.target_ppf,
            "target_keff": config.target_keff,
            "keff_tolerance": config.keff_tolerance,
            "objective_weights": {
                "ppf": config.ppf_weight,
                "keff": config.keff_weight,
            },
            "ppf_scale": config.ppf_scale,
            "target_reached": target_reached,
            "keff_target_reached": keff_target_reached,
            "history": history,
        }
        config.work_dir.mkdir(parents=True, exist_ok=True)
        summary_file = config.work_dir / "optimization_summary.json"

        parameter_columns = [f"parameter_{parameter.name}" for parameter in config.parameters]
        history_csv = config.work_dir / "optimization_history.csv"
        with history_csv.open("w", newline="", encoding="utf-8") as stream:
            fieldnames = [
                "generation", "best_ppf", "mean_ppf", "best_keff", "mean_keff",
                "population_size", "keff_within_tolerance", "best_objective_loss",
                "mean_objective_loss",
                *parameter_columns,
            ]
            writer = csv.DictWriter(stream, fieldnames=fieldnames)
            writer.writeheader()
            for row in history:
                csv_row = {
                    "generation": row["generation"],
                    "population_size": row["population_size"],
                    "best_ppf": row["best_ppf"],
                    "mean_ppf": row["mean_ppf"],
                    "best_keff": row["best_keff"],
                    "mean_keff": row["mean_keff"],
                    "keff_within_tolerance": row["keff_within_tolerance"],
                    "best_objective_loss": row["best_objective_loss"],
                    "mean_objective_loss": row["mean_objective_loss"],
                }
                csv_row.update({
                    f"parameter_{name}": value
                    for name, value in row["best_parameters"].items()
                })
                writer.writerow(csv_row)

        candidates_csv = config.work_dir / "optimization_candidates.csv"
        with candidates_csv.open("w", newline="", encoding="utf-8") as stream:
            fieldnames = [
                "generation",
                "population_index",
                "ppf",
                "keff",
                "keff_std_dev",
                "objective_loss",
                "assembly_count",
                "mean_assembly_power",
                "max_assembly_power",
                "max_power_assembly_index",
                "assembly_cell_ids",
                *parameter_columns,
            ]
            writer = csv.DictWriter(stream, fieldnames=fieldnames)
            writer.writeheader()
            for candidate_result in evaluator.cache.values():
                csv_row = {
                    key: candidate_result.get(key)
                    for key in fieldnames
                    if not key.startswith("parameter_")
                }
                csv_row["assembly_cell_ids"] = json.dumps(
                    candidate_result.get("assembly_cell_ids", [])
                )
                csv_row.update({
                    f"parameter_{name}": value
                    for name, value in candidate_result.get("candidate", {}).items()
                })
                writer.writerow(csv_row)

        summary["csv_files"] = {
            "history": str(history_csv),
            "unique_candidates": str(candidates_csv),
        }
        summary_file.write_text(json.dumps(summary, indent=2) + "\n")

        return OptimizationResult(
            best_parameters=best_parameters,
            best_ppf=float(best_result["ppf"]),
            best_keff=float(best_result["keff"]),
            best_keff_std_dev=float(best_result["keff_std_dev"]),
            best_objective_loss=best_result.get("objective_loss"),
            target_ppf=config.target_ppf,
            target_keff=config.target_keff,
            keff_tolerance=config.keff_tolerance,
            target_reached=target_reached,
            keff_target_reached=keff_target_reached,
            history=history,
            summary_file=summary_file,
            history_csv=history_csv,
            candidates_csv=candidates_csv,
        )
