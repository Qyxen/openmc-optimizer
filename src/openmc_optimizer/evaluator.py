"""Candidate preparation and OpenMC evaluation."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

from .config import ProjectConfig


class OpenMCEvaluator:
    def __init__(self, config: ProjectConfig):
        self.config = config
        self.cache: dict[str, dict] = {}
        self.config.work_dir.mkdir(parents=True, exist_ok=True)
        self.adapter = self._load_adapter(config.adapter)
        xml_files = list(config.input_dir.glob("*.xml"))
        if not xml_files:
            raise FileNotFoundError(f"No XML model files found in {config.input_dir}")
        self.xml_files = xml_files

    @staticmethod
    def _load_adapter(path: Path):
        spec = importlib.util.spec_from_file_location("openmc_candidate_adapter", path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Could not load candidate adapter {path}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        if not hasattr(module, "apply_parameters"):
            raise AttributeError("Candidate adapter must define apply_parameters(values, model_dir, options)")
        return module

    @staticmethod
    def _key(values: dict[str, float]) -> str:
        packed = json.dumps(values, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(packed.encode()).hexdigest()[:12]

    def __call__(self, values: dict[str, float], generation: int, index: int) -> dict:
        key = self._key(values)
        if key in self.cache:
            return self.cache[key]

        candidate_dir = self.config.work_dir / f"candidate_{key}"
        candidate_dir.mkdir(parents=True, exist_ok=True)
        for source in self.xml_files:
            shutil.copy2(source, candidate_dir / source.name)

        self.adapter.apply_parameters(values, candidate_dir, self.config.adapter_options)
        evaluate_script = Path(__file__).with_name("openmc_worker.py")
        command = [
            sys.executable,
            str(evaluate_script),
            "--directory", str(candidate_dir),
            "--statepoint", self.config.statepoint_file,
            "--assembly-cell-ids", (
                "" if self.config.assembly_cell_ids is None
                else ",".join(map(str, self.config.assembly_cell_ids))
            ),
            "--openmc-command", json.dumps(self.config.openmc_command),
            "--ppf-tally-name", self.config.ppf_tally_name or "",
            "--ppf-tally-score", self.config.ppf_tally_score,
        ]
        subprocess.run(command, check=True)
        ppf_path = candidate_dir / "ppf_result.json"
        result = json.loads(ppf_path.read_text())
        result["candidate"] = values
        result["generation"] = generation
        result["population_index"] = index
        (candidate_dir / "candidate.json").write_text(json.dumps(result, indent=2) + "\n")
        self.cache[key] = result
        return result


def evaluate_openmc(
    directory: Path,
    statepoint_name: str,
    assembly_cell_ids: list[int] | None,
    openmc_command: tuple[str, ...],
    ppf_tally_name: str | None = None,
    ppf_tally_score: str = "fission-q-recoverable",
) -> dict:
    """Run OpenMC and compute PPF from built-in or user-supplied tally bins."""
    try:
        import openmc
    except ImportError as error:
        raise RuntimeError("OpenMC Python bindings are required to evaluate candidates") from error

    tallies = openmc.Tallies.from_xml(path=str(directory / "tallies.xml")) if (directory / "tallies.xml").exists() else openmc.Tallies()
    if ppf_tally_name is None:
        geometry = openmc.Geometry.from_xml(path=str(directory / "geometry.xml"))
        geometry.determine_paths(instances_only=True)
        all_cells = geometry.get_all_cells()
        if assembly_cell_ids is None:
            assembly_cell_ids = infer_assembly_cell_ids(geometry, openmc)
        instances = []
        for cell_id in assembly_cell_ids:
            if cell_id not in all_cells:
                raise ValueError(f"Assembly root cell ID {cell_id} is absent from geometry.xml")
            cell = all_cells[cell_id]
            instances.extend((cell, instance) for instance in range(cell.num_instances))
        if len(instances) < 2:
            raise ValueError("PPF requires at least two assembly instances")

        # Tally's automatic ID counter is independent of IDs loaded from tallies.xml.
        # Assign an ID above the model's existing tallies to avoid collisions.
        next_tally_id = max((tally.id for tally in tallies), default=0) + 1
        tally = openmc.Tally(name="optimizer_assembly_power", tally_id=next_tally_id)
        tally.filters = [openmc.CellInstanceFilter(instances)]
        tally.scores = [ppf_tally_score]
        tally.estimator = "tracklength"
        tallies.append(tally)
        tally_name = "optimizer_assembly_power"
    else:
        tally_name = ppf_tally_name
        try:
            tally = tallies.get_by_name(tally_name)
        except KeyError as error:
            raise ValueError(
                f"User PPF tally {tally_name!r} was not found in tallies.xml"
            ) from error
        if ppf_tally_score not in tally.scores:
            raise ValueError(
                f"User PPF tally {tally_name!r} does not include score "
                f"{ppf_tally_score!r}"
            )
        region_filter_types = {
            "CellFilter",
            "CellInstanceFilter",
            "DistribcellFilter",
            "MaterialFilter",
            "MeshFilter",
            "UniverseFilter",
        }
        if (
            len(tally.filters) != 1
            or type(tally.filters[0]).__name__ not in region_filter_types
            or list(tally.nuclides) != ["total"]
        ):
            raise ValueError(
                "A custom PPF tally must have exactly one supported spatial "
                "region filter (cell, cell instance, distribcell, material, "
                "mesh, or universe) and one total nuclide bin. Each bin must "
                "represent one spatial region's integrated power."
            )

    tallies.export_to_xml(path=str(directory / "tallies.xml"))

    subprocess.run(list(openmc_command), cwd=directory, check=True)
    statepoint_path = directory / statepoint_name
    if not statepoint_path.exists():
        raise FileNotFoundError(f"OpenMC did not produce {statepoint_path}")
    with openmc.StatePoint(str(statepoint_path)) as statepoint:
        power_values = statepoint.get_tally(name=tally_name).get_values(
            scores=[ppf_tally_score], value="mean"
        ).reshape(-1)
        keff = getattr(statepoint, "keff", None)
        if keff is None:
            # Compatibility with OpenMC versions before StatePoint.keff.
            keff = statepoint.k_combined
    if len(power_values) < 2 or any(
        not math.isfinite(float(value)) or float(value) <= 0.0
        for value in power_values
    ):
        raise ValueError("PPF tally must contain at least two positive regional power values")
    mean_power = float(power_values.mean())
    max_index = int(power_values.argmax())
    return {
        "ppf": float(power_values[max_index] / mean_power),
        "keff": float(keff.nominal_value),
        "keff_std_dev": float(keff.std_dev),
        "assembly_count": len(power_values),
        "power_region_count": len(power_values),
        "ppf_tally_name": tally_name,
        "ppf_tally_score": ppf_tally_score,
        "mean_assembly_power": mean_power,
        "mean_region_power": mean_power,
        "max_assembly_power": float(power_values[max_index]),
        "max_region_power": float(power_values[max_index]),
        "max_power_assembly_index": max_index,
        "max_power_region_index": max_index,
        "assembly_cell_ids": assembly_cell_ids,
    }


def infer_assembly_cell_ids(geometry, openmc) -> list[int]:
    """Find cells filled by the deepest nested lattices in common core models."""
    candidates: list[tuple[int, int]] = []
    visited: set[tuple[int, int]] = set()

    def lattice_universes(items):
        if isinstance(items, openmc.Universe):
            yield items
        elif isinstance(items, (list, tuple)):
            for item in items:
                yield from lattice_universes(item)
        elif hasattr(items, "flat"):
            for item in items.flat:
                yield from lattice_universes(item)

    def visit_universe(universe, depth):
        key = (universe.id, depth)
        if key in visited:
            return
        visited.add(key)
        for cell in universe.cells.values():
            fill = cell.fill
            if isinstance(fill, openmc.Lattice):
                candidates.append((depth, cell.id))
                for child_universe in lattice_universes(fill.universes):
                    visit_universe(child_universe, depth + 1)
            elif isinstance(fill, openmc.Universe):
                visit_universe(fill, depth)

    visit_universe(geometry.root_universe, 0)
    if not candidates:
        raise ValueError(
            "Could not infer assembly cells from the geometry. Provide "
            "assembly_cell_ids in the Python API."
        )
    deepest_level = max(depth for depth, _ in candidates)
    if deepest_level == 0:
        raise ValueError(
            "The geometry has only one lattice level, so assembly cells are "
            "ambiguous. Provide assembly_cell_ids explicitly."
        )
    return sorted({cell_id for depth, cell_id in candidates if depth == deepest_level})
