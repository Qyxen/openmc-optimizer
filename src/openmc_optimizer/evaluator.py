"""Candidate preparation and OpenMC evaluation."""

from __future__ import annotations

import hashlib
import importlib.util
import json
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
) -> dict:
    """Load XML, add the optimizer-owned assembly power tally, run OpenMC, compute PPF."""
    try:
        import openmc
    except ImportError as error:
        raise RuntimeError("OpenMC Python bindings are required to evaluate candidates") from error

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
    # Assign an ID above the model's existing tallies to avoid collisions when
    # candidate models already define tally 1 (as this SMR model does).
    tallies = openmc.Tallies.from_xml(path=str(directory / "tallies.xml")) if (directory / "tallies.xml").exists() else openmc.Tallies()
    next_tally_id = max((tally.id for tally in tallies), default=0) + 1
    tally = openmc.Tally(name="optimizer_assembly_power", tally_id=next_tally_id)
    tally.filters = [openmc.CellInstanceFilter(instances)]
    tally.scores = ["fission-q-recoverable"]
    tally.estimator = "tracklength"
    tallies.append(tally)
    tallies.export_to_xml(path=str(directory / "tallies.xml"))

    subprocess.run(list(openmc_command), cwd=directory, check=True)
    statepoint_path = directory / statepoint_name
    if not statepoint_path.exists():
        raise FileNotFoundError(f"OpenMC did not produce {statepoint_path}")
    with openmc.StatePoint(str(statepoint_path)) as statepoint:
        values = statepoint.get_tally(name="optimizer_assembly_power").get_values(
            scores=["fission-q-recoverable"], value="mean"
        ).reshape(-1)
        keff = getattr(statepoint, "keff", None)
        if keff is None:
            # Compatibility with OpenMC versions before StatePoint.keff.
            keff = statepoint.k_combined
    if len(values) != len(instances) or any(value <= 0 for value in values):
        raise ValueError("Assembly tally is incomplete or contains nonpositive power values")
    mean_power = float(values.mean())
    max_index = int(values.argmax())
    return {
        "ppf": float(values[max_index] / mean_power),
        "keff": float(keff.nominal_value),
        "keff_std_dev": float(keff.std_dev),
        "assembly_count": len(values),
        "mean_assembly_power": mean_power,
        "max_assembly_power": float(values[max_index]),
        "max_power_assembly_index": max_index,
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
