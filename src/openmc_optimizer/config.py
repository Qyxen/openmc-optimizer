"""Configuration loading and validation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class Parameter:
    name: str
    lower: float
    upper: float


@dataclass(frozen=True)
class ProjectConfig:
    method: str
    source_file: Path
    input_dir: Path
    work_dir: Path
    adapter: Path
    adapter_options: dict[str, Any]
    openmc_command: tuple[str, ...]
    statepoint_file: str
    assembly_cell_ids: tuple[int, ...] | None
    target_ppf: float | None
    target_keff: float | None
    keff_tolerance: float
    ppf_weight: float
    keff_weight: float
    ppf_scale: float
    population_reduction_rate: float
    population_size: int
    generations: int
    elite_count: int
    mutation_rate: float
    mutation_scale: float
    seed: int
    parameters: tuple[Parameter, ...]
