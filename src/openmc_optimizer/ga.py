"""A small bounded real-valued genetic algorithm."""

from __future__ import annotations

import random
from collections.abc import Callable

from .config import Parameter

MIN_POPULATION_SIZE = 10


def optimize(
    parameters: tuple[Parameter, ...],
    evaluate: Callable[[dict[str, float], int, int], dict],
    *,
    population_size: int,
    generations: int,
    elite_count: int,
    mutation_rate: float,
    mutation_scale: float,
    seed: int,
    target_ppf: float | None,
    target_keff: float | None = None,
    keff_tolerance: float = 0.005,
    ppf_weight: float = 0.7,
    keff_weight: float = 0.3,
    ppf_scale: float = 0.5,
    population_reduction_rate: float = 0.20,
) -> tuple[dict[str, float], dict, list[dict]]:
    rng = random.Random(seed)
    names = [parameter.name for parameter in parameters]
    population = [
        {p.name: rng.uniform(p.lower, p.upper) for p in parameters}
        for _ in range(population_size)
    ]
    # Keep one center point so each run samples a predictable initial candidate.
    population[0] = {p.name: (p.lower + p.upper) / 2.0 for p in parameters}
    minimum_population_size = min(MIN_POPULATION_SIZE, population_size)
    history = []
    best_candidate = None
    best_result = None

    for generation in range(generations):
        _print_population(parameters, population, generation + 1, generations)
        scored = [(evaluate(candidate, generation, index), candidate)
                  for index, candidate in enumerate(population)]
        scored.sort(key=lambda item: _fitness_key(
            item[0], target_keff, keff_tolerance,
            ppf_weight, keff_weight, ppf_scale,
        ))
        for result, _ in scored:
            result["objective_loss"] = _weighted_loss(
                result, target_keff, keff_tolerance,
                ppf_weight, keff_weight, ppf_scale,
            )
        generation_result, generation_best = scored[0]
        feasible_losses = [
            result["objective_loss"] for result, _ in scored
            if result["objective_loss"] is not None
        ]
        keff_values = [result.get("keff") for result, _ in scored]
        history.append({
            "generation": generation,
            "population_size": len(population),
            "best_ppf": generation_result["ppf"],
            "mean_ppf": sum(result["ppf"] for result, _ in scored) / len(scored),
            "best_keff": generation_result.get("keff"),
            "mean_keff": (
                sum(keff_values) / len(keff_values)
                if all(value is not None for value in keff_values) else None
            ),
            "keff_within_tolerance": _keff_is_feasible(
                generation_result, target_keff, keff_tolerance
            ),
            "best_objective_loss": generation_result["objective_loss"],
            "mean_objective_loss": (
                sum(feasible_losses) / len(feasible_losses) if feasible_losses else None
            ),
            "best_parameters": generation_best,
        })
        if best_result is None or _fitness_key(
            generation_result, target_keff, keff_tolerance,
            ppf_weight, keff_weight, ppf_scale,
        ) < _fitness_key(
            best_result, target_keff, keff_tolerance,
            ppf_weight, keff_weight, ppf_scale,
        ):
            best_result = generation_result
            best_candidate = dict(generation_best)
        status = f"best PPF={generation_result['ppf']:.6f}"
        if target_keff is not None:
            status += f", keff={generation_result['keff']:.6f}"
            status += f" (target {target_keff:.6f} +/- {keff_tolerance:.6f})"
            if generation_result["objective_loss"] is not None:
                status += f", weighted loss={generation_result['objective_loss']:.6f}"
        print(
            f"Generation {generation + 1}/{generations} "
            f"(population={len(population)}): {status}"
        )

        ppf_reached = target_ppf is None or generation_result["ppf"] <= target_ppf
        keff_reached = _keff_is_feasible(
            generation_result, target_keff, keff_tolerance
        )
        if target_ppf is not None and ppf_reached and keff_reached:
            break
        if generation == generations - 1:
            break

        next_size = max(
            minimum_population_size,
            int(len(population) * (1.0 - population_reduction_rate)),
        )
        next_elite_count = min(elite_count, next_size - 1)
        next_population = [
            dict(candidate) for _, candidate in scored[:next_elite_count]
        ]

        def tournament() -> dict[str, float]:
            choices = rng.sample(scored, k=min(3, len(scored)))
            return dict(min(
                choices,
                key=lambda item: _fitness_key(
                    item[0], target_keff, keff_tolerance,
                    ppf_weight, keff_weight, ppf_scale,
                ),
            )[1])

        while len(next_population) < next_size:
            first = tournament()
            second = tournament()
            child = {}
            for parameter in parameters:
                name = parameter.name
                value = first[name] if rng.random() < 0.5 else second[name]
                if rng.random() < mutation_rate:
                    value += rng.gauss(0.0, mutation_scale * (parameter.upper - parameter.lower))
                child[name] = min(parameter.upper, max(parameter.lower, value))
            next_population.append(child)
        population = next_population

    return best_candidate, best_result, history


def _keff_is_feasible(result, target_keff, keff_tolerance):
    return target_keff is None or abs(result["keff"] - target_keff) <= keff_tolerance


def _fitness_key(result, target_keff, keff_tolerance,
                 ppf_weight, keff_weight, ppf_scale):
    """Require keff feasibility, then rank feasible candidates by weighted loss."""
    if target_keff is None:
        return (0, _weighted_loss(
            result, None, keff_tolerance, ppf_weight, keff_weight, ppf_scale
        ))
    deviation = abs(result["keff"] - target_keff)
    if deviation <= keff_tolerance:
        return (0, _weighted_loss(
            result, target_keff, keff_tolerance,
            ppf_weight, keff_weight, ppf_scale,
        ))
    return (1, deviation, result["ppf"])


def _weighted_loss(result, target_keff, keff_tolerance,
                   ppf_weight, keff_weight, ppf_scale):
    """Compute normalized weighted loss for feasible candidates."""
    ppf_score = max(0.0, result["ppf"] - 1.0) / ppf_scale
    if target_keff is None:
        return ppf_score
    deviation = abs(result["keff"] - target_keff)
    if deviation > keff_tolerance:
        return None
    keff_score = deviation / keff_tolerance if keff_tolerance > 0.0 else 0.0
    weight_sum = ppf_weight + keff_weight
    return (ppf_weight * ppf_score + keff_weight * keff_score) / weight_sum


def _print_population(parameters, population, generation, generations) -> None:
    """Print the generation's candidates before evaluating them."""
    columns = ["Candidate", *(parameter.name for parameter in parameters)]
    rows = [
        [f"{index + 1}", *(f"{candidate[parameter.name]:.6g}" for parameter in parameters)]
        for index, candidate in enumerate(population)
    ]
    widths = [
        max(len(column), *(len(row[index]) for row in rows))
        for index, column in enumerate(columns)
    ]
    border = "+-" + "-+-".join("-" * width for width in widths) + "-+"

    def format_row(row):
        return "| " + " | ".join(value.ljust(width) for value, width in zip(row, widths)) + " |"

    print(f"\nGeneration {generation}/{generations} population (fractions):")
    print(border)
    print(format_row(columns))
    print(border)
    for row in rows:
        print(format_row(row))
    print(border)
