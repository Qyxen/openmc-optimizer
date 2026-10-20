"""A small bounded real-valued genetic algorithm."""

from __future__ import annotations

import random
from collections.abc import Callable
from math import ceil

from .config import Parameter

MIN_POPULATION_SIZE = 10
SURROGATE_EVALUATION_FRACTION = 0.5


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
    surrogate_evaluation_fraction: float = SURROGATE_EVALUATION_FRACTION,
) -> tuple[dict[str, float], dict, list[dict]]:
    """Run a bounded GA, using a surrogate to screen later generations.

    The initial population is evaluated in full. Later generations evaluate a
    configurable fraction of surrogate-ranked proposals; previous measured
    candidates remain eligible for survival. Duplicate points reuse cached
    results, and failed evaluations receive a worst-case score so the search
    can continue when other candidates succeed.
    """
    if not 0.0 < surrogate_evaluation_fraction <= 1.0:
        raise ValueError("surrogate_evaluation_fraction must be in (0, 1]")
    rng = random.Random(seed)
    population = [
        {p.name: rng.uniform(p.lower, p.upper) for p in parameters}
        for _ in range(population_size)
    ]
    # Keep one center point so each run samples a predictable initial candidate.
    population[0] = {p.name: (p.lower + p.upper) / 2.0 for p in parameters}
    minimum_population_size = min(MIN_POPULATION_SIZE, population_size)
    history = []
    cache = {}
    scored = []
    stale_generations = 0
    previous_best_key = None
    successful_evaluations = 0

    for generation in range(generations):
        if generation == 0:
            proposals = population
            evaluation_budget = len(proposals)
            next_size = len(population)
            elite_candidates = []
        else:
            next_size = max(
                minimum_population_size,
                int(len(population) * (1.0 - population_reduction_rate)),
            )
            elite_candidates = [
                (dict(result), dict(candidate))
                for result, candidate in scored[:min(elite_count, next_size)]
            ]
            diversity = _population_diversity(population, parameters)
            rate, scale = _adaptive_mutation(
                mutation_rate, mutation_scale, diversity, stale_generations
            )
            proposal_count = max(next_size * 3, next_size + 4)
            proposals = _make_proposals(
                parameters, scored, proposal_count, rate, scale,
                diversity, rng, target_keff, keff_tolerance,
                ppf_weight, keff_weight, ppf_scale,
            )
            evaluation_budget = min(
                len(proposals), max(2, ceil(next_size * surrogate_evaluation_fraction))
            )
            proposals = _screen_proposals(
                parameters, proposals, cache, evaluation_budget,
                target_keff, keff_tolerance, ppf_weight, keff_weight,
                ppf_scale, seed + generation, generation + 1,
            )
        if generation == 0:
            print(
                f"Generation {generation + 1}/{generations}: initial generation, "
                f"evaluating {len(proposals)} candidates with OpenMC."
            )

        _print_population(parameters, proposals, generation + 1, generations)
        # Carry measured parents forward without rerunning them. New measured
        # candidates compete with the previous population for survival.
        current = [] if generation == 0 else [
            (dict(result), dict(candidate)) for result, candidate in scored
        ]
        sims_this_generation = 0
        for index, candidate in enumerate(proposals):
            key = _candidate_key(candidate, parameters)
            cached = key in cache
            if cached:
                result, error, _ = cache[key]
                result = dict(result)
            else:
                try:
                    result = dict(evaluate(candidate, generation, index))
                    if "ppf" not in result:
                        raise ValueError("evaluation result is missing 'ppf'")
                    if target_keff is not None and "keff" not in result:
                        raise ValueError("evaluation result is missing 'keff'")
                    error = None
                except Exception as exception:
                    error = f"{type(exception).__name__}: {exception}"
                    result = {
                        "ppf": float("inf"),
                        "keff": float("inf") if target_keff is not None else None,
                        "evaluation_error": error,
                    }
                cache[key] = (dict(result), error, dict(candidate))
                sims_this_generation += 1
                if error is None:
                    successful_evaluations += 1
            result["objective_loss"] = _weighted_loss(
                result, target_keff, keff_tolerance,
                ppf_weight, keff_weight, ppf_scale,
            )
            cache[key] = (dict(result), error, dict(candidate))
            if error:
                print(
                    f"Generation {generation + 1}/{generations}, "
                    f"candidate {index + 1}/{len(proposals)} failed: {error}"
                )
            else:
                print(
                    f"Generation {generation + 1}/{generations}, "
                    f"{'cached result' if cached else 'simulation'} "
                    f"{index + 1}/{len(proposals)}: "
                    f"keff={_format_metric(result.get('keff'))}, "
                    f"PPF={result['ppf']:.6f}, "
                    f"loss={_format_loss(result['objective_loss'])}"
                )
            current.append((result, dict(candidate)))
        current.sort(key=lambda item: _fitness_key(
            item[0], target_keff, keff_tolerance,
            ppf_weight, keff_weight, ppf_scale,
        ))
        # Keep the best measured candidates as the next genetic population.
        elite_keys = {
            _candidate_key(candidate, parameters)
            for _, candidate in elite_candidates
        }
        survivors = elite_candidates + [
            item for item in current
            if _candidate_key(item[1], parameters) not in elite_keys
        ][:max(0, next_size - len(elite_candidates))]
        scored = sorted(survivors, key=lambda item: _fitness_key(
            item[0], target_keff, keff_tolerance,
            ppf_weight, keff_weight, ppf_scale,
        ))
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
            "simulations_this_generation": sims_this_generation,
            "unique_evaluations": len(cache),
            "population_diversity": _population_diversity(
                [candidate for _, candidate in scored], parameters
            ),
        })
        generation_best_key = _fitness_key(
            generation_result, target_keff, keff_tolerance,
            ppf_weight, keff_weight, ppf_scale,
        )
        if previous_best_key is not None and generation_best_key >= previous_best_key:
            stale_generations += 1
        else:
            stale_generations = 0
        previous_best_key = generation_best_key
        status = f"best PPF={generation_result['ppf']:.6f}"
        if target_keff is not None:
            status += f", keff={generation_result['keff']:.6f}"
            status += f" (target {target_keff:.6f} +/- {keff_tolerance:.6f})"
        status += f", loss={_format_loss(generation_result['objective_loss'])}"
        print(
            f"Generation {generation + 1}/{generations} "
            f"(population={len(population)}, OpenMC runs={sims_this_generation}, "
            f"unique evaluations={len(cache)}, "
            f"diversity={history[-1]['population_diversity']:.3f}): {status}"
        )

        ppf_reached = target_ppf is None or generation_result["ppf"] <= target_ppf
        keff_reached = _keff_is_feasible(
            generation_result, target_keff, keff_tolerance
        )
        if target_ppf is not None and ppf_reached and keff_reached:
            break
        if generation == generations - 1:
            break
        population = [dict(candidate) for _, candidate in scored]

    if successful_evaluations == 0:
        failures = [result.get("evaluation_error") for result, _ in scored]
        raise RuntimeError(
            "All genetic algorithm evaluations failed. Last errors: "
            + "; ".join(error for error in failures if error)
        )
    best_result, best_candidate = min(
        ((result, candidate) for result, _, candidate in cache.values()
         if result.get("evaluation_error") is None),
        key=lambda item: _fitness_key(
            item[0], target_keff, keff_tolerance,
            ppf_weight, keff_weight, ppf_scale,
        ),
    )
    return dict(best_candidate), dict(best_result), history


def _keff_is_feasible(result, target_keff, keff_tolerance):
    return target_keff is None or abs(result["keff"] - target_keff) <= keff_tolerance


def _candidate_key(candidate, parameters):
    return tuple(round(float(candidate[p.name]), 12) for p in parameters)


def _population_diversity(population, parameters):
    if len(population) < 2 or not parameters:
        return 0.0
    spreads = []
    for parameter in parameters:
        values = [candidate[parameter.name] for candidate in population]
        spreads.append(
            min(1.0, (max(values) - min(values)) / (parameter.upper - parameter.lower))
        )
    return sum(spreads) / len(spreads)


def _adaptive_mutation(base_rate, base_scale, diversity, stale_generations):
    stagnation_boost = min(0.5, 0.08 * stale_generations)
    low_diversity_boost = 0.25 if diversity < 0.10 else 0.0
    rate = min(0.95, base_rate + stagnation_boost + low_diversity_boost)
    scale_multiplier = 1.0 + min(2.0, 0.25 * stale_generations)
    if diversity < 0.10:
        scale_multiplier *= 1.5
    return rate, min(0.5, base_scale * scale_multiplier)


def _make_proposals(parameters, scored, count, mutation_rate,
                    mutation_scale, diversity, rng, target_keff,
                    keff_tolerance, ppf_weight, keff_weight, ppf_scale):
    def tournament():
        contestants = rng.sample(scored, k=min(3, len(scored)))
        return dict(min(
            contestants,
            key=lambda item: _fitness_key(
                item[0], target_keff, keff_tolerance,
                ppf_weight, keff_weight, ppf_scale,
            ),
        )[1])

    immigrant_fraction = 0.10 + (0.20 if diversity < 0.10 else 0.0)
    immigrant_count = max(1, int(count * immigrant_fraction))
    proposals = []
    seen = set()
    # Keep random immigrants in the mix, increasing their share as diversity falls.
    for _ in range(immigrant_count):
        candidate = {
            p.name: rng.uniform(p.lower, p.upper) for p in parameters
        }
        key = _candidate_key(candidate, parameters)
        if key not in seen:
            proposals.append(candidate)
            seen.add(key)
    while len(proposals) < count:
        first = tournament()
        second = tournament()
        child = {}
        for parameter in parameters:
            name = parameter.name
            # Blend parents rather than copying one coordinate wholesale.
            alpha = rng.random()
            value = alpha * first[name] + (1.0 - alpha) * second[name]
            if rng.random() < mutation_rate:
                value += rng.gauss(
                    0.0, mutation_scale * (parameter.upper - parameter.lower)
                )
            child[name] = min(parameter.upper, max(parameter.lower, value))
        key = _candidate_key(child, parameters)
        if key not in seen:
            proposals.append(child)
            seen.add(key)
        # Avoid an unbounded loop when the continuous space is effectively
        # collapsed by narrow bounds or a zero mutation scale.
        if len(seen) < count and len(proposals) < count and len(seen) > count * 20:
            break
    return proposals


def _screen_proposals(parameters, proposals, cache, budget, target_keff,
                      keff_tolerance, ppf_weight, keff_weight, ppf_scale, seed,
                      generation):
    available = []
    seen = set()
    for candidate in proposals:
        key = _candidate_key(candidate, parameters)
        if key not in cache and key not in seen:
            available.append(candidate)
            seen.add(key)
    if not available:
        print(
            f"Generation {generation}: no new unique proposals; all candidates "
            "were already evaluated."
        )
        return available
    if len(available) <= budget:
        print(
            f"Generation {generation}: {len(available)} unique virtual proposals "
            f"available; sending all to OpenMC (budget {budget})."
        )
        return available

    observations = [
        (result, candidate)
        for result, error, candidate in cache.values()
        if error is None
    ]
    if len(observations) < 3:
        print(
            f"Generation {generation}: only {len(observations)} successful "
            "training evaluations; selecting proposals randomly."
        )
        rng = random.Random(seed)
        rng.shuffle(available)
        return available[:budget]

    try:
        import numpy as np
        from sklearn.ensemble import RandomForestRegressor

        def normalized(candidate):
            return [
                (candidate[p.name] - p.lower) / (p.upper - p.lower)
                for p in parameters
            ]

        x = np.asarray([normalized(candidate) for _, candidate in observations])
        y = np.asarray([
            _surrogate_loss(result, target_keff, keff_tolerance,
                            ppf_weight, keff_weight, ppf_scale)
            for result, _ in observations
        ])
        candidate_x = np.asarray([normalized(candidate) for candidate in available])
        model = RandomForestRegressor(
            n_estimators=64, min_samples_leaf=2, max_features=1.0,
            random_state=seed, n_jobs=1,
        )
        model.fit(x, y)
        tree_predictions = np.asarray([
            tree.predict(candidate_x) for tree in model.estimators_
        ])
        if not np.all(np.isfinite(tree_predictions)):
            raise ValueError("RandomForestRegressor produced non-finite predictions")
        mean = tree_predictions.mean(axis=0)
        uncertainty = tree_predictions.std(axis=0)
        acquisition = mean - 0.5 * uncertainty
        ranked_indices = list(np.argsort(acquisition))
        selected = _diverse_top_candidates(
            available, candidate_x, ranked_indices, budget
        )
        selected_model = "RandomForestRegressor"
        selected_keys = {_candidate_key(candidate, parameters) for candidate in selected}
        print(
            f"Generation {generation}: using {selected_model}; ranked "
            f"{len(available)} virtual candidates and selected {len(selected)} "
            "for OpenMC."
        )
        selected_indices = [
            index for index, candidate in enumerate(available)
            if _candidate_key(candidate, parameters) in selected_keys
        ]
        for rank, index in enumerate(
            sorted(selected_indices, key=lambda i: float(acquisition[i])), start=1
        ):
            parameters_text = ", ".join(
                f"{parameter.name}={available[index][parameter.name]:.6g}"
                for parameter in parameters
            )
            print(
                f"  OpenMC candidate {rank}/{len(selected)}: "
                f"predicted loss={mean[index]:.6g}, "
                f"uncertainty={uncertainty[index]:.6g}; {parameters_text}"
            )
        return selected
    except Exception as error:
        # Surrogate fitting is advisory; preserve GA progress if it fails.
        print(
            f"Generation {generation}: RandomForestRegressor failed "
            f"({type(error).__name__}: {error}); selecting proposals randomly."
        )
        rng = random.Random(seed)
        rng.shuffle(available)
        return available[:budget]


def _surrogate_loss(result, target_keff, keff_tolerance,
                    ppf_weight, keff_weight, ppf_scale):
    ppf_score = max(0.0, result["ppf"] - 1.0) / ppf_scale
    if target_keff is None:
        return ppf_score
    deviation = abs(result["keff"] - target_keff)
    keff_score = deviation / max(keff_tolerance, 1.0e-6)
    total_weight = ppf_weight + keff_weight
    return (ppf_weight * ppf_score + keff_weight * keff_score) / total_weight


def _diverse_top_candidates(candidates, points, ranked_indices, count):
    selected = []
    selected_indices = []
    for index in ranked_indices:
        if all(
            float(((points[index] - points[chosen]) ** 2).sum()) ** 0.5 >= 0.02
            for chosen in selected_indices
        ):
            selected.append(candidates[index])
            selected_indices.append(index)
            if len(selected) == count:
                return selected
    for index in ranked_indices:
        if index not in selected_indices:
            selected.append(candidates[index])
            if len(selected) == count:
                break
    return selected


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


def _format_loss(loss):
    return "n/a (keff outside tolerance)" if loss is None else f"{loss:.6f}"


def _format_metric(value):
    return "n/a" if value is None else f"{value:.6f}"


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
