"""Constrained Bayesian optimization for expensive OpenMC evaluations."""

from __future__ import annotations

import math
import random
from collections.abc import Callable

from .config import Parameter


def optimize_bayesian(
    parameters: tuple[Parameter, ...],
    evaluate: Callable[[dict[str, float], int, int], dict],
    *,
    initial_samples: int,
    iterations: int,
    seed: int,
    target_ppf: float | None,
    target_keff: float | None,
    keff_tolerance: float,
    ppf_weight: float,
    keff_weight: float,
    ppf_scale: float,
) -> tuple[dict[str, float], dict, list[dict]]:
    """Run a Latin-hypercube start followed by constrained GP proposals.

    ``initial_samples`` is the space-filling design size. ``iterations`` is
    the maximum number of additional, sequential OpenMC evaluations.
    """
    try:
        import numpy as np
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import Matern, WhiteKernel
    except ImportError as error:
        raise RuntimeError(
            "Bayesian optimization requires scikit-learn. Install the project "
            "dependencies with `python -m pip install -e .`."
        ) from error

    rng = random.Random(seed)
    dimension = len(parameters)
    initial_points = _latin_hypercube(initial_samples, dimension, rng)
    initial_points[0] = [0.5] * dimension
    evaluations: list[tuple[dict, dict[str, float], list[float]]] = []

    print("\nInitial Bayesian design (fractions):")
    _print_design(parameters, initial_points)
    for index, point in enumerate(initial_points):
        values = _to_parameters(parameters, point)
        result = evaluate(values, 0, index)
        evaluations.append((result, values, point))
        _set_objective_loss(
            result, target_keff, keff_tolerance,
            ppf_weight, keff_weight, ppf_scale,
        )
        print(
            f"Initial sample {index + 1}/{initial_samples}: "
            f"PPF={result['ppf']:.6f}, keff={result['keff']:.6f}"
        )

    history = [_history_entry(0, initial_samples, evaluations, target_keff, keff_tolerance)]
    best_result, best_parameters = _best(evaluations, target_keff, keff_tolerance)
    if _target_met(best_result, target_ppf, target_keff, keff_tolerance):
        return best_parameters, best_result, history

    candidate_pool_size = 4096
    for iteration in range(1, iterations + 1):
        next_point = _propose(
            parameters,
            evaluations,
            rng,
            candidate_pool_size,
            target_keff,
            keff_tolerance,
            ppf_weight,
            keff_weight,
            ppf_scale,
            np,
            GaussianProcessRegressor,
            Matern,
            WhiteKernel,
            seed + iteration,
        )
        values = _to_parameters(parameters, next_point)
        result = evaluate(values, iteration, 0)
        _set_objective_loss(
            result, target_keff, keff_tolerance,
            ppf_weight, keff_weight, ppf_scale,
        )
        evaluations.append((result, values, next_point))
        best_result, best_parameters = _best(
            evaluations, target_keff, keff_tolerance
        )
        history.append(_history_entry(
            iteration, 1, evaluations, target_keff, keff_tolerance
        ))
        feasibility = (
            "feasible" if _is_feasible(best_result, target_keff, keff_tolerance)
            else "closest keff candidate"
        )
        print(
            f"Bayesian iteration {iteration}/{iterations}: "
            f"{feasibility}, PPF={best_result['ppf']:.6f}, "
            f"keff={best_result['keff']:.6f}, "
            f"loss={best_result.get('objective_loss')}"
        )
        if _target_met(best_result, target_ppf, target_keff, keff_tolerance):
            break

    return best_parameters, best_result, history


def _latin_hypercube(size: int, dimension: int, rng: random.Random) -> list[list[float]]:
    columns = []
    for _ in range(dimension):
        strata = list(range(size))
        rng.shuffle(strata)
        columns.append([(strata[row] + rng.random()) / size for row in range(size)])
    return [[columns[column][row] for column in range(dimension)] for row in range(size)]


def _to_parameters(parameters: tuple[Parameter, ...], point: list[float]) -> dict[str, float]:
    return {
        parameter.name: parameter.lower + coordinate * (parameter.upper - parameter.lower)
        for parameter, coordinate in zip(parameters, point)
    }


def _loss(result, target_keff, keff_tolerance, ppf_weight, keff_weight, ppf_scale):
    ppf_score = max(0.0, result["ppf"] - 1.0) / ppf_scale
    if target_keff is None:
        return ppf_score
    deviation = abs(result["keff"] - target_keff)
    keff_score = deviation / max(keff_tolerance, 1.0e-6)
    weight_sum = ppf_weight + keff_weight
    return (ppf_weight * ppf_score + keff_weight * keff_score) / weight_sum


def _set_objective_loss(result, target_keff, keff_tolerance,
                        ppf_weight, keff_weight, ppf_scale):
    feasible = _is_feasible(result, target_keff, keff_tolerance)
    result["objective_loss"] = (
        _loss(result, target_keff, keff_tolerance,
              ppf_weight, keff_weight, ppf_scale)
        if feasible else None
    )


def _is_feasible(result, target_keff, keff_tolerance):
    return target_keff is None or abs(result["keff"] - target_keff) <= keff_tolerance


def _best(evaluations, target_keff, keff_tolerance):
    def rank(item):
        result, values, _ = item
        if _is_feasible(result, target_keff, keff_tolerance):
            return (0, result["objective_loss"])
        return (1, abs(result["keff"] - target_keff), result["ppf"])

    result, values, _ = min(evaluations, key=rank)
    return result, dict(values)


def _target_met(result, target_ppf, target_keff, keff_tolerance):
    return (
        target_ppf is not None
        and result["ppf"] <= target_ppf
        and _is_feasible(result, target_keff, keff_tolerance)
    )


def _history_entry(generation, population_size, evaluations,
                   target_keff, keff_tolerance):
    result, values = _best(evaluations, target_keff, keff_tolerance)
    results = [item[0] for item in evaluations]
    feasible_losses = [item["objective_loss"] for item in results
                       if item["objective_loss"] is not None]
    return {
        "generation": generation,
        "population_size": population_size,
        "evaluations_so_far": len(evaluations),
        "best_ppf": result["ppf"],
        "mean_ppf": sum(item["ppf"] for item in results) / len(results),
        "best_keff": result["keff"],
        "mean_keff": sum(item["keff"] for item in results) / len(results),
        "keff_within_tolerance": _is_feasible(result, target_keff, keff_tolerance),
        "best_objective_loss": result["objective_loss"],
        "mean_objective_loss": (
            sum(feasible_losses) / len(feasible_losses) if feasible_losses else None
        ),
        "best_parameters": values,
    }


def _propose(parameters, evaluations, rng, pool_size, target_keff,
             keff_tolerance, ppf_weight, keff_weight, ppf_scale, np,
             GaussianProcessRegressor, Matern, WhiteKernel, seed):
    dimension = len(parameters)
    observed_x = np.asarray([item[2] for item in evaluations], dtype=float)
    candidate_pool = np.asarray(
        [[rng.random() for _ in range(dimension)] for _ in range(pool_size)],
        dtype=float,
    )
    distances = np.linalg.norm(
        candidate_pool[:, None, :] - observed_x[None, :, :], axis=2
    ).min(axis=1)
    candidate_pool = candidate_pool[distances > 1.0e-6]
    if not len(candidate_pool):
        raise RuntimeError("Bayesian optimizer could not generate a new candidate")

    kernel = Matern(
        length_scale=np.ones(dimension),
        length_scale_bounds=(1.0e-2, 1.0e2),
        nu=2.5,
    ) + WhiteKernel(noise_level=1.0e-4, noise_level_bounds=(1.0e-8, 1.0e-1))

    losses = np.asarray([
        _loss(item[0], target_keff, keff_tolerance,
              ppf_weight, keff_weight, ppf_scale)
        for item in evaluations
    ], dtype=float)
    objective_noise = np.asarray([
        max((keff_weight / (ppf_weight + keff_weight)
             * item[0].get("keff_std_dev", 0.0)
             / max(keff_tolerance, 1.0e-6)) ** 2, 1.0e-6)
        if target_keff is not None else 1.0e-6
        for item in evaluations
    ])
    objective_gp = _fit_gp(
        GaussianProcessRegressor, kernel, observed_x, losses,
        objective_noise, seed, np,
    )
    objective_mean, objective_std = objective_gp.predict(candidate_pool, return_std=True)

    if target_keff is None:
        incumbent = min(losses)
        acquisition = _expected_improvement(objective_mean, objective_std, incumbent, np)
    else:
        constraint_values = np.asarray([
            abs(item[0]["keff"] - target_keff) - keff_tolerance
            for item in evaluations
        ])
        constraint_noise = np.asarray([
            max(item[0].get("keff_std_dev", 0.0) ** 2, 1.0e-8)
            for item in evaluations
        ])
        constraint_gp = _fit_gp(
            GaussianProcessRegressor, kernel, observed_x, constraint_values,
            constraint_noise, seed + 1, np,
        )
        constraint_mean, constraint_std = constraint_gp.predict(
            candidate_pool, return_std=True
        )
        feasible_losses = [item[0]["objective_loss"] for item in evaluations
                           if item[0]["objective_loss"] is not None]
        if feasible_losses:
            incumbent = min(feasible_losses)
            probability_feasible = _normal_cdf(
                -constraint_mean / np.maximum(constraint_std, 1.0e-12), np
            )
            acquisition = (
                _expected_improvement(objective_mean, objective_std, incumbent, np)
                * probability_feasible
            )
        else:
            # Before a feasible point exists, explore optimistically toward
            # the feasible boundary rather than minimizing an arbitrary penalty.
            acquisition = -(constraint_mean - 2.0 * constraint_std)

    return candidate_pool[int(np.argmax(acquisition))].tolist()


def _fit_gp(GaussianProcessRegressor, kernel, x, y, noise, seed, np):
    scale = max(float(np.std(y)), 1.0e-6)
    alpha = np.square(noise / scale) + 1.0e-8
    model = GaussianProcessRegressor(
        kernel=kernel,
        alpha=alpha,
        normalize_y=True,
        n_restarts_optimizer=1,
        random_state=seed,
    )
    model.fit(x, y)
    return model


def _expected_improvement(mean, standard_deviation, incumbent, np):
    improvement = incumbent - mean
    safe_std = [max(float(value), 1.0e-12) for value in standard_deviation]
    z = [float(delta) / sigma for delta, sigma in zip(improvement, safe_std)]
    return np.asarray([
        max(0.0, float(delta) * _normal_cdf(value, np)
            + sigma * _normal_pdf(value))
        for delta, value, sigma in zip(improvement, z, safe_std)
    ])


def _normal_cdf(values, np):
    return 0.5 * (1.0 + np.vectorize(math.erf)(values / math.sqrt(2.0)))


def _normal_pdf(value):
    return math.exp(-0.5 * value * value) / math.sqrt(2.0 * math.pi)


def _print_design(parameters, points):
    columns = ["Sample", *(parameter.name for parameter in parameters)]
    rows = [
        [str(index + 1), *(
            f"{parameter.lower + point[column] * (parameter.upper - parameter.lower):.6g}"
            for column, parameter in enumerate(parameters)
        )]
        for index, point in enumerate(points)
    ]
    widths = [max(len(column), *(len(row[index]) for row in rows))
              for index, column in enumerate(columns)]
    border = "+-" + "-+-".join("-" * width for width in widths) + "-+"
    print(border)
    print("| " + " | ".join(value.ljust(width) for value, width in zip(columns, widths)) + " |")
    print(border)
    for row in rows:
        print("| " + " | ".join(value.ljust(width) for value, width in zip(row, widths)) + " |")
    print(border)
