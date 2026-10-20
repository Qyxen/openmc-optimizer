from openmc_optimizer.bayesian import (
    _best,
    _is_feasible,
    _latin_hypercube,
    _local_candidate_pool,
    optimize_bayesian,
)
from openmc_optimizer.config import Parameter


def test_latin_hypercube_uses_each_stratum_once_per_parameter():
    import random

    design = _latin_hypercube(size=12, dimension=3, rng=random.Random(7))

    assert len(design) == 12
    for column in range(3):
        strata = {int(point[column] * 12) for point in design}
        assert strata == set(range(12))
        assert all(0.0 <= point[column] < 1.0 for point in design)


def test_local_candidate_pool_is_bounded_and_repeatable():
    import numpy as np

    evaluations = [
        ({"ppf": 1.1, "keff": 1.0, "objective_loss": 0.2}, {}, [0.5] * 12),
    ]

    first = _local_candidate_pool(
        evaluations, 128, 12, 1.0, 0.005, np, seed=17
    )
    second = _local_candidate_pool(
        evaluations, 128, 12, 1.0, 0.005, np, seed=17
    )

    assert first.shape == (128, 12)
    assert np.all((first >= 0.0) & (first <= 1.0))
    assert np.array_equal(first, second)
    assert np.any(np.linalg.norm(first - 0.5, axis=1) < 0.5)


def test_best_prefers_feasible_candidate_over_lower_ppf_infeasible_candidate():
    evaluations = [
        ({"ppf": 1.01, "keff": 1.02, "objective_loss": None}, {"x": 0.9}, [0.9]),
        ({"ppf": 1.20, "keff": 1.00, "objective_loss": 0.4}, {"x": 0.5}, [0.5]),
    ]

    result, parameters = _best(evaluations, target_keff=1.0, keff_tolerance=0.005)

    assert result["keff"] == 1.0
    assert parameters == {"x": 0.5}
    assert _is_feasible(result, target_keff=1.0, keff_tolerance=0.005)


def test_bayesian_search_returns_bounded_candidate_and_history():
    evaluated = []

    def evaluate(candidate, generation, index):
        evaluated.append((generation, index, dict(candidate)))
        x = candidate["x"]
        return {
            "ppf": 1.0 + (x - 0.3) ** 2,
            "keff": 1.0,
            "keff_std_dev": 0.001,
        }

    best_parameters, best_result, history = optimize_bayesian(
        (Parameter("x", 0.0, 1.0),),
        evaluate,
        initial_samples=10,
        iterations=2,
        seed=31,
        target_ppf=None,
        target_keff=1.0,
        keff_tolerance=0.005,
        ppf_weight=0.7,
        keff_weight=0.3,
        ppf_scale=0.5,
    )

    assert len(evaluated) == 12
    assert 0.0 <= best_parameters["x"] <= 1.0
    assert best_result["keff"] == 1.0
    assert best_result["objective_loss"] is not None
    assert [row["generation"] for row in history] == [0, 1, 2]
    assert history[-1]["evaluations_so_far"] == 12


def test_bayesian_prints_loss_for_each_initial_sample(capsys):
    def evaluate(candidate, generation, index):
        return {"ppf": 1.2, "keff": 1.0, "keff_std_dev": 0.001}

    optimize_bayesian(
        (Parameter("x", 0.0, 1.0),),
        evaluate,
        initial_samples=2,
        iterations=0,
        seed=3,
        target_ppf=None,
        target_keff=1.0,
        keff_tolerance=0.01,
        ppf_weight=0.7,
        keff_weight=0.3,
        ppf_scale=0.5,
    )

    output = capsys.readouterr().out
    assert output.count("loss=") == 2
    assert output.count("keff=1.000000, PPF=1.200000, loss=") == 2
