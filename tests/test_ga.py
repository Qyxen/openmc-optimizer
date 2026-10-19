from openmc_optimizer.config import Parameter
from openmc_optimizer.ga import optimize


def _run(seed):
    evaluated = []

    def evaluate(candidate, generation, index):
        evaluated.append((generation, index, dict(candidate)))
        value = candidate["x"]
        return {"ppf": 1.0 + (value - 0.23) ** 2}

    result = optimize(
        (Parameter("x", 0.0, 1.0),),
        evaluate,
        population_size=8,
        generations=4,
        elite_count=2,
        mutation_rate=0.5,
        mutation_scale=0.1,
        seed=seed,
        target_ppf=None,
    )
    return result, evaluated


def test_ga_is_repeatable_and_stays_within_bounds():
    first, first_evaluated = _run(seed=91)
    second, second_evaluated = _run(seed=91)

    assert first == second
    assert first_evaluated == second_evaluated
    best_parameters, best_result, history = first
    assert 0.0 <= best_parameters["x"] <= 1.0
    assert best_result["ppf"] >= 1.0
    assert len(history) == 4
    assert [row["generation"] for row in history] == [0, 1, 2, 3]


def test_ga_starts_with_parameter_midpoint_and_stops_on_target():
    evaluated = []

    def evaluate(candidate, generation, index):
        evaluated.append((generation, index, dict(candidate)))
        return {"ppf": 1.0 + abs(candidate["x"] - 0.5)}

    _, _, history = optimize(
        (Parameter("x", 0.0, 1.0),),
        evaluate,
        population_size=5,
        generations=10,
        elite_count=1,
        mutation_rate=0.0,
        mutation_scale=0.0,
        seed=4,
        target_ppf=1.0,
    )

    assert evaluated[0] == (0, 0, {"x": 0.5})
    assert len(evaluated) == 5  # It scores the complete generation before stopping.
    assert len(history) == 1
    assert history[0]["best_ppf"] == 1.0
