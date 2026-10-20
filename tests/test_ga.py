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


def test_ga_prints_metrics_for_each_simulation(capsys):
    def evaluate(candidate, generation, index):
        return {"ppf": 1.2, "keff": 1.0}

    optimize(
        (Parameter("x", 0.0, 1.0),),
        evaluate,
        population_size=2,
        generations=1,
        elite_count=1,
        mutation_rate=0.0,
        mutation_scale=0.0,
        seed=1,
        target_ppf=None,
        target_keff=1.0,
        keff_tolerance=0.01,
    )

    output = capsys.readouterr().out
    assert output.count("simulation ") == 2
    assert output.count("keff=1.000000, PPF=1.200000, loss=") == 2


def test_ga_surrogate_reduces_later_generation_evaluations():
    evaluated = []

    def evaluate(candidate, generation, index):
        evaluated.append((generation, index))
        x = candidate["x"]
        return {"ppf": 1.0 + (x - 0.31) ** 2}

    _, _, history = optimize(
        (Parameter("x", 0.0, 1.0),),
        evaluate,
        population_size=10,
        generations=3,
        elite_count=2,
        mutation_rate=0.2,
        mutation_scale=0.1,
        seed=29,
        target_ppf=None,
        population_reduction_rate=0.0,
        surrogate_evaluation_fraction=0.5,
    )

    assert history[0]["simulations_this_generation"] == 10
    assert all(row["simulations_this_generation"] <= 5 for row in history[1:])
    assert len(evaluated) < 10 * 3


def test_ga_continues_after_an_individual_evaluation_failure():
    def evaluate(candidate, generation, index):
        if index == 0:
            raise RuntimeError("simulated failure")
        return {"ppf": 1.2}

    best_parameters, best_result, history = optimize(
        (Parameter("x", 0.0, 1.0),),
        evaluate,
        population_size=2,
        generations=1,
        elite_count=1,
        mutation_rate=0.0,
        mutation_scale=0.0,
        seed=1,
        target_ppf=None,
    )

    assert 0.0 <= best_parameters["x"] <= 1.0
    assert best_result["ppf"] == 1.2
    assert history[0]["simulations_this_generation"] == 2


def test_ga_reports_random_forest_surrogate(capsys):
    def evaluate(candidate, generation, index):
        return {"ppf": 1.0 + (candidate["x"] - 0.4) ** 2}

    optimize(
        (Parameter("x", 0.0, 1.0),),
        evaluate,
        population_size=10,
        generations=2,
        elite_count=2,
        mutation_rate=0.2,
        mutation_scale=0.1,
        seed=9,
        target_ppf=None,
        surrogate_evaluation_fraction=0.5,
    )

    output = capsys.readouterr().out
    assert "using RandomForestRegressor" in output
