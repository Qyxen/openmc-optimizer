# OpenMC Optimizer

A small, project-style optimization driver for OpenMC models. The initial
search method is a bounded real-valued genetic algorithm. Its built-in
objective is to minimize assembly power peaking factor (PPF), optionally while
keeping k-effective within a user-set tolerance of a target. The design
variables are independent of the objective, so the same search can vary
enrichment or Gd₂O₃ rod loading using the built-in material adapter.

## What the project does

For each candidate, the driver copies the user's XML model into a separate
work directory, applies candidate values with its built-in material adapter,
adds its own `fission-q-recoverable` assembly-instance tally, runs OpenMC, and
calculates

```text
PPF = maximum assembly power / mean assembly power
```

The GA ranks candidates satisfying `target_keff` within `keff_tolerance`
(default `0.005`) by a weighted loss. The default weights are 0.7 for PPF and
0.3 for keff; candidates outside the tolerance remain infeasible regardless of
their loss. Until one is feasible, the search favors candidates closest to the
target keff. If `target_ppf` is also set, the search stops when both targets
are met. All candidate inputs, statepoints, and results are retained under the
configured work directory.

## Install

Use an environment with the OpenMC Python package and executable available.
Then install this project:

```bash
python -m pip install -e .
```

OpenMC and its nuclear data are provided by the user environment.

## Documentation

- [Python API reference](docs/api.rst)
- [Python API example](docs/example.rst)
- [Theory and methods](docs/theory.rst)

Sphinx documentation is built and deployed to GitHub Pages by
`.github/workflows/docs.yml`. In the repository's **Settings > Pages**, set
the publishing source to **GitHub Actions** to enable deployment.

## Tests

Run the unit tests locally with:

```bash
python -m pip install -e ".[test]"
python -m pytest -q
```

The ``.github/workflows/tests.yml`` workflow runs the suite on pushes and pull
requests with Python 3.10, 3.12, and 3.13. Tests cover the genetic algorithm,
input validation, and built-in material adapter without requiring OpenMC runs.

## Python API

The primary interface is Python. Put the model XML files in a directory and
write a small script with the material IDs and search bounds:

```python
from openmc_optimizer import MaterialFraction, OpenMCOptimizer

optimizer = OpenMCOptimizer(
    model_dir="model",
    variables=[
        MaterialFraction(material_id=11, component="U235", normalization_group=("U235", "U238"), bounds=(0.02, 0.05), name="A01_U235"),
        MaterialFraction(material_id=12, component="U235", normalization_group=("U235", "U238"), bounds=(0.02, 0.05), name="A02_U235"),
    ],
    target_ppf=1.35,
    target_keff=1.0,
    keff_tolerance=0.005,
    ppf_weight=0.7,
    keff_weight=0.3,
    ppf_scale=0.5,
    population_reduction_rate=0.20,
    population_size=16,
    generations=25,
)
result = optimizer.run()
print(result.best_ppf, result.best_keff, result.best_objective_loss, result.best_parameters)
```

`MaterialFraction` can vary a nuclide or element in any material. The selected
component is assigned the candidate fraction and every other composition
entry is scaled proportionally. By default the bounds refer to the whole
material. Set `normalization_group` to vary a component relative to a named
group while preserving the group's total in the material. Choose
`basis="wo"` for weight fractions or `basis="ao"` for atom fractions. For
example, vary U-235 relative to total uranium:

```python
from openmc_optimizer import MaterialFraction

MaterialFraction(
    material_id=20,
    component="U235",
    component_type="nuclide",
    normalization_group=("U235", "U238"),
    bounds=(0.015, 0.0495),
    basis="wo",
    name="U235_enrichment",
)
```

The normalization group must already be present in the material and include
the selected component. Without it, a U-235 fraction refers to the whole fuel
composition rather than the uranium isotope vector.

Use `VolumeFraction` when mixing two complete material recipes by physical
volume. The same `MaterialFraction` variable handles U-235 or any other
component; the other entries in that material scale proportionally.

The full runnable template is [examples/optimize.py](examples/optimize.py).
Users provide the model directory and a variable definition for each material
being optimized. The optimizer infers assembly root-cell IDs for common nested
lattice geometries; unusual geometries can pass `assembly_cell_ids` explicitly.
GA settings and target PPF are optional and have defaults. The API uses the
built-in material adapter and creates the PPF tally automatically.

## Prepare a model

Put the model files in a directory: `geometry.xml`, `materials.xml`,
`settings.xml`, and any other XML files required by the model. Set the
cross-section path through the environment or the user's normal OpenMC setup.
The optimizer attempts to identify assembly instances from the nested lattice
structure. If your geometry does not use an assembly-lattice hierarchy, set
`assembly_cell_ids` in the Python API. These must be cells filled by assembly
lattices, not pin cells.

Set the model directory, material IDs, bounds, and any optional GA settings in
your Python script. Replace the sample material IDs with IDs from your model.

Run the Python template with:

```bash
python examples/optimize.py
```

## Optimization variables

Each `MaterialFraction` variable carries the material ID, component, fraction
basis, and bounds. The built-in adapter supports any nuclide or element, and
`VolumeFraction` mixes two complete material recipes by physical volume. For
example, this Python variable varies a gadolinium isotope:

```python
MaterialFraction(
    material_id=11,
    component="Gd155",
    component_type="nuclide",
    basis="wo",
    bounds=(0.02, 0.05),
    name="Gd155_weight_fraction",
)
```

Material fraction bounds apply to the complete material composition. The
adapter keeps the other entries in their current proportions and scales them
to fill the remainder. Entries in the material must all use the selected
atom (`ao`) or weight (`wo`) basis. To vary the U-235 fraction, set
`component="U235"`; its bounds apply to the complete material composition.

For Gd₂O₃ volume fraction, specify the rod material and Gd₂O₃ component
material in Python:

```python
VolumeFraction(
    target_material_id=11,
    component_material_id=90,
    bounds=(0.0, 0.10),
    name="Gd2O3_volume_fraction",
)
```

The component material needs an explicit density in g/cm³ and a composition
in `materials.xml`. This variable treats the requested fraction as a physical
volume replacement: increasing Gd₂O₃ reduces the volume share of the target
fuel. Material IDs, isotope vectors, and component densities remain specific
to the user's model. The built-in variable types are defined in the Python
API, so users do not need to create an adapter file.

## PPF tally requirements

The worker loads `geometry.xml`, resolves the configured assembly root-cell
instances, creates a `CellInstanceFilter`, and scores
`fission-q-recoverable` with the tracklength estimator. Existing `tallies.xml`
content is retained and the optimizer tally is added. A run must produce the
configured statepoint file. Use a geometry with at least two assembly
instances and positive fission power in every tallied instance.

## Results

`work_dir/optimization_summary.json` contains the best parameter values,
best PPF, whether the target was reached, and per-generation history. The
optimizer also exports `optimization_history.csv` with the best and mean PPF
and best parameter values for each generation, plus
`optimization_candidates.csv` with the metrics and parameter values for each
unique candidate that required an OpenMC run. Each candidate directory
contains its input XML, OpenMC output, statepoint, and `ppf_result.json`.
Duplicate parameter vectors are cached within a run.

## Limitations in this initial implementation

- Candidate evaluations run sequentially; GA randomness is reproducible from
  the configured seed.
- The built-in PPF tally assumes assembly root-cell IDs can be supplied and
  that `fission-q-recoverable` is appropriate for the model.
- Runtime, memory, and tally precision depend on the user's OpenMC settings.
- Generic material fractions rescale all other entries proportionally; they
  do not model coupled recipe or chemistry constraints.
