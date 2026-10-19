Examples
========

This page shows a complete setup for optimizing assembly power peaking factor
(PPF) and k-effective with the OpenMC Optimizer Python API.
The API edits copies of the model XML files for each candidate, adds its PPF
tally, runs OpenMC, and reads PPF and keff from each statepoint.

Requirements
------------

Install the optimizer in the Python environment that has the OpenMC Python
package and executable available::

   python -m pip install -e .

Place the model XML files in one directory. At minimum, the model needs
``geometry.xml``, ``materials.xml``, and ``settings.xml``. Include any other
XML files needed by the model, such as ``tallies.xml``. Make sure the OpenMC
cross-section data path is configured in the environment.

Define and run an optimization
------------------------------

Create a Python script such as ``optimize.py``. Replace the sample material
IDs and bounds with values from your model::

   from pathlib import Path

   from openmc_optimizer import MaterialFraction, OpenMCOptimizer

   project_dir = Path(__file__).resolve().parent

   optimizer = OpenMCOptimizer(
       model_dir=project_dir / "model",
       variables=[
           MaterialFraction(
               material_id=11,
               component="U235",
               component_type="nuclide",
               basis="wo",
               normalization_group=("U235", "U238"),
               bounds=(0.02, 0.05),
               name="assembly_A_U235",
           ),
           MaterialFraction(
               material_id=12,
               component="U235",
               component_type="nuclide",
               basis="wo",
               normalization_group=("U235", "U238"),
               bounds=(0.02, 0.05),
               name="assembly_B_U235",
           ),
       ],
       target_ppf=1.35,
       target_keff=1.0,
       keff_tolerance=0.005,
       ppf_weight=0.7,
       keff_weight=0.3,
       ppf_scale=0.5,
       population_reduction_rate=0.20,
       work_dir=project_dir / "results",
       population_size=16,
       generations=25,
       elite_count=2,
       mutation_rate=0.20,
       mutation_scale=0.10,
       seed=14,
   )

   result = optimizer.run()
   print(f"Best PPF: {result.best_ppf:.5f}")
   print(f"Best keff: {result.best_keff:.5f} +/- {result.best_keff_std_dev:.5f}")
   print(f"Best weighted loss: {result.best_objective_loss:.5f}")
   print(f"Best parameter values: {result.best_parameters}")
   print(f"Target reached: {result.target_reached}")
   print(f"Summary file: {result.summary_file}")
   print(f"Generation history CSV: {result.history_csv}")
   print(f"Unique candidates CSV: {result.candidates_csv}")

Run the script from the environment configured for OpenMC::

   python optimize.py

Use Bayesian optimization
-------------------------

Select the surrogate-guided method in the same API call. Scikit-learn is
installed with the standard package install. Here ``population_size`` sets the
number of initial space-filling samples, and ``generations`` sets the maximum
number of additional sequential evaluations::

   optimizer = OpenMCOptimizer(
       model_dir=project_dir / "model",
       variables=[...],
       method="bayesian",
       population_size=16,
       generations=25,
       target_keff=1.0,
       keff_tolerance=0.005,
       seed=14,
   )

Bayesian optimization evaluates candidates sequentially and proposes each new
point using Gaussian-process predictions of the objective and keff constraint.
It can be useful when each OpenMC run is costly and the number of design
variables is modest. The runtime still depends on the configured number of
OpenMC evaluations and the cost of fitting the surrogate.

The optimizer infers assembly root-cell IDs for common nested lattice
geometries. If inference is unsuitable for a model, provide explicit IDs with
``assembly_cell_ids``::

   optimizer = OpenMCOptimizer(
       model_dir="model",
       variables=[...],
       assembly_cell_ids=[101, 102, 103, 104],
   )

Design variables
----------------

``MaterialFraction`` varies the fraction of one nuclide or element in a
material. Without ``normalization_group``, bounds refer to the **whole
material composition**. When a group is provided, bounds refer to the
selected component's share within that group; the group's total fraction in
the material remains fixed. In either case, the other applicable entries are
scaled proportionally to fill the remainder. Use ``basis="wo"`` for weight
fractions or ``basis="ao"`` for atom fractions. All composition entries in
the selected material must use the chosen basis.

For example, to vary the weight fraction of gadolinium-155 in material 90::

   MaterialFraction(
       material_id=90,
       component="Gd155",
       component_type="nuclide",
       basis="wo",
       normalization_group=("Gd155", "Gd157"),
       bounds=(0.0, 0.10),
       name="Gd155_fraction",
   )

``VolumeFraction`` mixes two complete material recipes by physical volume.
For example, to vary the volume fraction of a Gd2O3 recipe mixed into a fuel
material::

   from openmc_optimizer import VolumeFraction

   VolumeFraction(
       target_material_id=20,
       component_material_id=90,
       bounds=(0.0, 0.10),
       name="Gd2O3_loading",
   )

For volume mixing, both materials need explicit densities in ``g/cm3`` and
composition definitions in ``materials.xml``.

PPF and output files
--------------------

The optimizer uses a weighted loss combining normalized PPF and keff terms
for candidates that satisfy the configured keff tolerance. Set
``ppf_weight``, ``keff_weight``, and ``ppf_scale`` to control the tradeoff.
The PPF definition is:

.. math::

   \mathrm{PPF}=\frac{\max_i P_i}{\frac{1}{N}\sum_{i=1}^{N}P_i}.

The optimizer adds a ``fission-q-recoverable`` tally over assembly instances.
Use at least two assembly instances and ensure each tallied assembly has
positive fission power. Existing ``tallies.xml`` content is retained.

The work directory contains candidate model files and OpenMC outputs. The
``optimization_summary.json`` file records the best parameter values, best
PPF, keff, weighted loss, target status, and per-generation history. Each evaluated candidate also
gets a ``ppf_result.json`` result file. The optimizer exports
``optimization_history.csv`` with each generation's population size, best and
mean metrics, and the generation's best parameter values. It also writes
``optimization_candidates.csv`` with the metrics and parameter values for
each unique candidate that required an OpenMC run. CSV file paths are
available as ``result.history_csv`` and ``result.candidates_csv``.

The runnable project template is available at
``examples/optimize.py`` from the optimizer repository.
