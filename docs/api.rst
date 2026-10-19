Python API reference
====================

The public classes are available from the top-level ``openmc_optimizer``
package::

   from openmc_optimizer import (
       MaterialFraction,
       OpenMCOptimizer,
       OptimizationResult,
       VolumeFraction,
   )

``OpenMCOptimizer``
-------------------

Create an optimizer by passing the model directory and one or more design
variables to ``OpenMCOptimizer``. All arguments are keyword-only::

   OpenMCOptimizer(
       *,
       model_dir,
       variables,
       method="genetic",
       assembly_cell_ids=None,
       target_ppf=None,
       target_keff=None,
       keff_tolerance=0.005,
       ppf_weight=0.7,
       keff_weight=0.3,
       ppf_scale=0.5,
       population_reduction_rate=0.20,
       work_dir="results",
       population_size=20,
       generations=30,
       elite_count=2,
       mutation_rate=0.15,
       mutation_scale=0.12,
       seed=1,
       openmc_command=("openmc",),
       statepoint_file="statepoint.100.h5",
   )

Arguments
~~~~~~~~~

``method``
   Search method: ``"genetic"`` (default) or ``"bayesian"``. Bayesian
   optimization uses a Gaussian-process surrogate. Scikit-learn is installed
   as a required project dependency.

``model_dir``
   Required path to the directory containing the OpenMC XML model. It must
   contain at least ``geometry.xml``, ``materials.xml``, and ``settings.xml``.
   Relative paths are resolved from the current working directory.

``variables``
   Required non-empty list of ``MaterialFraction`` and/or ``VolumeFraction``
   instances. Each variable's ``name`` must be unique. When omitted, a
   variable's name is constructed from its component, basis, and material ID
   (or from the target and component IDs for volume mixing).

``assembly_cell_ids``
   Optional list or tuple of positive geometry cell IDs for assembly root
   cells. If omitted, the optimizer attempts to infer suitable cells from the
   nested lattice structure. Supply IDs explicitly if inference does not
   match the assemblies used by the model.

``target_ppf``
   Optional positive PPF threshold. The selected method stops early when its
   best candidate meets this PPF threshold and, if configured, the keff target
   tolerance.

``target_keff``
   Optional positive effective multiplication factor target. Candidates within
   ``keff_tolerance`` are considered feasible; among feasible candidates the
   optimizer minimizes the weighted loss. Both methods prioritize candidates
   closer to the target while no feasible candidate has been found.

``keff_tolerance``
   Absolute tolerance around ``target_keff``. Defaults to ``0.005``. A candidate
   is feasible when :math:`|k_{\mathrm{eff}}-k_{\mathrm{target}}|\le\tau_k`.
   This setting
   is ignored when ``target_keff`` is ``None``.

``ppf_weight`` and ``keff_weight``
   Non-negative weights for the PPF and keff terms in the loss. The weights
   are normalized by their sum. Defaults to ``0.7`` for PPF and ``0.3`` for
   keff. They affect ranking only after the candidate meets the keff tolerance.
   If ``target_keff`` is ``None``, only the PPF objective is active and
   ``ppf_weight`` must be positive.

``ppf_scale``
   Positive normalization scale for PPF excess above 1. Defaults to ``0.5``.
   The normalized PPF term is :math:`\max(\mathrm{PPF}-1,0)/s_{\mathrm{PPF}}`;
   the normalized keff term is :math:`|k_{\mathrm{eff}}-k_{\mathrm{target}}|/\tau_k`.

``work_dir``
   Directory for candidate runs and optimization outputs. Defaults to
   ``results`` relative to the current working directory.

``population_size``
   For the genetic method, number of candidates in the first generation. For
   Bayesian optimization, number of samples in the initial Latin-hypercube
   design. Must be at least 10.

``population_reduction_rate``
   Genetic algorithm only. Fraction by which each generation's population is
   reduced. Defaults to ``0.20`` (20%). Must be in ``[0, 1)``. Set to ``0`` to
   keep the population size constant; the minimum population size remains 10.

``generations``
   For the genetic method, maximum number of generations. For Bayesian
   optimization, maximum number of sequential surrogate-guided evaluations
   after the initial design. Must be at least 1.

``elite_count``
   Genetic algorithm only. Maximum number of top candidates copied unchanged to the next generation.
   Must be at least 1 and less than the initial ``population_size``. As the
   population shrinks, the elite count is capped to leave at least one slot
   for a new child.

``mutation_rate``
   Genetic algorithm only. Probability that each child parameter is mutated independently. Must be
   between 0 and 1 inclusive.

``mutation_scale``
   Genetic algorithm only. Standard deviation of a mutation as a fraction of that parameter's bounds
   range. Must be non-negative. A value of 0 disables mutation changes even
   when a parameter is selected for mutation.

``seed``
   Integer seed for the selected method's random choices. It does not remove
   OpenMC Monte Carlo uncertainty.

``openmc_command``
   Command used to run OpenMC in each candidate directory. Provide a sequence
   such as ``("openmc",)`` or ``("mpirun", "-n", "8", "openmc")``.

``statepoint_file``
   Name of the statepoint file OpenMC is expected to produce. Defaults to
   ``statepoint.100.h5``.

Methods
~~~~~~~

``run()``
   Executes the optimization and returns an ``OptimizationResult``. It raises
   an exception if the model directory, XML inputs, OpenMC execution, or PPF
   tally evaluation fails.

``MaterialFraction``
--------------------

Defines a variable for the fraction of one nuclide or element in a material::

   MaterialFraction(
       material_id,
       component,
       bounds,
       name=None,
       component_type="nuclide",
       basis="wo",
       normalization_group=None,
   )

``material_id``
   Positive OpenMC material ID to modify.

``component``
   Nuclide name such as ``"U235"`` or element name such as ``"O"``.

``bounds``
   Inclusive lower and upper fraction values, each within ``[0, 1]``. Without
   ``normalization_group``, the fractions describe the whole material
   composition. With a group, the values describe the selected component's
   share of that group. The applicable other entries are rescaled
   proportionally to fill the remainder.

``name``
   Optional unique parameter name used in result dictionaries and CSV column
   names. Defaults to ``"{component}_{basis}_{material_id}"``.

``component_type``
   Either ``"nuclide"`` or ``"element"``. Defaults to ``"nuclide"``.

``basis``
   Either ``"wo"`` for weight fraction or ``"ao"`` for atom fraction.
   Defaults to ``"wo"``. Every nuclide and element entry in the selected
   material must use the chosen basis in ``materials.xml``.

``normalization_group``
   Optional tuple of component names of the same ``component_type``. When
   supplied, ``bounds`` describe the selected component's fraction **within
   this group**. The total fraction of the group in the full material is
   preserved, as are all entries outside the group; the other group members
   are rescaled proportionally. Include ``component`` in the group. For
   example, ``("U235", "U238")`` makes a U-235 fraction vary relative to
   total uranium while keeping the material's total uranium amount fixed.
   Every named group component must already be present in the material.

``VolumeFraction``
------------------

Defines a physical-volume mixture between two materials::

   VolumeFraction(
       target_material_id,
       component_material_id,
       bounds,
       name=None,
   )

``target_material_id``
   Positive OpenMC material ID whose recipe and density will be modified in
   each candidate.

``component_material_id``
   Positive OpenMC material ID used as the component recipe. Its composition
   and density are read from ``materials.xml``.

``bounds``
   Inclusive lower and upper physical volume fractions in ``[0, 1]``.

``name``
   Optional unique parameter name. Defaults to
   ``"volume_fraction_{target_material_id}_{component_material_id}"``.

Both materials need explicit densities in ``g/cm3``. Increasing the volume
fraction replaces a corresponding share of the target material.

``OptimizationResult``
----------------------

``run()`` returns a frozen dataclass with these fields:

``best_parameters``
   Mapping from design variable name to the best value found.

``best_ppf``
   PPF of the best candidate found during the run.

``best_keff`` and ``best_keff_std_dev``
   Combined k-effective estimate and its standard deviation from the best
   candidate's OpenMC statepoint.

``best_objective_loss``
   Weighted normalized loss for the best feasible candidate. It is ``None``
   when a configured keff target was not met by any candidate.

``target_ppf``
   The stopping threshold passed to the optimizer, or ``None``.

``target_keff`` and ``keff_tolerance``
   The k-effective target and absolute tolerance, or ``None`` for the target.

``keff_target_reached``
   Whether the best candidate is within the requested keff tolerance, or
   ``None`` when no keff target was supplied.

``target_reached``
   Whether all configured targets were met. With only ``target_keff`` set, this
   reports whether the best candidate is within its tolerance. With neither
   target configured, it is ``None``.

``history``
   List of per-generation dictionaries. Each contains the zero-based
   ``generation`` index, best and mean PPF, best and mean keff (when available),
   keff feasibility, weighted objective loss, and ``best_parameters``.

``summary_file``
   Path to ``optimization_summary.json``.

``history_csv``
   Path to ``optimization_history.csv``. Its rows contain the generation's
   population size, best and mean PPF, best and mean keff, keff feasibility,
   best and mean objective loss, and best parameter values.

``candidates_csv``
   Path to ``optimization_candidates.csv``. It contains one row per unique
   candidate that required an OpenMC run, including its PPF, keff and keff
   uncertainty, tally summary, and parameter values. Repeated identical
   vectors are served from the in-memory cache and are not written as
   additional unique rows.

Example
-------

For a complete script with model setup, see :doc:`example`. For the behavior
of the PPF calculation and genetic algorithm, see :doc:`theory`.
