Theory and methods
==================

This page describes how the optimizer represents candidate designs, estimates
assembly power peaking factor (PPF), and searches the design space with a
genetic algorithm (GA).

Optimization problem
--------------------

Each candidate design is a vector of bounded real-valued parameters
``x = (x_1, x_2, ..., x_n)``. The OpenMC model is modified using those values,
then run to estimate assembly powers and k-effective. The objective is to
minimize a normalized weighted loss over candidates meeting the optional keff
constraint::

   minimize_x loss(PPF(x), keff(x))

subject to the lower and upper bounds supplied for every parameter. An
optional ``target_keff`` defines a feasibility constraint: a candidate is
feasible when its calculated keff lies within ``keff_tolerance`` of the
target. The GA ranks feasible candidates by the weighted loss below. Until it
finds one, it favors candidates with the smallest absolute keff deviation. The
default absolute tolerance is 0.005. An optional ``target_ppf`` is an early
stopping threshold; when both targets are supplied, both must be met before
early stopping.

Within the feasible set, candidates are ranked by normalized weighted loss::

   PPF_score = max(PPF - 1, 0) / ppf_scale
   keff_score = abs(keff - target_keff) / keff_tolerance
   loss = (ppf_weight * PPF_score + keff_weight * keff_score)
          / (ppf_weight + keff_weight)

The default weights are 0.7 for PPF and 0.3 for keff, and the default
``ppf_scale`` is 0.5. The scale makes PPF excess dimensionless and comparable
to keff deviation measured in tolerance units. Changing the weights changes
preferences among feasible candidates; it cannot make an out-of-tolerance
candidate feasible. Before a feasible candidate is found, the GA prioritizes
the smallest keff deviation.

Assembly power peaking factor
-----------------------------

For each assembly instance ``i``, the optimizer obtains an integrated
``fission-q-recoverable`` tally value ``P_i``. It calculates::

   mean_power = (1 / N) * sum(P_i for i = 1 ... N)
   PPF = max(P_i) / mean_power

Here ``N`` is the number of tallied assembly instances. Equal assembly
instances are weighted equally in the mean. Since every instance uses the
same tally score and normalization, this ratio is independent of the overall
power normalization. A uniform distribution has PPF 1; larger values indicate
more variation in assembly power.

The worker creates a ``CellInstanceFilter`` and scores
``fission-q-recoverable`` with the tracklength estimator. By default, it
searches nested lattices for assembly root cells. This inference is intended
for common core/assembly lattice hierarchies. If the inferred cells are not
the intended assemblies, provide ``assembly_cell_ids`` explicitly. Those IDs
should identify geometry cells filled by assembly lattices. The model should
have at least two tallied assembly instances and positive tally values in all
of them.

Material design variables
-------------------------

``MaterialFraction`` changes a selected nuclide or element fraction in a
material. Without a ``normalization_group``, the supplied value ``f`` is a
fraction of the **whole material composition**, measured on the configured
atom (``ao``) or weight (``wo``) basis. If the other entries have fractions
``x_i``, the updated composition is::

   selected component: f
   each other entry i: (1 - f) * x_i / sum(x_j for all other entries j)

Thus, the non-selected entries retain their relative proportions and share the
remaining fraction. If a component is not already listed, the adapter adds it
to the material and scales the existing entries to the remainder.

An optional ``normalization_group`` changes the selected fraction relative to
a named subset of nuclides or elements. The group must already be present in
the material and include the selected component. The group's total material
fraction and all entries outside it remain unchanged; the group's other
members scale proportionally to fill the remainder. For uranium enrichment,
``component="U235"`` and ``normalization_group=("U235", "U238")`` make the
candidate value the U-235 weight fraction of uranium while preserving total
uranium mass. Without this group, the U-235 fraction refers to the whole fuel
material instead. This mechanism supports any named group, but does not model
other coupled rules such as oxide stoichiometry unless they are represented by
the group.

``VolumeFraction`` mixes a target material with a component material using a
physical volume fraction ``v``. Given their initial densities ``rho_t`` and
``rho_c``, the mixture density is::

   rho_mix = (1 - v) * rho_t + v * rho_c

The mass contribution of each nuclide or element is based on its material's
volume share, density, and weight fraction. The adapter combines these mass
contributions and normalizes the resulting composition by total mass. Both
source materials therefore need valid compositions and explicit densities in
``g/cm3``. The target material's density and composition are updated in the
candidate model; the source component material is retained as the recipe
definition.

Genetic algorithm
-----------------

The genetic algorithm is a bounded, real-valued search. A candidate, also
called an individual, is a vector with one real number for each parameter::

   x = (x_1, x_2, ..., x_n)
   lower_j <= x_j <= upper_j

For example, if the design varies U-235 fractions in two materials, a
candidate contains one fraction for the first material and one for the
second. The candidate is applied to a fresh copy of the model XML, and its
fitness is calculated from the resulting OpenMC run. Out-of-tolerance keff
candidates rank behind feasible candidates; feasible candidates are compared
using the configured weighted loss.

Initial population
~~~~~~~~~~~~~~~~~~

The user sets ``population_size`` for the first generation; it must be at
least 10. The algorithm creates that many candidates. Every parameter in each
candidate is drawn independently from a uniform distribution over that
parameter's ``(lower, upper)`` interval. It then replaces the first candidate
with the vector of interval midpoints::

   x_j = (lower_j + upper_j) / 2

This ensures that the center point is included in every run. It does not imply
that the center is a good design; it is simply one reproducible starting
sample. The other initial values depend on ``seed``.

Fitness evaluation and ranking
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each candidate in a generation is evaluated, then candidates are sorted by
feasibility and weighted loss. The generation's population size, best and
mean PPF and keff, and objective loss are recorded in the history. The
optimizer tracks the best candidate found over the full run separately from
the current generation's best.

OpenMC evaluations are performed sequentially. Their cost usually dominates
the GA's arithmetic, so ``population_size`` controls how many candidates are
considered in the first generation. Repeated identical parameter
vectors are retrieved from an in-memory cache and do not trigger another
OpenMC run. This cache is exact-value based; near-equal candidates are still
evaluated separately.

Elitism and parent selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After a generation is evaluated, the next population is reduced by the
user-configured ``population_reduction_rate``, rounded down, with a minimum
size of 10. Elitism copies up to
``elite_count`` of the top candidates unchanged into that smaller population;
the count is capped to leave at least one slot for a new child. These
candidates do not undergo crossover or mutation. This keeps strong solutions
available for later generations and carries the best generation member
forward.

The remaining positions are filled with children. To choose each parent, the
algorithm samples up to three members of the current generation, without
replacement within that tournament, and selects the member ranked highest by
the fitness criterion. It runs this selection twice for each child. The two selections are
independent, so the same individual can be selected as both parents. A
three-member tournament gives good candidates a higher chance of reproduction
while still allowing weaker candidates to be selected. Candidate quality uses
the same feasibility and weighted-loss ranking as generation sorting.

Crossover and mutation
~~~~~~~~~~~~~~~~~~~~~~

The algorithm uses uniform, per-parameter crossover. For every parameter in
a child, it copies the value from the first or second parent with equal
probability. The choice is made independently for each parameter, so a child
can combine values from both parents. There is no separate crossover-rate
setting: this choice is made for every child gene.

After crossover, each parameter is mutated independently with probability
``mutation_rate``. When selected, a zero-mean Gaussian perturbation is added::

   sigma_j = mutation_scale * (upper_j - lower_j)
   x'_j = x_j + Normal(mean=0, standard_deviation=sigma_j)

The resulting value is clipped to the parameter's bounds. Thus,
``mutation_scale`` is relative to each parameter's own range, even when the
parameters use different units or interval widths. A larger mutation rate
makes more child values change; a larger mutation scale makes those changes
larger on average. Excessively small mutation can limit exploration, while
large mutation can disrupt useful parameter combinations. These controls are
heuristic and may need tuning for the model and number of design variables.

Next generation and termination
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each next generation contains the unchanged elite candidates plus newly
created children until it reaches the previous generation's size multiplied
by ``1 - population_reduction_rate``, rounded down to a minimum of 10. This
process repeats for at most
``generations`` generations. After evaluating a full generation,
the algorithm stops early if that generation's best candidate meets
``target_ppf`` and, when configured, the keff tolerance. Otherwise, it stops
after the requested final generation. The history uses zero-based generation
numbers, while the progress message prints generation counts starting at one.

``seed`` controls the GA's random sampling, parent selection, crossover, and
mutation, so it makes those choices reproducible. It does not remove OpenMC
Monte Carlo uncertainty. Different candidate runs can have noisy PPF
estimates, which can change their ranking and the subsequent search path. Set
adequate particles, batches, and inactive batches in the user's OpenMC model
for the desired precision. This GA does not estimate uncertainty, average
repeated candidate runs, or guarantee a global optimum.

The generation CSV records best and mean PPF and the best parameter vector
for each generation. The candidate CSV records unique candidate vectors that
required an OpenMC run, along with their PPF and tally summary values.

Further details and a runnable script are in :doc:`example`.
