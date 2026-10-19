Theory and methods
==================

This page describes how the optimizer represents candidate designs, estimates
assembly power peaking factor (PPF), and searches the design space with a
genetic algorithm (GA) or constrained Bayesian optimization (BO).

Optimization problem
--------------------

Each candidate design is a vector of bounded real-valued parameters
:math:`\mathbf{x}=(x_1,x_2,\ldots,x_n)`. The OpenMC model is modified using those values,
then run to estimate assembly powers and k-effective. The objective is to
minimize a normalized weighted loss over candidates meeting the optional keff
constraint:

.. math::

   \underset{\mathbf{x}}{\operatorname{minimize}}\;
   \mathcal{L}\bigl(\operatorname{PPF}(\mathbf{x}),k_{\mathrm{eff}}(\mathbf{x})\bigr)

subject to the lower and upper bounds supplied for every parameter. An
optional ``target_keff`` defines a feasibility constraint: a candidate is
feasible when its calculated keff lies within ``keff_tolerance`` of the
target. The GA ranks feasible candidates by the weighted loss below. Until it
finds one, it favors candidates with the smallest absolute keff deviation. The
default absolute tolerance is 0.005. An optional ``target_ppf`` is an early
stopping threshold; when both targets are supplied, both must be met before
early stopping.

Within the feasible set, candidates are ranked by normalized weighted loss:

.. math::

   \begin{aligned}
   S_{\mathrm{PPF}} &= \frac{\max(\mathrm{PPF}-1,\,0)}{s_{\mathrm{PPF}}}, \\
   S_k &= \frac{|k_{\mathrm{eff}}-k_{\mathrm{target}}|}{\tau_k}, \\
   \mathcal{L} &= \frac{w_{\mathrm{PPF}}S_{\mathrm{PPF}}+w_kS_k}
                        {w_{\mathrm{PPF}}+w_k}.
   \end{aligned}

The default weights are 0.7 for PPF and 0.3 for keff, and the default
``ppf_scale`` is 0.5. The scale makes PPF excess dimensionless and comparable
to keff deviation measured in tolerance units. Changing the weights changes
preferences among feasible candidates; it cannot make an out-of-tolerance
candidate feasible. Before a feasible candidate is found, the GA prioritizes
the smallest keff deviation.

Assembly power peaking factor
-----------------------------

For each power region ``i``, the optimizer obtains an integrated tally value
:math:`P_i` proportional to regional power. By default, these regions are
assembly instances tallied with ``fission-q-recoverable``. With a user tally,
they are the bins of its selected spatial filter. It calculates:

.. math::

   \overline{P}=\frac{1}{N}\sum_{i=1}^{N}P_i,
   \qquad
   \mathrm{PPF}=\frac{\max_i P_i}{\overline{P}}.

Here ``N`` is the number of regions returned by the tally. Each region is
weighted equally in the mean. Since the same tally score and normalization
are used for all bins, the ratio is independent of the overall power
normalization. A uniform distribution has PPF 1; larger values indicate more
variation in regional power.

By default, the worker creates a ``CellInstanceFilter`` and scores
``fission-q-recoverable`` with the tracklength estimator. The user may instead
name a compatible tally from ``tallies.xml`` with ``ppf_tally_name`` and
select its score with ``ppf_tally_score``. A custom tally must return one
integrated, power-proportional value per spatial region, with one region
filter and the total nuclide bin. The optimizer computes the maximum divided
by the arithmetic mean across those values; it cannot infer region boundaries
or combine extra filter, energy, or nuclide bins on the user's behalf. By
default, it searches nested lattices for assembly root cells. This inference is intended
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
:math:`x_i`, the updated composition is:

.. math::

   f_{\mathrm{selected}}'=f,
   \qquad
   f_i'=(1-f)\frac{x_i}{\sum_{j\ne\mathrm{selected}}x_j}.

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

``VolumeFraction`` modifies the target material to become the output mixture
used by the cell, adding a physical volume fraction :math:`v` of the component
material's recipe. For UO2--Gd2O3, the target ID is the UO2 fuel material
assigned to the cell (or a separate output mixture material assigned there),
and the component ID is the Gd2O3 source recipe. The geometry continues to
refer to the target material ID; the component ID is an input recipe and does
not need to be assigned to the cell. Given their initial densities
:math:`\rho_t` and :math:`\rho_c`, the resulting mixture density is:

.. math::

   \rho_{\mathrm{mix}}=(1-v)\rho_t+v\rho_c.

The mass contribution of each nuclide or element is based on its material's
volume share, density, and weight fraction. The adapter combines these mass
contributions and normalizes the resulting composition by total mass. Both
source materials therefore need valid compositions and explicit densities in
``g/cm3``. The target material's density and composition are updated in the
candidate model, so any cell that uses that target ID sees the new mixture.
The source component material is retained as the recipe definition. If the
original target recipe must also remain unchanged, use a distinct output
mixture material and make the cell reference it before optimization.

Genetic algorithm
-----------------

The genetic algorithm is a bounded, real-valued search. A candidate, also
called an individual, is a vector with one real number for each parameter:

.. math::

   \mathbf{x}=(x_1,x_2,\ldots,x_n),
   \qquad
   \ell_j\le x_j\le u_j.

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
with the vector of interval midpoints:

.. math::

   x_j=\frac{\ell_j+u_j}{2}.

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
``mutation_rate``. When selected, a zero-mean Gaussian perturbation is added:

.. math::

   \sigma_j=m(u_j-\ell_j),
   \qquad
   x_j'=x_j+\mathcal{N}(0,\sigma_j^2),

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

Bayesian optimization
---------------------

Bayesian optimization is available with ``method="bayesian"``; scikit-learn
is installed as a required project dependency. It is a sequential method: it first
evaluates a space-filling Latin-hypercube design, then fits Gaussian-process
(GP) surrogate models to the observed objective and keff constraint. Each
later OpenMC run is selected from a random candidate pool using the GP
predictions and their uncertainty. Thus, unlike the GA, it spends each new
evaluation where the current surrogate estimates a useful trade-off between
improvement and uncertainty.

When keff is not constrained, the method uses expected improvement (EI) over
the best observed objective. For a predicted objective mean :math:`\mu`,
standard deviation :math:`\sigma`, and best feasible observed loss
:math:`f_{\min}`, minimization EI is:

.. math::

   \operatorname{EI}(\mathbf{x})=(f_{\min}-\mu)\Phi(z)+\sigma\phi(z),
   \qquad z=\frac{f_{\min}-\mu}{\sigma},

where :math:`\Phi` and :math:`\phi` are the standard normal cumulative
distribution and density. With a keff target, a second GP predicts the signed
feasibility measure

.. math::

   g(\mathbf{x})=|k_{\mathrm{eff}}(\mathbf{x})-k_{\mathrm{target}}|-\tau_k.

Values at or below zero satisfy the tolerance. Once a feasible observation
exists, the acquisition score is EI multiplied by the GP-estimated
probability that :math:`g(\mathbf{x})\le 0`. Before feasibility is observed,
the optimizer chooses points with an optimistic predicted constraint value to
seek the feasible region. The keff tally standard deviation is used as an
approximation to constraint observation noise; PPF uncertainty is not
currently modeled.

``population_size`` controls the initial design size, and ``generations`` is
the maximum number of additional single-candidate evaluations. The optimizer
stops early if the best observed candidate reaches ``target_ppf`` and the keff
tolerance, when configured. GP models are fit on normalized parameter
coordinates and use a Matérn covariance with a white-noise term. The candidate
pool currently contains 4096 random points per proposal. BO is often a useful
choice when OpenMC evaluations are expensive and the parameter dimension is
modest, but model quality can degrade in high-dimensional spaces or with very
noisy tallies. The selected candidate is still the best measured OpenMC
result, not the GP's predicted optimum.

Surrogate model
~~~~~~~~~~~~~~~

Bayesian optimization is useful when evaluating the real model is expensive.
It does not use the GP to replace OpenMC or to report an unevaluated design as
the answer. Instead, it alternates between running OpenMC at selected points
and updating a surrogate that predicts how PPF, keff, and the objective loss
vary over the design space. The surrogate helps decide which point to run
next.

Before fitting a GP, each variable is mapped to the unit interval. For
parameter :math:`j`, with bounds :math:`[\ell_j,u_j]`, the normalized value is

.. math::

   z_j=\frac{x_j-\ell_j}{u_j-\ell_j}, \qquad 0\le z_j\le 1.

This makes distances comparable when parameters have different units or
bound widths. The initial design uses Latin-hypercube sampling: each
parameter's interval is divided into :math:`n` equal strata, one point is
placed in every stratum, and the strata are independently shuffled between
parameters. This spreads the initial points across each range. The first
point is replaced by the midpoint of all parameter bounds. ``population_size``
sets :math:`n`, the number of initial simulations.

After the initial runs, the implementation fits a Gaussian process (GP) using
a Matérn covariance kernel with :math:`\nu=2.5` plus a white-noise kernel.
The Matérn kernel describes how quickly nearby designs are expected to have
similar responses; its length scales are estimated from the observations.
The white-noise term allows for variation that should not be fit as a smooth
trend. A GP prediction has a mean :math:`\mu(\mathbf{x})` and standard
deviation :math:`\sigma(\mathbf{x})`. The mean is the predicted response;
the standard deviation represents the model's uncertainty at that point.
Predictions are typically more informed close to observed designs and less
certain in unexplored regions.

The GP is refit after every new OpenMC result. To select each next point, the
optimizer scores 4096 random points in normalized parameter space and
chooses the point with the highest acquisition score. Points extremely close
to previous samples are removed. This finite candidate pool makes the
acquisition search practical, but a narrow promising region can be missed if
the random pool does not sample it.

Expected improvement
~~~~~~~~~~~~~~~~~~~~

Without a keff constraint, the acquisition function is expected improvement
(EI) for minimization. Let :math:`f_{\min}` be the best measured objective
loss, and let :math:`\mu(\mathbf{x})` and :math:`\sigma(\mathbf{x})` be the
GP mean and standard deviation of loss. Define

.. math::

   z=\frac{f_{\min}-\mu(\mathbf{x})}{\sigma(\mathbf{x})}.

Then EI is

.. math::

   \operatorname{EI}(\mathbf{x})=
   (f_{\min}-\mu(\mathbf{x}))\Phi(z)+\sigma(\mathbf{x})\phi(z),

where :math:`\Phi` is the standard normal cumulative distribution function
and :math:`\phi` is its probability density function. The first term rewards
points whose predicted mean is below the best measured loss. The second term
rewards uncertainty when it leaves a chance of a large improvement. Thus, EI
balances exploitation of promising predictions with exploration of
uncertain regions. A small positive floor is applied to predicted standard
deviations for numerical stability.

Keff feasibility and weighted loss
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With ``target_keff``, the optimizer also fits a constraint GP to

.. math::

   g(\mathbf{x})=
   |k_{\mathrm{eff}}(\mathbf{x})-k_{\mathrm{target}}|-\tau_k.

A design is feasible when :math:`g(\mathbf{x})\le 0`, meaning its measured
keff is within ``keff_tolerance``. Before any feasible result exists, there
is no feasible incumbent for EI. During this phase, the optimizer chooses the
candidate with the most optimistic lower confidence bound,
:math:`\mu_g-2\sigma_g`, to search for the feasible region. The uncertainty
term means a candidate can be worth evaluating even when its predicted keff
deviation is not yet within tolerance.

After a feasible result exists, EI is based on the best feasible measured
loss. It is multiplied by the GP probability that a proposed point satisfies
the keff constraint:

.. math::

   P(g(\mathbf{x})\le 0)=
   \Phi\!\left(-\frac{\mu_g(\mathbf{x})}{\sigma_g(\mathbf{x})}\right),
   \qquad
   a(\mathbf{x})=\operatorname{EI}(\mathbf{x})P(g(\mathbf{x})\le 0).

This favors points that are promising for the objective and likely to meet
the constraint. It is a proposal score; the final ranking uses actual OpenMC
results, requiring feasibility before comparing objective loss. The
objective GP is trained on the scalar weighted-loss formula, including the
keff contribution for infeasible observations, while the separate constraint
GP guides feasibility. An infeasible result cannot be returned as the best
candidate just because its modeled loss is low.

Monte Carlo noise and practical limits
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The statepoint's keff standard deviation is used as an approximate noise
estimate for the constraint GP. The objective GP's noise estimate includes
the weighted keff uncertainty contribution when a keff target is active.
PPF tally uncertainty is not currently modeled. Consequently, the surrogate
may treat noisy PPF values as more precise than they are. More particles or
batches can improve tally precision, but increase the cost of every OpenMC
evaluation. Repeated runs at the same design are not currently used to
estimate PPF uncertainty.

``generations`` is the maximum number of sequential proposals after the
initial design, so the run can evaluate up to
``population_size + generations`` candidates. It stops early if the best
measured result reaches ``target_ppf`` and meets the keff constraint, when
configured. Each proposal waits for the previous OpenMC run and GP refit;
this implementation does not evaluate BO candidates in parallel. The
returned design is the best feasible measured candidate (or the closest-keff
measured candidate if none is feasible), not the GP's unverified prediction.

GP surrogates generally become less informative as the number of design
variables rises, because a finite number of samples covers less of the
growing design space. BO does not guarantee a global optimum. Increasing the
initial design and proposal budget or comparing multiple random seeds can
provide a broader search, at the cost of additional OpenMC runs.
