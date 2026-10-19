# OpenMC Optimizer

A small, project-style optimization driver for OpenMC models. The initial
search method is a bounded real-valued genetic algorithm. Its built-in
objective is to minimize assembly power peaking factor (PPF), optionally while
keeping k-effective within a user-set tolerance of a target. The design
variables are independent of the objective, so the same search can vary
enrichment or Gd₂O₃ rod loading using the built-in material adapter.

## Install

Use an environment with the OpenMC Python package and executable available.
Then install this project:

```bash
python -m pip install -e .
```

## Documentation

For detailed user guild please visit: [OpenMC Optimizer Documentation](https://qyxen.github.io/openmc-optimizer/)

## Limitations in this initial implementation

- Candidate evaluations run sequentially; GA randomness is reproducible from
  the configured seed.
- The built-in PPF tally assumes assembly root-cell IDs can be supplied and
  that `fission-q-recoverable` is appropriate for the model.
- Runtime, memory, and tally precision depend on the user's OpenMC settings.
- Generic material fractions rescale all other entries proportionally; they
  do not model coupled recipe or chemistry constraints.
