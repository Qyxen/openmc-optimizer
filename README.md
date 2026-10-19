# OpenMC Optimizer

OpenMC Optimizer is a Python toolkit for optimizing OpenMC reactor models. It
provides bounded genetic algorithm and constrained Bayesian optimization
methods to minimize the power peaking factor (PPF), with an optional
k-effective target and tolerance constraint. Design variables are independent
of the search method. The built-in material adapter can vary an element or
nuclide fraction in a material, or mix two material recipes by volume. These
generic operations support a wide range of material studies, including fuel
enrichment and absorber loading.

## Install

Use an environment with the OpenMC Python package and executable available.
Then install this project:

```bash
python -m pip install -e .
```

Install the released package:
```bash
pip install openmc-optimizer
```

Scikit-learn is installed as a required dependency, so both optimization
methods are available after the standard install.

## Documentation

For detailed user guild please visit: [OpenMC Optimizer Documentation](https://qyxen.github.io/openmc-optimizer/)

## Citation
Cite this project using the metadata in [CITATION.cff](CITATION.cff).