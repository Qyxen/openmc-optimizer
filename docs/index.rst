OpenMC Optimizer documentation
==============================

The OpenMC Optimizer runs either a genetic algorithm or constrained Bayesian
optimization over design variables in an OpenMC model. Both methods minimize
a weighted PPF and keff loss for candidates meeting the configured keff
tolerance. You can adjust the objective weights, search bounds, and method
settings in Python.

Start with the :doc:`example` for a runnable setup. The :doc:`api` page lists
all configuration options, and :doc:`theory` explains the loss and both
search methods.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   theory
   example
   vver_example
   api
