OpenMC Optimizer documentation
==============================

The OpenMC Optimizer runs a genetic algorithm over design variables in an
OpenMC model. It ranks candidates by a weighted PPF and keff loss, while
requiring candidates to meet the configured keff tolerance. You can adjust
the objective weights, search bounds, and generation settings in Python.

Start with the :doc:`example` for a runnable setup. The :doc:`api` page lists
all configuration options, and :doc:`theory` explains the loss and genetic
algorithm.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   theory
   example
   api

Publishing
----------

Documentation is built with Sphinx and published to GitHub Pages by the
workflow in ``.github/workflows/docs.yml``. To enable publishing for a
repository, open its **Settings > Pages** and set **Build and deployment >
Source** to **GitHub Actions**. Pushes to the repository's default branch
deploy the site; pull requests build the documentation without deploying it.
