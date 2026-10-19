OpenMC Optimizer documentation
==============================

The OpenMC Optimizer runs a genetic algorithm over design variables in an
OpenMC model. Its current objective is to minimize assembly power peaking
factor (PPF).

.. toctree::
   :maxdepth: 2
   :caption: User guide

   api
   example
   theory

Publishing
----------

Documentation is built with Sphinx and published to GitHub Pages by the
workflow in ``.github/workflows/docs.yml``. To enable publishing for a
repository, open its **Settings > Pages** and set **Build and deployment >
Source** to **GitHub Actions**. Pushes to the repository's default branch
deploy the site; pull requests build the documentation without deploying it.
