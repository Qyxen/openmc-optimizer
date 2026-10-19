# Releasing OpenMC Optimizer

Releases are created from version tags. The release workflow builds the
source distribution and wheel, then attaches both files to a GitHub Release.

1. Update the project version in `pyproject.toml` and `docs/conf.py`.
2. Update `CITATION.cff` to the same version.
3. Commit the release changes and push them to the repository.
4. Create and push a matching version tag, for example `v0.1.0`:

   ```bash
   git tag v0.1.0
   git push origin v0.1.0
   ```

GitHub Actions builds the package and creates a release with generated release
notes and the built distribution files. The tag version should match the
version in `pyproject.toml`.
