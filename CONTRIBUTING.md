# Contributing to ommflow

## Setting up

```bash
bash install.sh                      # conda env with automatic ligands
bash install.sh --pip-only           # venv, everything except GAFF ligands
```

`install.sh` reports what the environment can do and runs the test suite. See
[README](README.md#installation) for what the two modes install and why the
dependency stack is split.

## Working on a change

`main` is protected: it takes changes only through a pull request whose checks
have passed. Work on a branch and open a pull request against `main`.

```bash
git switch -c short-descriptive-branch-name
pytest tests -q                      # the same suite CI runs
python -m sphinx -b html -W docs docs/_build/html   # docs, warnings as errors
```

Pull requests are squash-merged, so the pull request title and body become the
single commit message on `main`. Write them for someone reading `git log` a
year from now; the template prompts for what belongs there.

## What the tests hold you to

Beyond correctness, the suite enforces some project conventions, so a change
can fail CI for a documentation reason:

- every configuration setting and command-line option must appear in both the
  README and `docs/`
- every TOML block in the documentation must parse and match the real defaults
- the Sphinx version must match `pyproject.toml`, the license holder must be
  consistent, and ReadTheDocs must build on a supported Python
- the installed OpenMM must provide every force field the package offers

Adding a setting therefore means updating the configuration template,
`final.toml`, the CLI reference, and both documentation surfaces.

## Repository settings

`scripts/configure-github.sh` applies the branch protection and merge settings
described above. It is idempotent and prints what it changed.
