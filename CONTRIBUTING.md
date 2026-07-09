# Contributing

Thank you for your interest in this repository. Because the code archived here corresponds to a peer-reviewed publication, the canonical source is **frozen at the version-of-record DOI on Zenodo**. This document describes how to propose improvements without changing the published artefact.

## Reporting issues

- Use the GitHub issue tracker for bug reports, reproducibility problems, or documentation gaps.
- Please include: operating system, Python version, package versions (`pip freeze`), the exact command you ran, the full traceback, and (if possible) a minimal example that reproduces the problem.

## Proposing fixes

1. Fork the repository.
2. Create a feature branch: `git checkout -b fix/short-description`.
3. Make your changes, keeping the analytical seeds (`RANDOM_STATE = 42`) and stage interfaces unchanged unless explicitly justified.
4. Run the full pipeline to confirm reproducibility of the published numbers, or document any deviation in `CHANGELOG.md`.
5. Open a pull request describing the change and its rationale.

## Coding standards

- Python ≥ 3.12.
- Follow [PEP 8](https://peps.python.org/pep-0008/); we use [`ruff`](https://github.com/astral-sh/ruff) for linting.
- Keep stage scripts self-contained: a stage may depend on outputs of earlier stages, but should not reach into a later stage's namespace.
- Documentation strings should follow NumPy style.

## Versioning

We use [Semantic Versioning](https://semver.org/). Patch releases (`x.y.Z`) for documentation-only or non-functional fixes, minor releases (`x.Y.0`) for additive changes that do not affect published numbers, major releases (`X.0.0`) for any change that alters reported results.

## Code of conduct

By participating, you agree to uphold a respectful and inclusive collaborative environment. Personal attacks, harassment, or discriminatory language will not be tolerated.
