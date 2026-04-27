# Repository Guidelines

## Project Structure & Module Organization

This repository is a research workspace rather than a single package. Root-level folders such as `H3S/`, `MH6/`, `MH9/`, and `MH10/` store material-specific calculation data. Core theory and analysis live in `docs/`, `EliashbergEquation/`, `scripts/`, and `reference/`. Three Python subprojects are maintained independently: `DeePTB/` (`dptb/` source), `dftio/` (`dftio/` source, `test/` tests), and `pythtb/` (`pythtb/` source, `tests/` tests). Keep new code close to the owning module; do not mix reusable package code into data directories.

## Build, Test, and Development Commands

Run commands from the relevant subproject root, not from `/`.

- `python "EliashbergEquation/eliashberg_solver.py" --input INPUT --alpha2f ALPHA2F.OUT`: run the standalone Eliashberg solver.
- `cd "DeePTB" && uv sync`: install DeePTB with locked dependencies.
- `cd "dftio" && uv sync --group dev`: install `dftio` plus pytest/doc tooling.
- `cd "pythtb" && pip install -e ".[dev,tests]"`: editable install for local development.
- `cd "dftio" && pytest`: run `dftio` tests.
- `cd "pythtb" && pytest`: run `pythtb` tests.

If a change touches only scripts or notebooks, run the smallest relevant command and document any unverified paths in the PR.

## Coding Style & Naming Conventions

Follow the existing style of the target module. Python uses 4-space indentation, `snake_case` for functions/files, `PascalCase` for classes, and concise module docstrings where needed. Respect local tooling when present: `ruff` and `black` are configured in `pythtb`; pytest settings are embedded in `dftio` and `pythtb` `pyproject.toml`. Prefer small, single-purpose functions and avoid introducing framework-wide abstractions into one-off research scripts.

## Testing Guidelines

Place tests beside the owning package conventions: `dftio/test/test_*.py` and `pythtb/tests/test_*.py`. Reuse existing pytest naming (`test_*`, `Test*`). Add regression tests for parser, IO, or numerical workflow changes; for data-heavy research updates, provide a minimal reproducible input or a before/after result summary.

## Commit & Pull Request Guidelines

Recent history uses short Chinese summaries such as `添加了dftio子模块`, `修改了H6.nb模型`, and `删除了不必要的内容以及追加了.gitignore`. Keep commit messages brief, specific, and action-oriented; one logical change per commit. PRs should include scope, affected paths, validation commands, and representative plots/screenshots when changing analysis output or documentation figures.

## Data & Safety Notes

Large calculation outputs and vendor code snapshots already exist in-tree. Avoid bulk reformatting, renaming dataset folders, or committing regenerated binaries unless the change is intentional and documented. Record external data sources and parameter assumptions in the nearest `README.md` or docs file.
