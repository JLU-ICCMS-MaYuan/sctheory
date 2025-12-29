# GEMINI.md - Project Context & Instructions

This project is a research environment focused on the theoretical analysis of high-pressure hydride superconductors (e.g., $LaH_{10}$, $CeH_{10}$, $H_3S$). It combines Density Functional Theory (DFT), Eliashberg theory, and Group Theory for symmetry analysis and chemical bonding investigation.

## Project Overview

The primary research goal is to understand the microscopic mechanisms governing superconductivity in hydrides, specifically focusing on:
- **Comparative Studies:** Explaining $T_c$ differences (e.g., $LaH_{10}$ vs $CeH_{10}$) via orbital inversion and electron-phonon coupling (EPC) mechanisms.
- **Spectral Analysis:** Solving Eliashberg equations to calculate $T_c$ and superconducting gaps.
- **Symmetry & Bonding:** Using `GTPack` (Mathematica) and Crystal Orbital Hamilton Population (COHP) to analyze bonding environments.

## Technical Stack

- **Simulations:** VASP (5.4+), Quantum ESPRESSO (6.7+), EPW, LOBSTER.
- **Analysis:** Python (3.7+) with `numpy`, `scipy`, `matplotlib`, `pymatgen`.
- **Symmetry Analysis:** Mathematica with `GTPack` (v1.4).
- **Visualization:** VESTA, Phonopy, FermiSurfer.

## Key Workflows

### 1. Eliashberg Equation Solver
Located in `EliashbergEquation/`, this Python solver calculates superconducting properties from the spectral function $\alpha^2F(\omega)$.
- **Input:** `ALPHA2F.OUT` (spectral function), `INPUT` (contains $\mu^*$ and temperature steps).
- **Run:** `python eliashberg_solver.py --input INPUT --alpha2f ALPHA2F.OUT`
- **Output:** $T_c$, McMillan parameters, and frequency-dependent gaps.

### 2. $LaH_{10}/CeH_{10}$ Comparative Analysis
A structured workflow in `scripts/LaH10_CeH10_analysis/` follows the research plan in `docs/04_LaH10_CeH10_research_plan.md`:
1.  **Magnetism Test:** Determine Ce 4f character (`01_magnetic_test`).
2.  **Band Analysis:** Identify orbital inversion and density of states at Fermi level (`04_band_analysis`).
3.  **Frozen Phonon:** Directly probe EPC strength by monitoring band splitting vs atomic displacement (`06_frozen_phonon`).

### 3. Symmetry Analysis (Mathematica)
Using `GTPack` for:
- Constructing tight-binding models.
- Irreducible representation decomposition of orbitals (e.g., f-orbitals in $D_{6h}$).
- Symmetry-adapted wavefunctions.

## Directory Structure

- `docs/`: Comprehensive theoretical guides and research plans.
  - `QUICKSTART.md`: Concise steps for the $LaH_{10}/CeH_{10}$ workflow.
  - `03_theory_comprehensive.md`: Deep dive into Eliashberg theory and EPC.
  - `04_LaH10_CeH10_research_plan.md`: Detailed methodology for current research.
- `EliashbergEquation/`: Core Python solver and sample I/O.
- `scripts/`: Task-specific Python and Shell scripts.
- `GroupTheory-V1_4/`: The `GTPack` Mathematica package.
- `MH10/`, `H3S/`, `MH9/`: Project-specific data, VASP files, and Mathematica notebooks.
- `reference/`: Key literature (PDFs/extracted text) supporting the research hypothesis.
- `Book_Tasks/`: Educational tasks for learning `GTPack`.

## Development Conventions

- **Data Tracking:** Analysis scripts should output logs and plots to their respective subdirectories.
- **Python Style:** Use type hints and docstrings for research scripts to maintain clarity.
- **Mathematica:** Keep notebooks organized by system name (e.g., `LaH10_orbit_projection.nb`).

## Key Files
- `EliashbergEquation/eliashberg_solver.py`: Main numerical tool for $T_c$ calculation.
- `docs/04_LaH10_CeH10_research_plan.md`: The "Source of Truth" for current research directions.
- `request.md`: Log of previous AI interactions and problem-solving history.
- `reference/2024 化学键决定两个看似完全相同超导体的临界温度差 Robert.pdf`: The primary paper inspiring the current $LaH_{10}$ vs $CeH_{10}$ study.
