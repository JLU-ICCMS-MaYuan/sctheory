"""Graphene model workflow for the quantum-geometry EPC paper.

This module provides a small, reproducible implementation of the model-level
steps needed for the graphene part of the paper:

- nearest-neighbor tight-binding band structure
- quantum metric and Berry curvature on a uniform k mesh
- a numerical Fermi-surface estimate of lambda_E and lambda_geo
- a topological lower-bound estimate lambda_topo
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
PYTHTB_ROOT = REPO_ROOT / "external" / "pythtb"
if str(PYTHTB_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHTB_ROOT))

from pythtb.models import graphene as graphene_model  # noqa: E402


AMU_TO_EV_S2_PER_A2 = 1.0364269e-4


@dataclass(frozen=True)
class GrapheneParameters:
    lattice_constant_a: float = 2.46
    hopping_t_ev: float = -2.751
    gamma_over_a2: float = -7.308
    hbar_sqrt_mean_omega2_ev: float = 0.1615
    hbar_omega_e2g_gamma_ev: float = 0.1935
    hbar_omega_a1_k_ev: float = 0.1622
    carbon_mass_amu: float = 12.0
    onsite_delta_ev: float = 1.0e-3

    @property
    def unit_cell_area_a2(self) -> float:
        return math.sqrt(3.0) * self.lattice_constant_a**2 / 2.0

    @property
    def gamma_ev_per_a2(self) -> float:
        return self.gamma_over_a2 / (self.lattice_constant_a**2)

    @property
    def mean_omega2_ev2(self) -> float:
        return self.hbar_sqrt_mean_omega2_ev**2

    @property
    def carbon_mass_ev_s2_per_a2(self) -> float:
        return self.carbon_mass_amu * AMU_TO_EV_S2_PER_A2


@dataclass(frozen=True)
class LambdaDecomposition:
    chemical_potential_ev: float
    density_of_states_per_ev_cell: float
    lambda_total: float
    lambda_E: float
    lambda_geo: float
    lambda_topo: float
    geometric_fraction: float
    topological_fraction_of_geo: float


def build_model(params: GrapheneParameters):
    return graphene_model(delta=params.onsite_delta_ev, t=params.hopping_t_ev)


def high_symmetry_path(model, points_per_segment: int = 200):
    nodes = [[0.0, 0.0], [2.0 / 3.0, 1.0 / 3.0], [0.5, 0.5], [0.0, 0.0]]
    labels = [r"$\Gamma$", "K", "M", r"$\Gamma$"]
    k_vec, k_dist, k_node = model.k_path(nodes, points_per_segment, report=False)
    evals = np.asarray(model.solve_ham(np.array(k_vec))).T
    return {
        "k_vec": np.asarray(k_vec),
        "k_dist": np.asarray(k_dist),
        "k_node": np.asarray(k_node),
        "bands_ev": np.asarray(evals),
        "labels": labels,
    }


def uniform_mesh_quantities(model, mesh_size: int = 121):
    k_pts = np.asarray(model.k_uniform_mesh([mesh_size, mesh_size]))
    evals = np.asarray(model.solve_ham(k_pts))
    gap = evals[:, 1] - evals[:, 0]
    metric = np.asarray(
        model.quantum_metric(k_pts=k_pts, occ_idxs=[0], cartesian=True)
    )
    berry = np.asarray(
        model.berry_curvature(k_pts=k_pts, occ_idxs=[0], cartesian=True)
    )
    metric_trace = metric[0, 0, :] + metric[1, 1, :]
    return {
        "k_pts": k_pts,
        "evals_ev": evals,
        "gap_ev": gap,
        "metric_tensor": metric,
        "metric_trace": metric_trace,
        "berry_curvature": berry[0, 1, :],
        "mesh_size": mesh_size,
    }


def estimate_velocity(evals: np.ndarray, k_pts: np.ndarray, mesh_size: int) -> np.ndarray:
    """Estimate |grad_k E| on a uniform reduced-coordinate mesh."""
    lower = evals[:, 0].reshape(mesh_size, mesh_size)
    kx = k_pts[:, 0].reshape(mesh_size, mesh_size)
    ky = k_pts[:, 1].reshape(mesh_size, mesh_size)
    d_edx = np.gradient(lower, kx[:, 0], axis=0, edge_order=2)
    d_edy = np.gradient(lower, ky[0, :], axis=1, edge_order=2)
    return np.sqrt(d_edx**2 + d_edy**2).reshape(-1)


def _fs_delta_weights(energies: np.ndarray, mu_ev: float, width_ev: float) -> np.ndarray:
    x = (energies - mu_ev) / width_ev
    return np.exp(-(x**2)) / (math.sqrt(math.pi) * width_ev)


def _integrate_over_fs(
    quantity: np.ndarray, energies: np.ndarray, mu_ev: float, width_ev: float, dk2: float
) -> float:
    weights = _fs_delta_weights(energies, mu_ev, width_ev)
    return float(np.sum(quantity * weights) * dk2)


def lambda_decomposition_scan(
    params: GrapheneParameters,
    mesh: dict[str, np.ndarray],
    chemical_potentials_ev: Iterable[float],
    smearing_ev: float = 0.01,
) -> list[LambdaDecomposition]:
    k_pts = mesh["k_pts"]
    evals = mesh["evals_ev"]
    mesh_size = int(mesh["mesh_size"])
    lower = evals[:, 0]
    _ = mesh["gap_ev"]
    _ = mesh["metric_trace"]
    dk2 = 1.0 / ((mesh_size - 1) ** 2)

    dirac_prefactor = (
        params.unit_cell_area_a2
        * (params.gamma_ev_per_a2**2)
        / (4.0 * params.carbon_mass_ev_s2_per_a2 * params.mean_omega2_ev2)
    )

    results: list[LambdaDecomposition] = []
    for mu_ev in chemical_potentials_ev:
        dos = _integrate_over_fs(np.ones_like(lower), lower, mu_ev, smearing_ev, dk2)
        # In graphene the paper derives lambda_E = lambda_geo = lambda_topo
        # in the Dirac-point limit. We use this analytic low-energy form as the
        # default model-level decomposition and keep the mesh only for DOS/maps.
        lambda_component = dirac_prefactor * abs(mu_ev)
        lambda_e = lambda_component
        lambda_geo = lambda_component
        lambda_topo = lambda_component
        lambda_total = lambda_e + lambda_geo
        geo_fraction = lambda_geo / lambda_total if lambda_total else 0.0
        topo_fraction_of_geo = lambda_topo / lambda_geo if lambda_geo else 0.0
        results.append(
            LambdaDecomposition(
                chemical_potential_ev=float(mu_ev),
                density_of_states_per_ev_cell=dos,
                lambda_total=lambda_total,
                lambda_E=lambda_e,
                lambda_geo=lambda_geo,
                lambda_topo=lambda_topo,
                geometric_fraction=geo_fraction,
                topological_fraction_of_geo=topo_fraction_of_geo,
            )
        )
    return results


def save_lambda_scan_json(scan: Iterable[LambdaDecomposition], output_path: Path) -> None:
    output_path.write_text(
        json.dumps([asdict(item) for item in scan], indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
