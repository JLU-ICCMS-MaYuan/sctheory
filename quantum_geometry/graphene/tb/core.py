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

import matplotlib.pyplot as plt
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
    # Keep a tiny sublattice offset only as a numerical regulator.
    # The graphene Dirac-limit formulas in the paper assume the gapless case.
    onsite_delta_ev: float = 1.0e-4

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


def reduced_k_to_cartesian(k_pts: np.ndarray, recip_lat_vecs: np.ndarray) -> np.ndarray:
    """Convert reduced reciprocal coordinates to Cartesian reciprocal coordinates."""
    return np.asarray(k_pts) @ np.asarray(recip_lat_vecs)


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
    k_pts = np.asarray(model.k_uniform_mesh([mesh_size, mesh_size], include_endpoints=False))
    evals = np.asarray(model.solve_ham(k_pts))
    gap = evals[:, 1] - evals[:, 0]
    metric = np.asarray(
        model.quantum_metric(k_pts=k_pts, occ_idxs=[0], cartesian=True)
    )
    berry = np.asarray(
        model.berry_curvature(k_pts=k_pts, occ_idxs=[0], cartesian=True)
    )
    metric_trace = metric[0, 0, :] + metric[1, 1, :]
    recip_lat_vecs = np.asarray(model.recip_lat_vecs)
    k_pts_cart = reduced_k_to_cartesian(k_pts, recip_lat_vecs)
    energy_grad_cart = estimate_energy_gradient_cartesian(
        evals[:, 0], mesh_size, recip_lat_vecs
    )
    velocity_norm = np.linalg.norm(energy_grad_cart, axis=1)
    bz_area = abs(np.linalg.det(recip_lat_vecs))
    dk_area = bz_area / (mesh_size**2)
    return {
        "k_pts": k_pts,
        "k_pts_cart": k_pts_cart,
        "evals_ev": evals,
        "gap_ev": gap,
        "metric_tensor": metric,
        "metric_trace": metric_trace,
        "berry_curvature": berry[0, 1, :],
        "mesh_size": mesh_size,
        "recip_lat_vecs": recip_lat_vecs,
        "energy_grad_cart": energy_grad_cart,
        "velocity_norm": velocity_norm,
        "dk_area": dk_area,
    }


def estimate_energy_gradient_cartesian(
    lower_band_energies: np.ndarray,
    mesh_size: int,
    recip_lat_vecs: np.ndarray,
) -> np.ndarray:
    """Estimate Cartesian energy gradients on a periodic reduced-coordinate mesh."""
    lower = lower_band_energies.reshape(mesh_size, mesh_size)
    du = 1.0 / mesh_size
    dv = 1.0 / mesh_size
    d_edu = (np.roll(lower, -1, axis=0) - np.roll(lower, 1, axis=0)) / (2.0 * du)
    d_edv = (np.roll(lower, -1, axis=1) - np.roll(lower, 1, axis=1)) / (2.0 * dv)
    grad_reduced = np.stack([d_edu, d_edv], axis=-1).reshape(-1, 2)
    transform = np.linalg.inv(np.asarray(recip_lat_vecs))
    return grad_reduced @ transform.T


def extract_fermi_surface_segments(lower_band: np.ndarray, mu_ev: float, mesh_size: int) -> list[np.ndarray]:
    """Extract Fermi-surface contour segments in reduced reciprocal coordinates."""
    u = np.arange(mesh_size) / mesh_size
    v = np.arange(mesh_size) / mesh_size
    uu, vv = np.meshgrid(u, v, indexing="ij")
    fig, ax = plt.subplots()
    contour = ax.contour(uu, vv, lower_band.reshape(mesh_size, mesh_size), levels=[mu_ev])
    segments = [seg for seg in contour.allsegs[0] if len(seg) >= 2]
    plt.close(fig)
    return segments


def evaluate_band_quantities_on_points(
    model,
    points_reduced: np.ndarray,
    recip_lat_vecs: np.ndarray,
    reduced_step: float = 1.0e-5,
) -> dict[str, np.ndarray]:
    """Evaluate energies, gap, metric trace, and Cartesian velocity on arbitrary reduced k points."""
    wrapped = np.mod(points_reduced, 1.0)
    evals = np.asarray(model.solve_ham(wrapped))
    gap = evals[:, 1] - evals[:, 0]
    metric = np.asarray(model.quantum_metric(k_pts=wrapped, occ_idxs=[0], cartesian=True))
    metric_trace = metric[0, 0, :] + metric[1, 1, :]

    shift_u = np.array([reduced_step, 0.0])
    shift_v = np.array([0.0, reduced_step])
    evals_u_plus = np.asarray(model.solve_ham(np.mod(wrapped + shift_u, 1.0)))
    evals_u_minus = np.asarray(model.solve_ham(np.mod(wrapped - shift_u, 1.0)))
    evals_v_plus = np.asarray(model.solve_ham(np.mod(wrapped + shift_v, 1.0)))
    evals_v_minus = np.asarray(model.solve_ham(np.mod(wrapped - shift_v, 1.0)))
    d_edu = (evals_u_plus[:, 0] - evals_u_minus[:, 0]) / (2.0 * reduced_step)
    d_edv = (evals_v_plus[:, 0] - evals_v_minus[:, 0]) / (2.0 * reduced_step)
    grad_reduced = np.column_stack([d_edu, d_edv])
    transform = np.linalg.inv(np.asarray(recip_lat_vecs))
    grad_cart = grad_reduced @ transform.T
    velocity_norm = np.linalg.norm(grad_cart, axis=1)
    return {
        "lower_band": evals[:, 0],
        "gap_ev": gap,
        "metric_trace": metric_trace,
        "velocity_norm": velocity_norm,
    }


def integrate_fermi_surface_observables(
    model,
    segments_reduced: list[np.ndarray],
    recip_lat_vecs: np.ndarray,
) -> dict[str, float]:
    """Integrate all Fermi-surface observables needed for lambda decomposition."""
    if not segments_reduced:
        return {"dos_integral": 0.0, "lambda_e_integral": 0.0, "lambda_geo_integral": 0.0, "fs_inv_gap_integral": 0.0}

    total_dos = 0.0
    total_lambda_e = 0.0
    total_lambda_geo = 0.0
    total_fs_inv_gap = 0.0
    for seg in segments_reduced:
        seg_cart = reduced_k_to_cartesian(seg, recip_lat_vecs)
        mid_reduced = 0.5 * (seg[:-1] + seg[1:])
        ds = np.linalg.norm(seg_cart[1:] - seg_cart[:-1], axis=1)
        sampled = evaluate_band_quantities_on_points(model, mid_reduced, recip_lat_vecs)
        safe_velocity = np.maximum(sampled["velocity_norm"], 1.0e-12)
        safe_gap_sq = np.maximum(sampled["gap_ev"] ** 2, 1.0e-12)
        total_dos += float(np.sum(ds / safe_velocity))
        total_lambda_e += float(np.sum(ds * safe_velocity))
        total_lambda_geo += float(
            np.sum(ds * safe_gap_sq * sampled["metric_trace"] / safe_velocity)
        )
        total_fs_inv_gap += float(np.sum(ds * safe_velocity / safe_gap_sq))
    return {
        "dos_integral": total_dos,
        "lambda_e_integral": total_lambda_e,
        "lambda_geo_integral": total_lambda_geo,
        "fs_inv_gap_integral": total_fs_inv_gap,
    }


def estimate_dirac_winding_numbers(model, loop_radius: float = 0.04, num_points: int = 721) -> dict[str, int]:
    """Estimate winding numbers from the phase of the off-diagonal graphene Hamiltonian."""
    centers = {
        "K": np.array([2.0 / 3.0, 1.0 / 3.0]),
        "K_prime": np.array([1.0 / 3.0, 2.0 / 3.0]),
    }
    theta = np.linspace(0.0, 2.0 * math.pi, num_points, endpoint=True)
    windings: dict[str, int] = {}
    for label, center in centers.items():
        loop = np.column_stack(
            [
                center[0] + loop_radius * np.cos(theta),
                center[1] + loop_radius * np.sin(theta),
            ]
        )
        ham = np.asarray(model.hamiltonian(np.mod(loop, 1.0)))
        offdiag = ham[:, 0, 1]
        phases = np.unwrap(np.angle(offdiag))
        winding = int(round((phases[-1] - phases[0]) / (2.0 * math.pi)))
        windings[label] = winding
    return windings


def _solve_dirac_pocket_radius(
    model,
    center_cart: np.ndarray,
    direction_cart: np.ndarray,
    mu_ev: float,
    recip_lat_vecs: np.ndarray,
    initial_upper: float = 0.02,
    max_upper: float = 1.0,
    tol: float = 1.0e-8,
) -> float | None:
    """Find the radial distance from a Dirac point to the Fermi contour along a direction."""
    inv_recip = np.linalg.inv(np.asarray(recip_lat_vecs))

    def lower_energy(rho: float) -> float:
        point_cart = center_cart + rho * direction_cart
        point_reduced = np.mod(point_cart @ inv_recip, 1.0)
        evals = np.asarray(model.solve_ham(point_reduced))
        return float(evals[0])

    f0 = lower_energy(0.0) - mu_ev
    if f0 <= 0.0:
        return None

    upper = initial_upper
    f_upper = lower_energy(upper) - mu_ev
    while f_upper > 0.0 and upper < max_upper:
        upper *= 1.5
        f_upper = lower_energy(upper) - mu_ev
    if f_upper > 0.0:
        return None

    lower = 0.0
    for _ in range(80):
        mid = 0.5 * (lower + upper)
        f_mid = lower_energy(mid) - mu_ev
        if abs(f_mid) < tol:
            return mid
        if f_mid > 0.0:
            lower = mid
        else:
            upper = mid
    return 0.5 * (lower + upper)


def extract_dirac_pocket_segments(
    model,
    recip_lat_vecs: np.ndarray,
    mu_ev: float,
    num_angles: int = 720,
) -> list[np.ndarray]:
    """Construct high-precision Fermi pockets around K and K' by radial root finding."""
    centers_reduced = {
        "K": np.array([2.0 / 3.0, 1.0 / 3.0]),
        "K_prime": np.array([1.0 / 3.0, 2.0 / 3.0]),
    }
    centers_cart = {
        label: reduced_k_to_cartesian(center[None, :], recip_lat_vecs)[0]
        for label, center in centers_reduced.items()
    }

    theta = np.linspace(0.0, 2.0 * math.pi, num_angles, endpoint=False)
    segments: list[np.ndarray] = []
    for label, center_cart in centers_cart.items():
        points_cart = []
        for angle in theta:
            direction = np.array([math.cos(angle), math.sin(angle)])
            rho = _solve_dirac_pocket_radius(
                model,
                center_cart=center_cart,
                direction_cart=direction,
                mu_ev=mu_ev,
                recip_lat_vecs=recip_lat_vecs,
            )
            if rho is None:
                points_cart = []
                break
            points_cart.append(center_cart + rho * direction)
        if points_cart:
            points_cart_arr = np.asarray(points_cart)
            points_cart_arr = np.vstack([points_cart_arr, points_cart_arr[0]])
            inv_recip = np.linalg.inv(np.asarray(recip_lat_vecs))
            points_reduced = np.mod(points_cart_arr @ inv_recip, 1.0)
            segments.append(points_reduced)
    return segments


def lambda_decomposition_scan(
    model,
    params: GrapheneParameters,
    mesh: dict[str, np.ndarray],
    chemical_potentials_ev: Iterable[float],
    winding_sum_abs: int = 2,
    dirac_local_cutoff_ev: float = 0.10,
) -> list[LambdaDecomposition]:
    evals = mesh["evals_ev"]
    lower = evals[:, 0]
    mesh_size = int(mesh["mesh_size"])
    recip_lat_vecs = mesh["recip_lat_vecs"]

    common_prefactor = (
        params.unit_cell_area_a2
        * (params.gamma_ev_per_a2**2)
        / ((2.0 * math.pi) ** 2 * params.carbon_mass_ev_s2_per_a2 * params.mean_omega2_ev2)
    )
    topo_prefactor = (
        params.unit_cell_area_a2
        * (params.gamma_ev_per_a2**2)
        / (4.0 * params.carbon_mass_ev_s2_per_a2 * params.mean_omega2_ev2)
    )
    results: list[LambdaDecomposition] = []
    for mu_ev in chemical_potentials_ev:
        if abs(mu_ev - params.onsite_delta_ev) <= dirac_local_cutoff_ev:
            segments = extract_dirac_pocket_segments(model, recip_lat_vecs, mu_ev)
        else:
            segments = extract_fermi_surface_segments(lower, mu_ev, mesh_size)
        integrals = integrate_fermi_surface_observables(model, segments, recip_lat_vecs)
        dos = (
            params.unit_cell_area_a2
            / ((2.0 * math.pi) ** 2)
            * integrals["dos_integral"]
        )
        lambda_e = common_prefactor * integrals["lambda_e_integral"]
        lambda_geo = common_prefactor * integrals["lambda_geo_integral"]
        fs_inv_gap_integral = integrals["fs_inv_gap_integral"]
        lambda_topo = (
            topo_prefactor * (winding_sum_abs**2) / fs_inv_gap_integral
            if fs_inv_gap_integral > 0.0
            else 0.0
        )
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
