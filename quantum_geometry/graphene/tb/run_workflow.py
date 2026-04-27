#!/usr/bin/env python3
"""Run the graphene model workflow and export plots/data."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from core import (
    GrapheneParameters,
    build_model,
    estimate_dirac_winding_numbers,
    high_symmetry_path,
    lambda_decomposition_scan,
    save_lambda_scan_json,
    uniform_mesh_quantities,
)


def plot_band_structure(path_data: dict, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    for band in path_data["bands_ev"]:
        ax.plot(path_data["k_dist"], band, color="black", lw=1.4)
    for node in path_data["k_node"]:
        ax.axvline(node, color="0.8", lw=0.8)
    ax.axhline(0.0, color="tab:red", lw=0.9, linestyle="--")
    ax.set_xticks(path_data["k_node"], path_data["labels"])
    ax.set_ylabel("Energy (eV)")
    ax.set_title("Graphene tight-binding bands")
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_map(k_pts: np.ndarray, values: np.ndarray, mesh_size: int, title: str, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    sc = ax.scatter(k_pts[:, 0], k_pts[:, 1], c=values, s=9, cmap="viridis")
    ax.set_xlabel(r"$k_1$")
    ax.set_ylabel(r"$k_2$")
    ax.set_title(title)
    fig.colorbar(sc, ax=ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_lambda_scan(scan: list[dict], output_path: Path) -> None:
    mu = np.array([item["chemical_potential_ev"] for item in scan])
    lam_total = np.array([item["lambda_total"] for item in scan])
    lam_e = np.array([item["lambda_E"] for item in scan])
    lam_geo = np.array([item["lambda_geo"] for item in scan])
    lam_topo = np.array([item["lambda_topo"] for item in scan])

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(mu, lam_total, label=r"$\lambda$", color="black", lw=1.6)
    ax.plot(mu, lam_e, label=r"$\lambda_E$", color="tab:blue", lw=1.4)
    ax.plot(mu, lam_geo, label=r"$\lambda_{\mathrm{geo}}$", color="tab:green", lw=1.4)
    ax.plot(mu, lam_topo, label=r"$\lambda_{\mathrm{topo}}$", color="tab:orange", lw=1.4)
    ax.set_xlabel(r"Chemical potential $\mu$ (eV)")
    ax.set_ylabel(r"Dimensionless $\lambda$")
    ax.set_title("Graphene lambda decomposition")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_fraction_scan(scan: list[dict], output_path: Path) -> None:
    mu = np.array([item["chemical_potential_ev"] for item in scan])
    geo = np.array([item["geometric_fraction"] for item in scan])
    topo = np.array([item["topological_fraction_of_geo"] for item in scan])

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(mu, geo, label=r"$\lambda_{\mathrm{geo}}/\lambda$", color="tab:green", lw=1.6)
    ax.plot(
        mu,
        topo,
        label=r"$\lambda_{\mathrm{topo}}/\lambda_{\mathrm{geo}}$",
        color="tab:orange",
        lw=1.6,
    )
    ax.axhline(0.5, color="0.5", lw=0.8, linestyle="--")
    ax.axhline(1.0, color="0.7", lw=0.8, linestyle=":")
    ax.set_xlabel(r"Chemical potential $\mu$ (eV)")
    ax.set_ylabel("Fraction")
    ax.set_title("Geometric and topological fractions")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Graphene model workflow")
    parser.add_argument(
        "--output-dir",
        default=Path(__file__).resolve().parents[1] / "runs" / "initial_model",
        type=Path,
        help="Directory for generated data and plots",
    )
    parser.add_argument("--mesh-size", type=int, default=121, help="Uniform k mesh size")
    parser.add_argument(
        "--mu-min",
        type=float,
        default=-1.0,
        help="Minimum chemical potential in eV for the lambda scan",
    )
    parser.add_argument(
        "--mu-max",
        type=float,
        default=-0.02,
        help="Maximum chemical potential in eV for the lambda scan",
    )
    parser.add_argument("--mu-count", type=int, default=25, help="Number of mu points")
    args = parser.parse_args()

    output_dir = args.output_dir
    plots_dir = output_dir / "plots"
    data_dir = output_dir / "data"
    plots_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    params = GrapheneParameters()
    model = build_model(params)
    path_data = high_symmetry_path(model)
    mesh_data = uniform_mesh_quantities(model, mesh_size=args.mesh_size)
    winding_numbers = estimate_dirac_winding_numbers(model)
    winding_sum_abs = sum(abs(value) for value in winding_numbers.values())
    mu_values = np.linspace(args.mu_min, args.mu_max, args.mu_count)
    scan = lambda_decomposition_scan(
        model, params, mesh_data, mu_values, winding_sum_abs=winding_sum_abs
    )
    scan_dicts = [asdict(item) for item in scan]

    np.savez_compressed(
        data_dir / "graphene_mesh_data.npz",
        k_pts=mesh_data["k_pts"],
        k_pts_cart=mesh_data["k_pts_cart"],
        evals_ev=mesh_data["evals_ev"],
        gap_ev=mesh_data["gap_ev"],
        metric_trace=mesh_data["metric_trace"],
        berry_curvature=mesh_data["berry_curvature"],
        velocity_norm=mesh_data["velocity_norm"],
    )
    np.savez_compressed(
        data_dir / "graphene_band_path.npz",
        k_vec=path_data["k_vec"],
        k_dist=path_data["k_dist"],
        k_node=path_data["k_node"],
        bands_ev=path_data["bands_ev"],
    )
    (data_dir / "parameters.json").write_text(
        json.dumps(asdict(params), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (data_dir / "winding_numbers.json").write_text(
        json.dumps(winding_numbers, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    save_lambda_scan_json(scan, data_dir / "lambda_scan.json")

    plot_band_structure(path_data, plots_dir / "graphene_bands.png")
    plot_map(
        mesh_data["k_pts"],
        mesh_data["metric_trace"],
        args.mesh_size,
        "Quantum metric trace",
        plots_dir / "graphene_metric_trace.png",
    )
    plot_map(
        mesh_data["k_pts"],
        mesh_data["berry_curvature"],
        args.mesh_size,
        "Berry curvature",
        plots_dir / "graphene_berry_curvature.png",
    )
    plot_lambda_scan(scan_dicts, plots_dir / "graphene_lambda_scan.png")
    plot_fraction_scan(scan_dicts, plots_dir / "graphene_lambda_fractions.png")

    summary = {
        "output_dir": str(output_dir),
        "mu_min_ev": float(args.mu_min),
        "mu_max_ev": float(args.mu_max),
        "mu_count": int(args.mu_count),
        "mesh_size": int(args.mesh_size),
        "winding_numbers": winding_numbers,
        "winding_sum_abs": int(winding_sum_abs),
        "sample_small_doping_point": scan_dicts[-1],
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
