#!/usr/bin/env python3
"""
根据 Szałowski, Phys. Rev. B 74, 094501 (2006) 的正交紧束缚模型，
通过蒙特卡洛采样重复 MgB₂ 体材料的 DOS 与两带 BCS 参数。

主要步骤
---------
1. 在第一布里渊区均匀随机采样 k 点，计算 π、σ 能带能量。
2. 求解费米能量，使四个子带总占据数 n_e ≈ 2.9（文中给出的体相值）。
3. 在获得的费米能附近估算单自旋态密度 gσ、gπ。
4. 调整两带 BCS 耦合 λ = V·g，使 Δσ(0)=7.1 meV、Δπ(0)=2.80 meV、Tc=37.6 K。

注意：这里的 σ、π 命名沿用原文（σ → sp² 衍生空穴带，π → p_z 带）。
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import Tuple

import numpy as np

K_B_EV = 8.617333262145e-5  # eV/K


@dataclass(frozen=True)
class Lattice:
    a: float  # Å
    c: float  # Å

    @property
    def direct_vectors(self) -> np.ndarray:
        a1 = np.array([self.a, 0.0, 0.0])
        a2 = np.array([self.a / 2.0, math.sqrt(3) * self.a / 2.0, 0.0])
        a3 = np.array([0.0, 0.0, self.c])
        return np.column_stack([a1, a2, a3])

    @property
    def reciprocal_vectors(self) -> np.ndarray:
        return 2.0 * math.pi * np.linalg.inv(self.direct_vectors).T

    @property
    def bz_volume(self) -> float:
        return abs(np.linalg.det(self.reciprocal_vectors))


def honeycomb_neighbors(a: float) -> np.ndarray:
    """两个子晶格之间的最近邻位矢"""
    return np.array(
        [
            [a / math.sqrt(3), 0.0, 0.0],
            [-a / (2.0 * math.sqrt(3)), a / 2.0, 0.0],
            [-a / (2.0 * math.sqrt(3)), -a / 2.0, 0.0],
        ]
    )


def sample_k_points(n_samples: int, recip: np.ndarray, seed: int | None = None) -> np.ndarray:
    rng = np.random.default_rng(seed)
    coeffs = rng.random((n_samples, 3)) - 0.5  # [-0.5, 0.5)
    return coeffs @ recip.T


def pi_band_energies(k: np.ndarray, params: dict, nn: np.ndarray) -> np.ndarray:
    """p_z 轨道构成的 π 带（键合/反键合）"""
    phase = k @ nn.T
    f_val = np.exp(1j * phase).sum(axis=1)
    cos_term = np.cos(k[:, 2] * params["c"])
    prefactor = params["e"] + 2.0 * params["t_perp"] * cos_term
    magnitude = params["t_par"] * np.abs(f_val)
    bonding = prefactor - magnitude
    antibonding = prefactor + magnitude
    return np.stack([bonding, antibonding], axis=1)


def sigma_band_energies(k: np.ndarray, params: dict) -> np.ndarray:
    """
    σ 带（sp² 派生重/轻空穴）。沿用文中在 A 点附近的抛物线近似。
    """
    kx, ky, kz = k[:, 0], k[:, 1], k[:, 2]
    k_perp_sq = kx * kx + ky * ky
    cos_term = np.cos(kz * params["c"])
    base = params["e"] + 2.0 * (params["t_par1"] + params["t_par2"]) - 2.0 * params["t_perp"] * cos_term
    light = base - params["t_par1"] * (params["a"] ** 2) * k_perp_sq / 4.0
    heavy = base - 3.0 * params["t_par2"] * (params["a"] ** 2) * k_perp_sq / 4.0
    return np.stack([light, heavy], axis=1)


def total_occupancy(pi_energy: np.ndarray, sigma_energy: np.ndarray, fermi_level: float) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    以单位胞单自旋为计数，返回
    total: 四个子带占据数之和
    occ_sigma: σ 重/轻空穴占据
    occ_pi: π 键合/反键合占据
    """
    occ_sigma = (sigma_energy <= fermi_level).mean(axis=0)
    occ_pi = (pi_energy <= fermi_level).mean(axis=0)
    return float(occ_sigma.sum() + occ_pi.sum()), occ_sigma, occ_pi


def find_fermi_level(pi_energy: np.ndarray, sigma_energy: np.ndarray, target_ne: float, tol: float = 1e-5) -> Tuple[float, np.ndarray, np.ndarray]:
    """用二分法求解满足目标总电子数的费米面"""
    e_min = min(float(pi_energy.min()), float(sigma_energy.min())) - 5.0
    e_max = max(float(pi_energy.max()), float(sigma_energy.max())) + 5.0
    total_low, _, _ = total_occupancy(pi_energy, sigma_energy, e_min)
    total_high, occ_sigma, occ_pi = total_occupancy(pi_energy, sigma_energy, e_max)
    if not (total_low <= target_ne <= total_high):
        raise RuntimeError("目标电子数不在能量扫描范围内，请扩大 e_min/e_max。")
    lo, hi = e_min, e_max
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        total, occ_sigma, occ_pi = total_occupancy(pi_energy, sigma_energy, mid)
        if abs(total - target_ne) < tol:
            return mid, occ_sigma, occ_pi
        if total > target_ne:
            hi = mid
        else:
            lo = mid
    return mid, occ_sigma, occ_pi


def estimate_dos(energies_rel: np.ndarray, bz_volume: float, window: float, bins: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    通过直方图估算 DOS。输入需为相对费米能量。
    返回 (能量中心, DOS)，单位 states/(eV·cell·spin)。
    """
    flat = energies_rel.ravel()
    mask = np.abs(flat) <= window
    selected = flat[mask]
    hist, edges = np.histogram(selected, bins=bins, range=(-window, window))
    bin_width = edges[1] - edges[0]
    # 体元权重：BZ 体积 / 采样点数
    weight = bz_volume / energies_rel.shape[0]
    dos = hist * weight / ((2.0 * math.pi) ** 3 * bin_width)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, dos


def interpolate_at_zero(x: np.ndarray, y: np.ndarray) -> float:
    """线性插值估算 y(x=0)"""
    idx = np.argsort(np.abs(x))
    x0, x1 = x[idx[0]], x[idx[1]]
    y0, y1 = y[idx[0]], y[idx[1]]
    if abs(x1 - x0) < 1e-12:
        return float((y0 + y1) / 2.0)
    return float(y0 + (0.0 - x0) * (y1 - y0) / (x1 - x0))


def solve_couplings(
    g_sigma: float,
    g_pi: float,
    delta_sigma: float,
    delta_pi: float,
    theta_d: float,
    target_tc: float,
) -> np.ndarray:
    """
    求解 λ = V·g，满足 0K 及 Tc 方程。
    λ_ij = g_j V_ij；返回 2×2 λ。
    """
    e_d = K_B_EV * theta_d
    f_sigma = math.asinh(e_d / delta_sigma)
    f_pi = math.asinh(e_d / delta_pi)

    a = delta_sigma * f_sigma
    b = delta_pi * f_pi
    rhs_sigma = delta_sigma / 2.0
    rhs_pi = delta_pi / 2.0

    def integral_tanh(tc: float) -> float:
        # Simpson 积分
        n_steps = 2000
        x_max = e_d
        xs = np.linspace(0.0, x_max, n_steps + 1)
        xs[0] = 1e-12  # 避免除零
        integrand = np.tanh(xs / (2.0 * K_B_EV * tc)) / xs
        dx = x_max / n_steps
        value = integrand[0] + integrand[-1] + 4.0 * integrand[1:-1:2].sum() + 2.0 * integrand[2:-2:2].sum()
        return value * dx / 3.0

    lambda_target = 1.0 / integral_tanh(target_tc)

    lsp_min, lsp_max = 1e-4, 5.0  # 只考虑正值
    for _ in range(200):
        lsp = 0.5 * (lsp_min + lsp_max)
        lss = (rhs_sigma - b * lsp) / a
        lpp = (rhs_pi - a * lsp) / b
        matrix = np.array([[lss, lsp], [lsp, lpp]])
        eigvals = np.linalg.eigvals(matrix)
        max_eig = float(np.max(eigvals.real))
        if max_eig > lambda_target:
            lsp_max = lsp
        else:
            lsp_min = lsp
    lsp = 0.5 * (lsp_min + lsp_max)
    lss = (rhs_sigma - b * lsp) / a
    lpp = (rhs_pi - a * lsp) / b
    return np.array([[lss, lsp], [lsp, lpp]])


def main() -> None:
    parser = argparse.ArgumentParser(description="蒙特卡洛复现 MgB₂ 紧束缚 DOS 及两带 BCS 参数")
    parser.add_argument("-n", "--samples", type=int, default=200_000, help="随机 k 点数量（默认 200000）")
    parser.add_argument("-b", "--bins", type=int, default=240, help="DOS 直方图分箱数（默认 240）")
    parser.add_argument("-w", "--window", type=float, default=0.18, help="DOS 评估窗口（默认 ±0.18 eV）")
    parser.add_argument("--seed", type=int, default=42, help="随机数种子")
    args = parser.parse_args()

    lattice = Lattice(a=3.083, c=3.521)
    recip = lattice.reciprocal_vectors
    bz_volume = lattice.bz_volume
    nn = honeycomb_neighbors(lattice.a)

    k_points = sample_k_points(args.samples, recip, args.seed)

    pi_params = {"e": 0.04, "t_perp": 0.92, "t_par": 1.60, "c": lattice.c}
    sigma_params = {"e": -12.62, "t_perp": 0.094, "t_par1": 5.69, "t_par2": 0.91, "a": lattice.a, "c": lattice.c}

    pi_energy = pi_band_energies(k_points, pi_params, nn)
    sigma_energy = sigma_band_energies(k_points, sigma_params)

    target_ne = 2.9  # 文中体相单自旋电子数
    fermi_level, occ_sigma, occ_pi = find_fermi_level(pi_energy, sigma_energy, target_ne)

    print(f"求得费米能 E_F ≈ {fermi_level:.6f} eV")
    print("各子带占据数（单自旋）：")
    print(f"  σ: heavy {occ_sigma[1]:.3f}, light {occ_sigma[0]:.3f}")
    print(f"  π: bonding {occ_pi[0]:.3f}, antibonding {occ_pi[1]:.3f}")
    print(f"  总计 n_e ≈ {occ_sigma.sum() + occ_pi.sum():.3f}\n")

    pi_rel = pi_energy - fermi_level
    sigma_rel = sigma_energy - fermi_level

    centers_pi, dos_pi = estimate_dos(pi_rel, bz_volume, args.window, args.bins)
    centers_sigma, dos_sigma = estimate_dos(sigma_rel, bz_volume, args.window, args.bins)

    g_pi = interpolate_at_zero(centers_pi, dos_pi)
    g_sigma = interpolate_at_zero(centers_sigma, dos_sigma)

    print("单自旋 DOS（states / (eV·cell·spin））：")
    print(f"  g_sigma(E_F) ≈ {g_sigma:.4f}")
    print(f"  g_pi(E_F)    ≈ {g_pi:.4f}\n")

    delta_sigma = 7.1e-3  # eV
    delta_pi = 2.80e-3  # eV
    theta_d = 1050.0  # K
    target_tc = 37.6  # K

    lambda_matrix = solve_couplings(
        g_sigma=g_sigma,
        g_pi=g_pi,
        delta_sigma=delta_sigma,
        delta_pi=delta_pi,
        theta_d=theta_d,
        target_tc=target_tc,
    )

    print("λ = V·g（无量纲）矩阵：")
    print(lambda_matrix)

    v_matrix = np.array(
        [
            [lambda_matrix[0, 0] / g_sigma, lambda_matrix[0, 1] / g_pi],
            [lambda_matrix[1, 0] / g_sigma, lambda_matrix[1, 1] / g_pi],
        ]
    )
    print("\n对应的 BCS 相互作用常数 V_ij（eV·cell·spin）：")
    print(v_matrix)


if __name__ == "__main__":
    main()

