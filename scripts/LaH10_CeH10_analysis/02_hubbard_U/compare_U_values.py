#!/usr/bin/env python3
"""
对比不同Hubbard U值的计算结果
分析DFT+U对CeH10的影响
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

def extract_energy_from_outcar(outcar_file):
    """从OUTCAR提取总能量"""
    with open(outcar_file, 'r') as f:
        for line in reversed(f.readlines()):
            if 'free  energy   TOTEN' in line:
                return float(line.split()[-2])
    return None


def extract_dos_at_fermi(doscar_file):
    """
    从DOSCAR提取费米能处的态密度

    DOSCAR格式：
    - 第6行包含：Emax Emin NEDOS Efermi weight
    - 之后是能量和DOS数据
    """
    with open(doscar_file, 'r') as f:
        lines = f.readlines()

    # 读取参数
    params = lines[5].split()
    efermi = float(params[3])
    nedos = int(params[2])

    # 读取DOS数据（从第7行开始）
    energies = []
    dos_values = []
    for i in range(6, 6 + nedos):
        parts = lines[i].split()
        energies.append(float(parts[0]))
        dos_values.append(float(parts[1]))  # 总DOS（上自旋）

    energies = np.array(energies)
    dos_values = np.array(dos_values)

    # 找到最接近费米能的点
    idx_fermi = np.argmin(np.abs(energies - efermi))
    dos_at_ef = dos_values[idx_fermi]

    return efermi, dos_at_ef


def extract_bandgap_near_fermi(eigenval_file, efermi, window=1.0):
    """
    从EIGENVAL检查费米能附近的能隙

    参数:
        eigenval_file: EIGENVAL文件路径
        efermi: 费米能级（eV）
        window: 检查费米能±window范围内的能隙

    返回:
        gap: 能隙大小（eV），如果无能隙则返回0
    """
    if not os.path.exists(eigenval_file):
        return 0.0

    with open(eigenval_file, 'r') as f:
        lines = f.readlines()

    # 解析头部信息
    nelect, nkpts, nbands = map(int, lines[5].split()[:3])

    # 读取所有k点的能带能量
    all_bands = []
    idx = 6  # 数据从第7行开始
    for ikpt in range(nkpts):
        # 跳过k点坐标行
        idx += 2
        bands = []
        for ib in range(nbands):
            energy = float(lines[idx].split()[1])  # 第二列是能量
            bands.append(energy - efermi)
            idx += 1
        all_bands.append(bands)

    all_bands = np.array(all_bands)  # shape: (nkpts, nbands)

    # 找到费米能附近的最高占据态和最低未占据态
    # 假设占据态<0，未占据态>0
    occupied = all_bands[all_bands < 0]
    unoccupied = all_bands[all_bands > 0]

    if len(occupied) > 0 and len(unoccupied) > 0:
        vbm = np.max(occupied)  # 价带顶
        cbm = np.min(unoccupied)  # 导带底
        gap = cbm - vbm
        return gap if gap > 0 else 0.0
    else:
        return 0.0  # 金属


def main():
    """主程序"""
    print("=" * 70)
    print("Hubbard U扫描结果分析")
    print("=" * 70)

    # 定义U值列表（需要根据实际计算修改）
    # 假设目录结构：U0/, U3/, U5/, U7/
    U_dirs = {
        "0": "U0",
        "3": "U3",
        "5": "U5",
        "7": "U7"
    }

    results = {}

    # 提取每个U值的结果
    for U_str, dirname in U_dirs.items():
        if not os.path.exists(dirname):
            print(f"警告：目录不存在 {dirname}，跳过")
            continue

        U = float(U_str)
        results[U] = {}

        # 能量
        outcar = os.path.join(dirname, "OUTCAR")
        if os.path.exists(outcar):
            E = extract_energy_from_outcar(outcar)
            results[U]["energy"] = E
            print(f"\nU = {U} eV:")
            print(f"  能量: {E:.6f} eV")
        else:
            print(f"警告：{outcar} 不存在")
            results[U]["energy"] = None

        # 费米能处的态密度
        doscar = os.path.join(dirname, "DOSCAR")
        if os.path.exists(doscar):
            efermi, dos_ef = extract_dos_at_fermi(doscar)
            results[U]["efermi"] = efermi
            results[U]["dos_ef"] = dos_ef
            print(f"  费米能: {efermi:.3f} eV")
            print(f"  N(EF): {dos_ef:.3f} states/eV")
        else:
            print(f"警告：{doscar} 不存在")
            results[U]["efermi"] = None
            results[U]["dos_ef"] = None

        # 能隙
        eigenval = os.path.join(dirname, "EIGENVAL")
        if os.path.exists(eigenval) and results[U]["efermi"] is not None:
            gap = extract_bandgap_near_fermi(eigenval, results[U]["efermi"])
            results[U]["bandgap"] = gap
            print(f"  能隙: {gap:.4f} eV")
        else:
            results[U]["bandgap"] = 0.0

    # 数据准备（用于绘图）
    U_values = sorted([U for U in results.keys()])
    energies = [results[U]["energy"] for U in U_values if results[U]["energy"] is not None]
    dos_at_ef = [results[U]["dos_ef"] for U in U_values if results[U]["dos_ef"] is not None]
    bandgaps = [results[U]["bandgap"] for U in U_values if results[U]["bandgap"] is not None]

    # 归一化能量（相对于U=0）
    if len(energies) > 0:
        E_ref = energies[0]
        energies_rel = [(E - E_ref) * 1000 for E in energies]  # meV
    else:
        energies_rel = []

    # 绘图
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # 子图1：能量 vs U
    if len(energies_rel) > 0:
        axes[0].plot(U_values, energies_rel, 'o-', linewidth=2, markersize=8)
        axes[0].axhline(0, color='gray', linestyle='--', linewidth=1)
        axes[0].set_xlabel('Hubbard U (eV)', fontsize=14)
        axes[0].set_ylabel('Relative Energy (meV)', fontsize=14)
        axes[0].set_title('Total Energy vs U', fontsize=16)
        axes[0].grid(alpha=0.3)

    # 子图2：态密度 vs U
    if len(dos_at_ef) > 0:
        axes[1].plot(U_values, dos_at_ef, 'o-', linewidth=2, markersize=8, color='green')
        axes[1].set_xlabel('Hubbard U (eV)', fontsize=14)
        axes[1].set_ylabel('N(EF) (states/eV)', fontsize=14)
        axes[1].set_title('DOS at Fermi Level vs U', fontsize=16)
        axes[1].grid(alpha=0.3)

    # 子图3：能隙 vs U
    if len(bandgaps) > 0:
        axes[2].plot(U_values, bandgaps, 'o-', linewidth=2, markersize=8, color='red')
        axes[2].set_xlabel('Hubbard U (eV)', fontsize=14)
        axes[2].set_ylabel('Band Gap (eV)', fontsize=14)
        axes[2].set_title('Gap at Fermi Level vs U', fontsize=16)
        axes[2].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('U_scan_results.png', dpi=300)
    print(f"\n图像已保存：U_scan_results.png")

    # 分析与建议
    print("\n" + "=" * 70)
    print("分析与建议：")
    print("=" * 70)

    if len(dos_at_ef) >= 2:
        dos_change = abs(dos_at_ef[-1] - dos_at_ef[0]) / dos_at_ef[0] * 100
        print(f"\nN(EF)变化率（U=0 → U={U_values[-1]}）：{dos_change:.1f}%")

        if dos_change < 10:
            print("  → 4f对费米能处态密度影响很小")
            print("  → 建议：可以使用标准DFT（GGA），无需DFT+U")

        elif dos_change < 30:
            print("  → 4f有一定影响，但不强")
            print("  → 建议：使用中等U值（U=3 eV）作为修正")

        else:
            print("  → 4f显著影响费米能附近的电子结构")
            print("  → 建议：必须使用DFT+U，选择U=5-7 eV")

    if len(energies_rel) >= 2:
        E_change = abs(energies_rel[-1] - energies_rel[0])
        print(f"\n能量变化（U=0 → U={U_values[-1]}）：{E_change:.1f} meV/atom")

        if E_change < 50:
            print("  → 能量对U不敏感，系统可能不是强关联")

    # 能隙分析
    if any(g > 0.1 for g in bandgaps):
        print(f"\n⚠ 警告：某些U值下出现能隙（{max(bandgaps):.3f} eV）")
        print("  → 这可能不对：CeH10应该是金属！")
        print("  → 检查：")
        print("     1. K点网格是否够密？")
        print("     2. 是否有磁性导致能隙？")
        print("     3. U值是否过大？")

    print("\n" + "=" * 70)


if __name__ == "__main__":
    main()
