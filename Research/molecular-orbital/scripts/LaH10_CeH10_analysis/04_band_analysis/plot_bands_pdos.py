#!/usr/bin/env python3
"""
绘制能带结构+投影态密度（PDOS）
适用于LaH10和CeH10的对比分析
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os


def read_eigenval(filename):
    """
    读取EIGENVAL文件，提取能带结构

    返回:
        kpoints: k点坐标数组 (nkpts, 3)
        energies: 能带能量数组 (nkpts, nbands)
        efermi: 费米能级（从DOSCAR读取，这里先设为0）
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # 解析头部
    # 第6行：unknown, nkpts, nbands
    nelect, nkpts, nbands = map(int, lines[5].split())

    kpoints = []
    energies = []

    idx = 7  # 数据从第8行开始
    for ikpt in range(nkpts):
        # k点坐标和权重
        kpt_line = lines[idx].split()
        kpt = [float(x) for x in kpt_line[:3]]
        kpoints.append(kpt)
        idx += 1

        # 能带能量
        bands_at_kpt = []
        for ib in range(nbands):
            energy = float(lines[idx].split()[1])  # 第2列是能量
            bands_at_kpt.append(energy)
            idx += 1

        energies.append(bands_at_kpt)

    return np.array(kpoints), np.array(energies)


def read_doscar(filename):
    """
    读取DOSCAR文件，提取总态密度和费米能级

    返回:
        energies: 能量数组 (nedos,)
        dos: 态密度数组 (nedos,)
        efermi: 费米能级（eV）
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # 第6行：Emax Emin NEDOS Efermi weight
    params = lines[5].split()
    efermi = float(params[3])
    nedos = int(params[2])

    # DOS数据（从第7行开始）
    energies = []
    dos_total = []
    for i in range(6, 6 + nedos):
        parts = lines[i].split()
        energies.append(float(parts[0]))
        dos_total.append(float(parts[1]) + float(parts[2]))  # 上自旋+下自旋

    return np.array(energies), np.array(dos_total), efermi


def read_procar(filename, nbands, nkpts, natoms):
    """
    读取PROCAR文件，提取原子和轨道投影

    返回:
        projections: 投影数组 (nkpts, nbands, natoms, norbitals)
    """
    # PROCAR格式复杂，这里简化处理
    # 实际使用时建议用pymatgen或其他库
    print("警告：PROCAR解析功能未实现，请使用pymatgen")
    return None


def calculate_kpath_distance(kpoints):
    """
    计算k路径的累积距离

    参数:
        kpoints: k点坐标 (nkpts, 3)

    返回:
        distances: 累积距离 (nkpts,)
    """
    distances = [0.0]
    for i in range(1, len(kpoints)):
        dk = np.linalg.norm(kpoints[i] - kpoints[i-1])
        distances.append(distances[-1] + dk)
    return np.array(distances)


def plot_band_structure(kpoints, energies, efermi, material_name="Material"):
    """
    绘制能带结构

    参数:
        kpoints: k点坐标 (nkpts, 3)
        energies: 能带能量 (nkpts, nbands)
        efermi: 费米能级（eV）
        material_name: 材料名称（用于标题）
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # 计算k路径距离
    distances = calculate_kpath_distance(kpoints)

    # 将能量相对于费米能级
    energies_rel = energies - efermi

    # 绘制所有能带
    nkpts, nbands = energies.shape
    for ib in range(nbands):
        ax.plot(distances, energies_rel[:, ib], 'k-', linewidth=0.8, alpha=0.7)

    # 标记费米能级
    ax.axhline(0, color='r', linestyle='--', linewidth=2, label='$E_F$')

    # 高对称点标记（需要手动指定）
    # 这里假设是FCC的Γ-X-M-Γ-R-X路径
    # 实际使用时需要根据KPOINTS文件调整
    high_sym_labels = ['$\\Gamma$', 'X', 'M', '$\\Gamma$', 'R', 'X']
    high_sym_positions = [distances[0], distances[len(distances)//5],
                          distances[2*len(distances)//5], distances[3*len(distances)//5],
                          distances[4*len(distances)//5], distances[-1]]

    for pos in high_sym_positions:
        ax.axvline(pos, color='gray', linestyle='-', linewidth=0.5)

    ax.set_xticks(high_sym_positions)
    ax.set_xticklabels(high_sym_labels)
    ax.set_ylabel('Energy - $E_F$ (eV)', fontsize=14)
    ax.set_ylim(-6, 6)
    ax.set_title(f'{material_name} Band Structure', fontsize=16, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{material_name}_band_structure.png', dpi=300)
    print(f"能带结构图已保存：{material_name}_band_structure.png")

    # 统计费米面交叉次数
    crossings = 0
    for ib in range(nbands):
        for ik in range(nkpts - 1):
            if energies_rel[ik, ib] * energies_rel[ik+1, ib] < 0:
                crossings += 1

    print(f"\n{material_name}费米面交叉次数：{crossings}")

    return crossings


def plot_dos(energies, dos, efermi, material_name="Material"):
    """
    绘制态密度

    参数:
        energies: 能量数组 (nedos,)
        dos: 态密度数组 (nedos,)
        efermi: 费米能级（eV）
        material_name: 材料名称
    """
    fig, ax = plt.subplots(figsize=(6, 8))

    # 相对费米能级
    energies_rel = energies - efermi

    # 绘制DOS
    ax.plot(dos, energies_rel, 'b-', linewidth=2, label='Total DOS')
    ax.axhline(0, color='r', linestyle='--', linewidth=2, label='$E_F$')

    ax.set_xlabel('DOS (states/eV)', fontsize=14)
    ax.set_ylabel('Energy - $E_F$ (eV)', fontsize=14)
    ax.set_ylim(-6, 6)
    ax.set_title(f'{material_name} Density of States', fontsize=16, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(f'{material_name}_dos.png', dpi=300)
    print(f"态密度图已保存：{material_name}_dos.png")

    # 计算费米能处的态密度
    idx_fermi = np.argmin(np.abs(energies_rel))
    N_EF = dos[idx_fermi]
    print(f"\n费米能处的态密度 N(EF): {N_EF:.3f} states/eV")

    return N_EF


def main():
    """主程序"""
    if len(sys.argv) < 2:
        print("用法：python plot_bands_pdos.py <材料名称> [EIGENVAL路径] [DOSCAR路径]")
        print("示例：python plot_bands_pdos.py LaH10")
        print("      python plot_bands_pdos.py CeH10 ./CeH10/EIGENVAL ./CeH10/DOSCAR")
        sys.exit(1)

    material_name = sys.argv[1]
    eigenval_file = sys.argv[2] if len(sys.argv) > 2 else "EIGENVAL"
    doscar_file = sys.argv[3] if len(sys.argv) > 3 else "DOSCAR"

    print("=" * 70)
    print(f"{material_name}电子结构分析")
    print("=" * 70)

    # 检查文件
    if not os.path.exists(eigenval_file):
        print(f"错误：文件不存在 {eigenval_file}")
        sys.exit(1)

    if not os.path.exists(doscar_file):
        print(f"错误：文件不存在 {doscar_file}")
        sys.exit(1)

    # 读取能带结构
    print(f"\n读取能带结构：{eigenval_file}")
    kpoints, energies = read_eigenval(eigenval_file)
    print(f"  K点数：{len(kpoints)}")
    print(f"  能带数：{energies.shape[1]}")

    # 读取态密度
    print(f"\n读取态密度：{doscar_file}")
    dos_energies, dos_values, efermi = read_doscar(doscar_file)
    print(f"  费米能级：{efermi:.4f} eV")
    print(f"  DOS能量点数：{len(dos_energies)}")

    # 绘制能带结构
    print("\n绘制能带结构...")
    crossings = plot_band_structure(kpoints, energies, efermi, material_name)

    # 绘制态密度
    print("\n绘制态密度...")
    N_EF = plot_dos(dos_energies, dos_values, efermi, material_name)

    # 保存统计信息
    with open(f"{material_name}_electronic_summary.txt", 'w') as f:
        f.write(f"材料: {material_name}\n")
        f.write(f"费米能级: {efermi:.4f} eV\n")
        f.write(f"费米面交叉次数: {crossings}\n")
        f.write(f"费米能处态密度 N(EF): {N_EF:.3f} states/eV\n")

    print(f"\n统计信息已保存：{material_name}_electronic_summary.txt")

    print("\n" + "=" * 70)
    print("分析完成！")
    print("=" * 70)


if __name__ == "__main__":
    main()
