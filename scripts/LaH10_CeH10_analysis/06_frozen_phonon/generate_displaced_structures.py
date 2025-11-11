#!/usr/bin/env python3
"""
生成冻结声子的位移结构
用于研究声子-电子耦合机制
"""
import numpy as np
import sys
import os


def read_poscar(filename):
    """
    读取VASP的POSCAR文件

    返回:
        lattice: 晶格向量 (3, 3)
        species: 元素列表
        atom_counts: 每种元素的原子数
        positions: 原子坐标 (natoms, 3)，分数坐标
        coord_type: 坐标类型（'Direct' 或 'Cartesian'）
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # 第1行：注释
    comment = lines[0].strip()

    # 第2行：缩放因子
    scale = float(lines[1].strip())

    # 第3-5行：晶格向量
    lattice = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)])
    lattice *= scale

    # 第6行：元素列表
    species = lines[5].split()

    # 第7行：每种元素的原子数
    atom_counts = [int(x) for x in lines[6].split()]
    natoms = sum(atom_counts)

    # 第8行：坐标类型
    coord_type = lines[7].strip()[0].upper()  # 'D' for Direct, 'C' for Cartesian

    # 第9行开始：原子坐标
    positions = []
    for i in range(8, 8 + natoms):
        pos = [float(x) for x in lines[i].split()[:3]]
        positions.append(pos)
    positions = np.array(positions)

    # 如果是笛卡尔坐标，转换为分数坐标
    if coord_type == 'C':
        positions = positions @ np.linalg.inv(lattice.T)

    return lattice, species, atom_counts, positions, comment


def write_poscar(filename, lattice, species, atom_counts, positions, comment="Displaced structure"):
    """
    写入VASP的POSCAR文件

    参数:
        filename: 输出文件名
        lattice: 晶格向量 (3, 3)
        species: 元素列表
        atom_counts: 每种元素的原子数
        positions: 原子坐标 (natoms, 3)，分数坐标
        comment: 注释行
    """
    with open(filename, 'w') as f:
        f.write(f"{comment}\n")
        f.write("1.0\n")  # 缩放因子

        # 晶格向量
        for i in range(3):
            f.write(f"  {lattice[i, 0]:18.14f}  {lattice[i, 1]:18.14f}  {lattice[i, 2]:18.14f}\n")

        # 元素和数量
        f.write("  " + "  ".join(species) + "\n")
        f.write("  " + "  ".join(map(str, atom_counts)) + "\n")

        # 坐标类型
        f.write("Direct\n")

        # 原子坐标
        for pos in positions:
            f.write(f"  {pos[0]:18.14f}  {pos[1]:18.14f}  {pos[2]:18.14f}\n")


def normalize_eigenvector(eigenvector):
    """
    归一化声子本征矢

    参数:
        eigenvector: (natoms, 3)数组

    返回:
        normalized: 归一化的本征矢
    """
    norm = np.linalg.norm(eigenvector)
    if norm < 1e-10:
        print("警告：本征矢的模接近零！")
        return eigenvector
    return eigenvector / norm


def generate_displaced_structures(
    poscar_equilibrium,
    eigenvector,
    amplitudes,
    output_dir="displaced_structures",
    mode_number=0
):
    """
    生成一系列位移结构

    参数:
        poscar_equilibrium: 平衡结构的POSCAR文件路径
        eigenvector: 声子本征矢 (natoms, 3)，单位：归一化
        amplitudes: 位移振幅列表（Å）
        output_dir: 输出目录
        mode_number: 声子模式编号（用于文件命名）
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)

    # 读取平衡结构
    lattice, species, atom_counts, positions_frac, comment = read_poscar(poscar_equilibrium)
    natoms = len(positions_frac)

    # 检查本征矢维度
    if eigenvector.shape != (natoms, 3):
        print(f"错误：本征矢维度 {eigenvector.shape} 与原子数 {natoms} 不匹配！")
        return

    # 归一化本征矢
    eigenvector_norm = normalize_eigenvector(eigenvector.flatten()).reshape(natoms, 3)

    print(f"生成冻结声子结构：模式 #{mode_number}")
    print(f"  原子数：{natoms}")
    print(f"  振幅范围：{min(amplitudes):.3f} - {max(amplitudes):.3f} Å")
    print(f"  输出目录：{output_dir}")

    # 将分数坐标转换为笛卡尔坐标（用于位移）
    positions_cart = positions_frac @ lattice

    # 对每个振幅生成位移结构
    for amp in amplitudes:
        # 计算位移（笛卡尔坐标）
        displacement_cart = amp * eigenvector_norm  # (natoms, 3)

        # 新的笛卡尔坐标
        positions_new_cart = positions_cart + displacement_cart

        # 转换回分数坐标
        positions_new_frac = positions_new_cart @ np.linalg.inv(lattice)

        # 周期性边界条件（确保在0-1范围内）
        positions_new_frac = positions_new_frac % 1.0

        # 保存POSCAR
        output_file = os.path.join(output_dir, f"POSCAR_mode{mode_number}_amp{amp:.4f}")
        new_comment = f"{comment} | Mode {mode_number}, Amplitude {amp:.4f} Å"
        write_poscar(output_file, lattice, species, atom_counts, positions_new_frac, new_comment)

        print(f"  已生成：{output_file}")

    print(f"\n总共生成 {len(amplitudes)} 个位移结构")
    print(f"下一步：将这些POSCAR文件用于VASP能带计算")


def read_eigenvector_from_phonopy(phonopy_yaml, mode_index):
    """
    从phonopy的输出文件读取声子本征矢

    参数:
        phonopy_yaml: phonopy的band.yaml或mesh.yaml文件
        mode_index: 声子模式索引（从0开始）

    返回:
        eigenvector: (natoms, 3)数组
        frequency: 声子频率（THz）
    """
    # 这里简化处理，实际使用时建议用phonopy的python API
    print("警告：从phonopy读取本征矢功能未实现")
    print("建议使用phonopy的python API:")
    print("  from phonopy import load")
    print("  phonon = load(phonopy_yaml='phonopy_disp.yaml')")
    print("  eigenvectors = phonon.get_band_structure().eigenvectors")
    return None, None


def estimate_zpe_amplitude(frequency_cm, mass_amu):
    """
    估算零点振幅（ZPE）

    公式：A_ZPE = sqrt(ℏ / (2 * M * ω))

    参数:
        frequency_cm: 声子频率（cm⁻¹）
        mass_amu: 有效质量（amu）

    返回:
        A_ZPE: 零点振幅（Å）
    """
    # 常数
    hbar = 1.054571817e-34  # J·s
    amu_to_kg = 1.66053906660e-27  # kg
    c = 2.99792458e10  # cm/s

    # 转换频率到SI单位（rad/s）
    omega_SI = 2 * np.pi * frequency_cm * c

    # 计算ZPE振幅
    M_SI = mass_amu * amu_to_kg
    A_ZPE_m = np.sqrt(hbar / (2 * M_SI * omega_SI))
    A_ZPE_angstrom = A_ZPE_m * 1e10  # 转换为Å

    return A_ZPE_angstrom


def main():
    """主程序"""
    print("=" * 70)
    print("冻结声子位移结构生成器")
    print("=" * 70)

    if len(sys.argv) < 2:
        print("\n用法：")
        print("  python generate_displaced_structures.py <POSCAR> [选项]")
        print("\n选项：")
        print("  --mode <N>           声子模式编号（默认：0）")
        print("  --amplitudes <A1,A2,...>  振幅列表（Å，逗号分隔，默认：0.0,0.05,0.1,0.15,0.2）")
        print("  --eigenvector <file>     本征矢文件（numpy格式，.npy）")
        print("  --output <dir>           输出目录（默认：displaced_structures）")
        print("\n示例：")
        print("  python generate_displaced_structures.py POSCAR --mode 25 \\")
        print("         --amplitudes 0.0,0.05,0.1,0.15,0.2 \\")
        print("         --eigenvector mode25_eigenvector.npy")
        sys.exit(1)

    # 解析参数
    poscar_file = sys.argv[1]
    mode_number = 0
    amplitudes = [0.0, 0.05, 0.1, 0.15, 0.2]
    eigenvector_file = None
    output_dir = "displaced_structures"

    i = 2
    while i < len(sys.argv):
        if sys.argv[i] == "--mode":
            mode_number = int(sys.argv[i+1])
            i += 2
        elif sys.argv[i] == "--amplitudes":
            amplitudes = [float(x) for x in sys.argv[i+1].split(',')]
            i += 2
        elif sys.argv[i] == "--eigenvector":
            eigenvector_file = sys.argv[i+1]
            i += 2
        elif sys.argv[i] == "--output":
            output_dir = sys.argv[i+1]
            i += 2
        else:
            print(f"未知选项：{sys.argv[i]}")
            i += 1

    # 检查文件
    if not os.path.exists(poscar_file):
        print(f"错误：POSCAR文件不存在 {poscar_file}")
        sys.exit(1)

    if eigenvector_file is None:
        print("错误：必须提供本征矢文件（--eigenvector）")
        print("提示：可以从phonopy或VASP的OUTCAR中提取")
        sys.exit(1)

    if not os.path.exists(eigenvector_file):
        print(f"错误：本征矢文件不存在 {eigenvector_file}")
        sys.exit(1)

    # 读取本征矢
    print(f"\n读取本征矢：{eigenvector_file}")
    eigenvector = np.load(eigenvector_file)
    print(f"  本征矢形状：{eigenvector.shape}")

    # 生成位移结构
    generate_displaced_structures(
        poscar_file,
        eigenvector,
        amplitudes,
        output_dir,
        mode_number
    )

    # ZPE估算（假设主要是H原子振动）
    print("\n" + "=" * 70)
    print("零点振幅（ZPE）估算：")
    print("=" * 70)

    freq_example = 1000  # cm⁻¹（示例值，需要从声子谱读取）
    mass_H = 1.008  # amu

    A_ZPE = estimate_zpe_amplitude(freq_example, mass_H)
    print(f"\n假设：频率 = {freq_example} cm⁻¹，质量 = {mass_H} amu （H原子）")
    print(f"零点振幅 A_ZPE = {A_ZPE:.4f} Å")
    print(f"\n如果能级交叉发生在振幅 < {A_ZPE:.4f} Å → 交叉是物理可达的")
    print(f"如果能级交叉发生在振幅 > {A_ZPE:.4f} Å → 交叉不可达（仅理论）")

    print("\n" + "=" * 70)
    print("完成！")
    print("=" * 70)


if __name__ == "__main__":
    main()
