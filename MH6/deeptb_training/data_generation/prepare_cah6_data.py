"""
生成CaH6的DeePTB训练数据
包含Ca和H两种元素，模拟Ca电子填充效应
"""
import numpy as np
import json
from pathlib import Path
from ase.io import read, write
import sys
sys.path.append('.')
from generate_tb_eigenvalues import generate_deeptb_data

def prepare_cah6_dataset(output_dir: str = '../cah6_full/data/set.0'):
    """
    准备CaH6的DeePTB训练数据集

    Parameters:
    -----------
    output_dir : str
        输出目录路径
    """
    print("=" * 60)
    print("开始生成CaH6 DeePTB训练数据")
    print("=" * 60)

    # 创建输出目录
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # 输入文件路径
    poscar_file = '../../std_Ca2H12_Im-3m_229_.vasp'  # CaH6完整结构

    print(f"\n1. 读取结构文件: {poscar_file}")
    structure = read(poscar_file)
    print(f"   - 晶格常数: a = {structure.cell.lengths()[0]:.4f} Å")
    print(f"   - 原子组成: {structure.get_chemical_formula()}")
    print(f"   - 总原子数: {len(structure)}")

    # 紧束缚参数设置
    print("\n2. 设置紧束缚参数:")
    # 这里需要扩展hopping_params来区分Ca-H和H-H键
    # 简化版本：统一跳跃参数
    hopping_params = {
        'default': -1.5  # 简化的平均跳跃积分
    }
    onsite_energies = {
        'Ca': -2.0,  # Ca的onsite能级，模拟价电子填充
        'H': 0.0     # H的onsite能级
    }
    cutoff = 1.27  # 增大截断距离以包含Ca-H作用

    print(f"   - 跳跃积分: t = {hopping_params['default']} eV (平均值)")
    print(f"   - Ca onsite: e_Ca = {onsite_energies['Ca']} eV")
    print(f"   - H onsite: e_H = {onsite_energies['H']} eV")
    print(f"   - 截断距离: {cutoff} Å")

    # 生成能带数据
    print("\n3. 计算紧束缚能带...")
    eigenvalues, kpoints = generate_deeptb_data(
        structure_file=poscar_file,
        element_list=['Ca', 'H'],
        hopping_params=hopping_params,
        onsite_energies=onsite_energies,
        cutoff=cutoff,
        n_k_per_segment=100
    )

    print(f"   - eigenvalues shape: {eigenvalues.shape}")  # (1, 401, 14)
    print(f"   - kpoints shape: {kpoints.shape}")          # (401, 3)
    print(f"   - 能带范围: {eigenvalues.min():.3f} ~ {eigenvalues.max():.3f} eV")

    # 保存numpy数据
    print("\n4. 保存numpy数据文件...")
    np.save(output_path / 'eigenvalues.npy', eigenvalues)
    np.save(output_path / 'kpoints.npy', kpoints)
    print(f"   ✓ eigenvalues.npy")
    print(f"   ✓ kpoints.npy")

    # 保存ASE轨迹文件
    print("\n5. 保存ASE轨迹文件...")
    write(output_path / 'xdat.traj', structure)
    print(f"   ✓ xdat.traj")

    # 创建info.json
    print("\n6. 创建info.json元数据...")
    info = {
        "nframes": 1,
        "pos_type": "ase",
        "pbc": True,
        "bandinfo": {
            "band_min": 0,
            "band_max": 13,  # 14条能带，索引0-13
            "emin": None,
            "emax": None
        },
        "description": "CaH6 with Ca and H cage tight-binding model",
        "parameters": {
            "t_avg": hopping_params['default'],
            "e_Ca": onsite_energies['Ca'],
            "e_H": onsite_energies['H'],
            "cutoff": cutoff
        }
    }

    with open(output_path / 'info.json', 'w') as f:
        json.dump(info, f, indent=4)
    print(f"   ✓ info.json")

    # 分析Gamma点能级
    print("\n7. Gamma点能级分析:")
    gamma_energies = eigenvalues[0, 0, :]  # 第一个k点是Gamma点
    sorted_energies = np.sort(gamma_energies)

    print(f"   所有能级 (按能量排序):")
    for i, e in enumerate(sorted_energies):
        if e < -1.5:
            label = "Ca-dominated"
        elif e < 0:
            label = "H成键"
        else:
            label = "H反键"
        print(f"     能带 {i+1}: {e:.4f} eV  ({label})")

    print(f"\n   费米能级附近 (假设Ca贡献2个电子):")
    print(f"   → 前2条能带被占据 (Ca的2个价电子)")
    print(f"   → 费米能级预期位于能带2和能带3之间")

    print("\n" + "=" * 60)
    print("✓ CaH6数据生成完成！")
    print(f"✓ 输出目录: {output_path.absolute()}")
    print("=" * 60)

    return output_path


if __name__ == "__main__":
    # 执行数据生成
    output_dir = '../cah6_full/data/set.0'
    prepare_cah6_dataset(output_dir)

    print("\n下一步：")
    print("  1. 对比H6和CaH6的Gamma点能级")
    print("  2. 配置 cah6_full/input.json")
    print("  3. 开始训练: dptb train cah6_full/input.json -o cah6_full/output")
    print("  4. 运行分析脚本对比两个体系")
