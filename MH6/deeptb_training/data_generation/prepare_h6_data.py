"""
生成纯H6笼子的DeePTB训练数据
使用简单的紧束缚模型生成伪数据用于验证流程
"""
import numpy as np
import json
from pathlib import Path
from ase.io import read, write
from ase import Atoms
import sys
sys.path.append('.')  # 添加当前目录到路径
from generate_tb_eigenvalues import generate_deeptb_data

def prepare_h6_dataset(output_dir: str = '../h6_pure/data/set.0'):
    """
    准备纯H6的DeePTB训练数据集

    Parameters:
    -----------
    output_dir : str
        输出目录路径
    """
    print("=" * 60)
    print("开始生成纯H6 DeePTB训练数据")
    print("=" * 60)

    # 创建输出目录
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # 输入文件路径
    poscar_file = '../../POSCAR'  # 纯H6结构

    print(f"\n1. 读取结构文件: {poscar_file}")
    structure = read(poscar_file)
    print(f"   - 晶格常数: a = {structure.cell.lengths()[0]:.4f} Å")
    print(f"   - 原子数量: {len(structure)} H atoms")

    # 紧束缚参数设置
    print("\n2. 设置紧束缚参数:")
    hopping_params = {
        'default': -2.0  # H-H跳跃积分 (eV)
    }
    onsite_energies = {
        'H': 0.0  # H的onsite能级 (eV)
    }
    cutoff = 1.3  # 截断距离 (Å)，覆盖最近邻
    n_k_per_segment = 100  # 每段k路径100个点

    print(f"   - H-H跳跃: t = {hopping_params['default']} eV")
    print(f"   - H onsite: e₀ = {onsite_energies['H']} eV")
    print(f"   - 截断距离: {cutoff} Å")
    print(f"   - k点采样: {n_k_per_segment} 点/段")

    # 生成能带数据
    print("\n3. 计算紧束缚能带...")
    eigenvalues, kpoints = generate_deeptb_data(
        structure_file=poscar_file,
        element_list=['H'],
        hopping_params=hopping_params,
        onsite_energies=onsite_energies,
        cutoff=cutoff,
        n_k_per_segment=n_k_per_segment
    )

    print(f"   - eigenvalues shape: {eigenvalues.shape}")  # (1, 401, 12)
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
            "band_max": 11,  # 12条能带，索引0-11
            "emin": None,
            "emax": None
        },
        "description": "Pure H6 cage tight-binding model",
        "parameters": {
            "t_HH": hopping_params['default'],
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
    print(f"   成键态 (负能量):")
    bonding = gamma_energies[gamma_energies < 0]
    for i, e in enumerate(bonding):
        print(f"     能带 {i+1}: {e:.4f} eV")

    print(f"\n   反键态 (正能量):")
    antibonding = gamma_energies[gamma_energies >= 0]
    for i, e in enumerate(antibonding):
        print(f"     能带 {len(bonding)+i+1}: {e:.4f} eV")

    print("\n" + "=" * 60)
    print("✓ 纯H6数据生成完成！")
    print(f"✓ 输出目录: {output_path.absolute()}")
    print("=" * 60)

    return output_path


if __name__ == "__main__":
    # 执行数据生成
    output_dir = '../h6_pure/data/set.0'
    prepare_h6_dataset(output_dir)

    print("\n下一步：")
    print("  1. 检查生成的数据文件")
    print("  2. 运行 prepare_cah6_data.py 生成CaH6数据")
    print("  3. 配置 h6_pure/input.json")
    print("  4. 开始训练: dptb train h6_pure/input.json -o h6_pure/output")
