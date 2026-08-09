"""
提取和分析纯H6与CaH6的Gamma点能级
检测能级顺序变化和可能的轨道反转
"""
import numpy as np
from pathlib import Path

def load_gamma_levels(data_dir):
    """从训练数据中提取Gamma点能级"""
    eigenval_file = Path(data_dir) / 'eigenvalues.npy'
    eigenvalues = np.load(eigenval_file)
    # shape: (nframes, nkpoints, nbands)
    # Gamma点是第一个k点（索引0）
    gamma_levels = eigenvalues[0, 0, :]
    return gamma_levels

def analyze_level_changes(h6_levels, cah6_levels):
    """分析能级顺序变化"""
    print("\n" + "="*70)
    print("Gamma点能级顺序分析")
    print("="*70)

    # 纯H6能级分析
    print("\n【纯H6笼子】")
    print(f"能带数量: {len(h6_levels)}")
    print(f"能量范围: {h6_levels.min():.4f} ~ {h6_levels.max():.4f} eV")
    print("\n能级列表:")
    for i, e in enumerate(np.sort(h6_levels)):
        char = "成键" if e < 0 else "反键"
        print(f"  能带 {i+1:2d}: {e:8.4f} eV  ({char})")

    # CaH6能级分析
    print("\n【CaH6完整体系】")
    print(f"能带数量: {len(cah6_levels)}")
    print(f"能量范围: {cah6_levels.min():.4f} ~ {cah6_levels.max():.4f} eV")
    print("\n能级列表:")
    sorted_cah6 = np.sort(cah6_levels)
    for i, e in enumerate(sorted_cah6):
        if e < -1.5:
            char = "Ca主导"
        elif e < 0:
            char = "H成键"
        else:
            char = "H反键"
        print(f"  能带 {i+1:2d}: {e:8.4f} eV  ({char})")

    # 能级反转检测
    print("\n" + "-"*70)
    print("能级顺序变化检测")
    print("-"*70)

    # 比较H6和CaH6中H相关的能级（排除Ca的最低2条能带）
    h6_sorted = np.sort(h6_levels)
    cah6_h_levels = sorted_cah6[2:]  # 假设前2条是Ca主导

    print(f"\n比较策略：")
    print(f"  - 纯H6: 全部{len(h6_sorted)}条能带")
    print(f"  - CaH6: 排除Ca主导的前2条，剩余{len(cah6_h_levels)}条H相关能带")

    # 检查能量差异
    print(f"\n能量偏移分析:")
    if len(cah6_h_levels) == len(h6_sorted):
        energy_shifts = cah6_h_levels - h6_sorted
        avg_shift = np.mean(energy_shifts)
        max_shift = np.max(np.abs(energy_shifts))

        print(f"  平均能量偏移: {avg_shift:+.4f} eV")
        print(f"  最大能量偏移: {max_shift:.4f} eV")

        # 检测反转（符号变化）
        sign_changes = []
        for i in range(len(energy_shifts)-1):
            if energy_shifts[i] * energy_shifts[i+1] < 0:
                sign_changes.append((i, i+1))

        if sign_changes:
            print(f"\n⚠ 检测到{len(sign_changes)}处能级顺序可能反转:")
            for i, j in sign_changes:
                print(f"    能带 {i+1} 与 {j+1} 之间")
                print(f"      H6:  {h6_sorted[i]:.4f} eV vs {h6_sorted[j]:.4f} eV")
                print(f"      CaH6: {cah6_h_levels[i]:.4f} eV vs {cah6_h_levels[j]:.4f} eV")
        else:
            print(f"\n✓ 未检测到明显的能级顺序反转")
            print(f"  → Ca主要影响能量偏移，未改变轨道顺序")

    # 费米面分析
    print(f"\n" + "-"*70)
    print("费米能级附近分析")
    print("-"*70)

    # 纯H6（假设费米能在0 eV）
    h6_near_ef = h6_sorted[(h6_sorted > -1) & (h6_sorted < 1)]
    print(f"\n纯H6费米能附近能级 (±1 eV):")
    for e in h6_near_ef:
        print(f"  {e:8.4f} eV")

    # CaH6（考虑2个Ca电子填充）
    print(f"\nCaH6费米能附近 (考虑Ca的2个价电子):")
    print(f"  → 前2条能带被Ca电子占据")
    print(f"  → 实际费米能级在第3条能带附近")
    cah6_fermi_region = sorted_cah6[1:5]  # 第2-5条能带
    for i, e in enumerate(cah6_fermi_region, start=2):
        print(f"  能带{i}: {e:8.4f} eV")

    print("\n" + "="*70)

    return {
        'h6_levels': h6_sorted,
        'cah6_levels': sorted_cah6,
        'energy_shifts': energy_shifts if 'energy_shifts' in locals() else None
    }

def save_results(results, output_file='gamma_levels_comparison.txt'):
    """保存分析结果到文件"""
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("Gamma点能级对比分析结果\n")
        f.write("="*70 + "\n\n")

        f.write("纯H6能级:\n")
        for i, e in enumerate(results['h6_levels']):
            f.write(f"  能带{i+1:2d}: {e:8.4f} eV\n")

        f.write("\nCaH6能级:\n")
        for i, e in enumerate(results['cah6_levels']):
            f.write(f"  能带{i+1:2d}: {e:8.4f} eV\n")

        if results['energy_shifts'] is not None:
            f.write("\n能量偏移:\n")
            for i, shift in enumerate(results['energy_shifts']):
                f.write(f"  能带{i+1:2d}: {shift:+.4f} eV\n")

    print(f"\n✓ 结果已保存到: {output_file}")

if __name__ == "__main__":
    # 数据路径
    h6_data_dir = '../h6_pure/data/set.0'
    cah6_data_dir = '../cah6_full/data/set.0'

    print("正在加载数据...")
    h6_levels = load_gamma_levels(h6_data_dir)
    cah6_levels = load_gamma_levels(cah6_data_dir)

    # 分析
    results = analyze_level_changes(h6_levels, cah6_levels)

    # 保存结果
    save_results(results)

    print("\n分析完成！")
    print("下一步：运行 plot_level_comparison.py 绘制对比图")
