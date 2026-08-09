#!/usr/bin/env python3
"""
分析冻结声子计算中的能级交叉和能带分裂
这是判断电声耦合强度的关键证据
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import glob


def read_eigenval_at_gamma(eigenval_file):
    """
    从EIGENVAL文件提取Γ点的能带能量

    参数:
        eigenval_file: EIGENVAL文件路径

    返回:
        energies_at_gamma: Γ点的能带能量数组 (nbands,)
    """
    with open(eigenval_file, 'r') as f:
        lines = f.readlines()

    # 解析头部
    nelect, nkpts, nbands = map(int, lines[5].split())

    # Γ点通常是第一个k点（index=0）
    # 跳过k点坐标行（第8行）
    idx = 8

    # 读取能带能量
    energies = []
    for ib in range(nbands):
        energy = float(lines[idx].split()[1])  # 第2列
        energies.append(energy)
        idx += 1

    return np.array(energies)


def extract_fermi_energy(doscar_file):
    """从DOSCAR提取费米能级"""
    with open(doscar_file, 'r') as f:
        lines = f.readlines()
    params = lines[5].split()
    efermi = float(params[3])
    return efermi


def track_bands(energies_list, threshold=0.5):
    """
    追踪能带连续性（简单版本）

    参数:
        energies_list: 列表的列表，每个元素是某个振幅下的能带能量
        threshold: 能量差阈值（eV），超过此值认为能带交叉

    返回:
        tracked_bands: (n_amplitudes, nbands)数组，每行是追踪后的能带
    """
    n_amplitudes = len(energies_list)
    nbands = len(energies_list[0])

    tracked_bands = np.zeros((n_amplitudes, nbands))
    tracked_bands[0, :] = sorted(energies_list[0])

    for i in range(1, n_amplitudes):
        current_energies = sorted(energies_list[i])

        # 简单匹配：贪心算法，找最近的能带
        used = [False] * nbands
        for j in range(nbands):
            prev_energy = tracked_bands[i-1, j]
            best_match = -1
            best_diff = float('inf')

            for k in range(nbands):
                if not used[k]:
                    diff = abs(current_energies[k] - prev_energy)
                    if diff < best_diff:
                        best_diff = diff
                        best_match = k

            tracked_bands[i, j] = current_energies[best_match]
            used[best_match] = True

            # 检查是否发生交叉
            if best_diff > threshold:
                print(f"  警告：振幅索引 {i-1}→{i}，能带 {j} 能量跳变 {best_diff:.3f} eV（可能交叉）")

    return tracked_bands


def main():
    """主程序"""
    print("=" * 80)
    print("冻结声子能级交叉分析")
    print("=" * 80)

    if len(sys.argv) < 2:
        print("\n用法：")
        print("  python analyze_band_splitting.py <模式编号> [振幅列表]")
        print("\n说明：")
        print("  - 假设目录结构：mode<N>_amp<A>/EIGENVAL")
        print("  - 例如：mode25_amp0.000/EIGENVAL, mode25_amp0.050/EIGENVAL, ...")
        print("\n示例：")
        print("  python analyze_band_splitting.py 25")
        print("  python analyze_band_splitting.py 25 0.0,0.05,0.1,0.15,0.2")
        sys.exit(1)

    mode_num = int(sys.argv[1])

    # 振幅列表
    if len(sys.argv) > 2:
        amplitudes = [float(x) for x in sys.argv[2].split(',')]
    else:
        # 自动检测
        pattern = f"mode{mode_num}_amp*/EIGENVAL"
        files = glob.glob(pattern)
        if len(files) == 0:
            print(f"错误：未找到匹配的文件 {pattern}")
            sys.exit(1)

        amplitudes = []
        for f in files:
            # 从路径提取振幅
            dirname = os.path.dirname(f)
            amp_str = dirname.split('amp')[-1]
            amplitudes.append(float(amp_str))
        amplitudes = sorted(amplitudes)

    print(f"\n模式编号：{mode_num}")
    print(f"振幅列表：{amplitudes}")

    # 提取所有振幅的能带能量
    print("\n提取能带能量...")
    all_energies = []
    efermi = None

    for amp in amplitudes:
        dirname = f"mode{mode_num}_amp{amp:.3f}"
        eigenval_file = os.path.join(dirname, "EIGENVAL")
        doscar_file = os.path.join(dirname, "DOSCAR")

        if not os.path.exists(eigenval_file):
            print(f"警告：文件不存在 {eigenval_file}，跳过")
            continue

        # 提取费米能级（从第一个计算）
        if efermi is None and os.path.exists(doscar_file):
            efermi = extract_fermi_energy(doscar_file)
            print(f"  费米能级：{efermi:.4f} eV")

        # 提取Γ点能带
        energies = read_eigenval_at_gamma(eigenval_file)
        all_energies.append(energies)
        print(f"  振幅 {amp:.3f} Å：读取 {len(energies)} 条能带")

    if len(all_energies) == 0:
        print("错误：未读取到任何数据！")
        sys.exit(1)

    # 相对费米能级
    if efermi is None:
        print("警告：未找到DOSCAR，假设费米能级=0")
        efermi = 0.0

    all_energies_rel = [(e - efermi) for e in all_energies]

    # 追踪能带
    print("\n追踪能带连续性...")
    tracked_bands = track_bands(all_energies_rel, threshold=0.3)

    # 绘图：能带演化
    print("\n绘制能带演化图...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    # === 左图：所有能带 ===
    ax1 = axes[0]
    nbands = tracked_bands.shape[1]

    # 只显示费米能附近的能带（±3 eV）
    energy_window = (-3, 3)
    colors = plt.cm.viridis(np.linspace(0, 1, nbands))

    for ib in range(nbands):
        band_energies = tracked_bands[:, ib]
        if any((e > energy_window[0]) and (e < energy_window[1]) for e in band_energies):
            ax1.plot(amplitudes, band_energies, 'o-', color=colors[ib],
                     linewidth=1.5, markersize=5, alpha=0.7)

    ax1.axhline(0, color='r', linestyle='--', linewidth=2, label='$E_F$')
    ax1.set_xlabel('Phonon Amplitude (Å)', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Energy - $E_F$ (eV)', fontsize=14, fontweight='bold')
    ax1.set_ylim(energy_window)
    ax1.set_title(f'Mode {mode_num}: Band Evolution (All bands)', fontsize=16, fontweight='bold')
    ax1.legend(fontsize=12)
    ax1.grid(alpha=0.4)

    # === 右图：费米能附近放大（±0.5 eV）===
    ax2 = axes[1]
    energy_window_zoom = (-0.5, 0.5)

    for ib in range(nbands):
        band_energies = tracked_bands[:, ib]
        if any((e > energy_window_zoom[0]) and (e < energy_window_zoom[1]) for e in band_energies):
            ax2.plot(amplitudes, band_energies, 'o-', color=colors[ib],
                     linewidth=2, markersize=6, alpha=0.8, label=f'Band {ib+1}')

    ax2.axhline(0, color='r', linestyle='--', linewidth=2, label='$E_F$')
    ax2.set_xlabel('Phonon Amplitude (Å)', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Energy - $E_F$ (eV)', fontsize=14, fontweight='bold')
    ax2.set_ylim(energy_window_zoom)
    ax2.set_title(f'Mode {mode_num}: Near Fermi Level (Zoom)', fontsize=16, fontweight='bold')
    ax2.grid(alpha=0.4)

    plt.tight_layout()
    output_file = f'mode{mode_num}_band_evolution.png'
    plt.savefig(output_file, dpi=300)
    print(f"  图像已保存：{output_file}")

    # === 分析能级交叉 ===
    print("\n" + "=" * 80)
    print("能级交叉分析：")
    print("=" * 80)

    crossings_detected = []
    for i in range(len(amplitudes) - 1):
        amp1, amp2 = amplitudes[i], amplitudes[i+1]

        # 检查能带顺序是否变化
        order1 = np.argsort(tracked_bands[i, :])
        order2 = np.argsort(tracked_bands[i+1, :])

        if not np.array_equal(order1, order2):
            n_changes = np.sum(order1 != order2)
            print(f"\n振幅 {amp1:.3f} → {amp2:.3f} Å：")
            print(f"  检测到 {n_changes} 处能带顺序变化（可能的能级交叉）")

            # 找出具体哪些能带交叉
            for j in range(len(order1)):
                if order1[j] != order2[j]:
                    E1 = tracked_bands[i, j]
                    E2 = tracked_bands[i+1, j]
                    print(f"    能带 {j+1}：{E1:.4f} eV → {E2:.4f} eV")

            crossings_detected.append((amp1, amp2))

    # === 简并破缺分析 ===
    print("\n" + "=" * 80)
    print("简并破缺分析：")
    print("=" * 80)

    # 在平衡位置（振幅=0）检查简并
    idx_eq = amplitudes.index(0.0) if 0.0 in amplitudes else 0
    energies_eq = tracked_bands[idx_eq, :]

    # 寻找简并态（能量差<10 meV）
    degeneracy_threshold = 0.01  # eV
    degenerate_groups = []
    visited = [False] * nbands

    for i in range(nbands):
        if visited[i]:
            continue
        group = [i]
        visited[i] = True
        for j in range(i+1, nbands):
            if abs(energies_eq[i] - energies_eq[j]) < degeneracy_threshold:
                group.append(j)
                visited[j] = True
        if len(group) > 1:
            degenerate_groups.append(group)

    if len(degenerate_groups) > 0:
        print(f"\n平衡位置（振幅={amplitudes[idx_eq]:.3f} Å）发现 {len(degenerate_groups)} 组简并态：")
        for ig, group in enumerate(degenerate_groups):
            E_avg = np.mean(energies_eq[group])
            print(f"  组 {ig+1}：能带 {[g+1 for g in group]}，平均能量 {E_avg:.4f} eV")

            # 检查简并是否在后续振幅中破缺
            splitting = []
            for i in range(len(amplitudes)):
                energies_at_amp = tracked_bands[i, group]
                E_max = np.max(energies_at_amp)
                E_min = np.min(energies_at_amp)
                splitting.append(E_max - E_min)

            # 找到分裂最大的振幅
            max_split = np.max(splitting)
            max_split_amp = amplitudes[np.argmax(splitting)]

            print(f"    最大分裂：{max_split*1000:.1f} meV @ 振幅 {max_split_amp:.3f} Å")

            if max_split > 0.05:  # 50 meV
                print(f"    → 简并显著破缺！强电声耦合的标志")
            else:
                print(f"    → 简并保持稳定，弱电声耦合")

    # === ZPE判据 ===
    print("\n" + "=" * 80)
    print("零点振幅（ZPE）判据：")
    print("=" * 80)

    # 假设主要是H原子振动
    freq_example = 1000  # cm⁻¹（需要从声子谱读取）
    mass_H = 1.008  # amu

    hbar = 1.054571817e-34  # J·s
    amu_to_kg = 1.66053906660e-27
    c = 2.99792458e10  # cm/s
    omega_SI = 2 * np.pi * freq_example * c
    M_SI = mass_H * amu_to_kg
    A_ZPE = np.sqrt(hbar / (2 * M_SI * omega_SI)) * 1e10  # Å

    print(f"\n假设：频率 = {freq_example} cm⁻¹，质量 = {mass_H} amu")
    print(f"零点振幅 A_ZPE ≈ {A_ZPE:.4f} Å")

    if len(crossings_detected) > 0:
        first_crossing_amp = crossings_detected[0][0]
        print(f"\n第一个能级交叉发生在 ≈ {first_crossing_amp:.3f} Å")

        if first_crossing_amp <= A_ZPE:
            print(f"  → {first_crossing_amp:.3f} Å < {A_ZPE:.4f} Å：交叉在ZPE范围内")
            print(f"  → 结论：强电声耦合（类似YB6）")
        else:
            print(f"  → {first_crossing_amp:.3f} Å > {A_ZPE:.4f} Å：交叉超出ZPE")
            print(f"  → 结论：弱电声耦合（类似LaB6）")
    else:
        print("\n未检测到明显的能级交叉")
        print(f"  → 结论：弱电声耦合，简并态稳定")

    # 保存结果
    output_summary = f"mode{mode_num}_crossing_summary.txt"
    with open(output_summary, 'w') as f:
        f.write(f"模式编号：{mode_num}\n")
        f.write(f"振幅范围：{min(amplitudes):.3f} - {max(amplitudes):.3f} Å\n")
        f.write(f"零点振幅（估算）：{A_ZPE:.4f} Å\n")
        f.write(f"检测到的能级交叉：{len(crossings_detected)}\n")
        if len(crossings_detected) > 0:
            f.write(f"第一个交叉振幅：{crossings_detected[0][0]:.3f} Å\n")
            if crossings_detected[0][0] <= A_ZPE:
                f.write(f"结论：强电声耦合（交叉在ZPE内）\n")
            else:
                f.write(f"结论：弱电声耦合（交叉超出ZPE）\n")
        else:
            f.write(f"结论：弱电声耦合（无交叉）\n")

    print(f"\n结果已保存：{output_summary}")

    print("\n" + "=" * 80)
    print("分析完成！")
    print("=" * 80)


if __name__ == "__main__":
    main()
