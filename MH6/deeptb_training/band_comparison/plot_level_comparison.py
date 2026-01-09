"""
绘制纯H6与CaH6的Gamma点能级对比图
可视化能级顺序变化
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_gamma_levels(data_dir):
    """从训练数据中提取Gamma点能级"""
    eigenval_file = Path(data_dir) / 'eigenvalues.npy'
    eigenvalues = np.load(eigenval_file)
    gamma_levels = eigenvalues[0, 0, :]
    return np.sort(gamma_levels)

def plot_level_comparison(h6_levels, cah6_levels, output_file='gamma_levels_plot.png'):
    """绘制能级棒图对比"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    # 绘制纯H6能级（左侧）
    x_h6 = 0.8
    for i, e in enumerate(h6_levels):
        color = 'blue' if e < 0 else 'red'
        ax.hlines(e, x_h6-0.15, x_h6+0.15, colors=color, linewidth=3,
                 label='成键' if (i==0 and e<0) else ('反键' if (i==len(h6_levels)//2 and e>=0) else ''))

    # 绘制CaH6能级（右侧）
    x_cah6 = 2.2
    for i, e in enumerate(cah6_levels):
        if e < -1.5:
            color = 'green'
            label_text = 'Ca主导' if i == 0 else ''
        elif e < 0:
            color = 'blue'
            label_text = ''
        else:
            color = 'red'
            label_text = ''

        ax.hlines(e, x_cah6-0.15, x_cah6+0.15, colors=color, linewidth=3, label=label_text)

    # 连接相同顺序的能级（排除Ca的前2条）
    cah6_h_levels = cah6_levels[2:]
    if len(cah6_h_levels) == len(h6_levels):
        for i, (e1, e2) in enumerate(zip(h6_levels, cah6_h_levels)):
            ax.plot([x_h6+0.15, x_cah6-0.15], [e1, e2], 'k--', alpha=0.2, linewidth=0.5)

    # 费米能级参考线
    ax.axhline(0, color='black', linestyle=':', linewidth=2, label='费米能级(参考)')

    # 设置图表属性
    ax.set_xlim(0.3, 2.7)
    ax.set_ylim(cah6_levels.min()-0.5, h6_levels.max()+0.5)
    ax.set_xticks([x_h6, x_cah6])
    ax.set_xticklabels(['纯H6', 'CaH6'], fontsize=14)
    ax.set_ylabel('能量 (eV)', fontsize=14)
    ax.set_title('Gamma点能级对比：纯H6 vs CaH6', fontsize=16, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(axis='y', alpha=0.3)

    # 标注关键信息
    textstr = f'纯H6: {len(h6_levels)}条能带\nCaH6: {len(cah6_levels)}条能带\n(绿色=Ca主导)'
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ 能级对比图已保存: {output_file}")

def plot_energy_shifts(h6_levels, cah6_levels, output_file='energy_shifts.png'):
    """绘制能量偏移图"""
    cah6_h_levels = cah6_levels[2:]  # 排除Ca主导的前2条

    if len(cah6_h_levels) != len(h6_levels):
        print("⚠ 能带数量不匹配，跳过能量偏移图")
        return

    energy_shifts = cah6_h_levels - h6_levels
    band_indices = np.arange(1, len(h6_levels)+1)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.bar(band_indices, energy_shifts, color='steelblue', alpha=0.7, edgecolor='black')
    ax.axhline(0, color='red', linestyle='--', linewidth=1)

    ax.set_xlabel('能带序号', fontsize=12)
    ax.set_ylabel('能量偏移 (CaH6 - H6, eV)', fontsize=12)
    ax.set_title('Ca原子引起的能量偏移', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ 能量偏移图已保存: {output_file}")

if __name__ == "__main__":
    # 数据路径
    h6_data_dir = '../h6_pure/data/set.0'
    cah6_data_dir = '../cah6_full/data/set.0'

    print("正在加载数据...")
    h6_levels = load_gamma_levels(h6_data_dir)
    cah6_levels = load_gamma_levels(cah6_data_dir)

    print(f"纯H6: {len(h6_levels)}条能带")
    print(f"CaH6: {len(cah6_levels)}条能带")

    # 绘图
    plot_level_comparison(h6_levels, cah6_levels)
    plot_energy_shifts(h6_levels, cah6_levels)

    print("\n✓ 所有图表生成完成！")
