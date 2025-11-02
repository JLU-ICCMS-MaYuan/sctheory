#!/usr/bin/env python3
"""
分析声子模式对电声耦合的贡献
"""
import numpy as np
import matplotlib.pyplot as plt

# 常数
RYDBERG_TO_CM = 219474.63  # 1 Ry = 219474.63 cm⁻¹
RYDBERG_TO_EV = 13.6057

def read_a2F_dos(filename):
    """读取a2F.dosX文件"""
    # 跳过注释行（#开头）和末尾的总结行
    # 使用comments参数跳过#开头的行
    # 使用invalid_raise=False让numpy忽略格式不一致的行
    try:
        data = np.loadtxt(filename, comments='#')
    except ValueError:
        # 如果还是报错（因为末尾的lambda行），手动处理
        with open(filename, 'r') as f:
            lines = f.readlines()
        # 只保留数据行（非注释，非空行，包含科学计数法的行）
        data_lines = []
        for line in lines:
            line = line.strip()
            if line and not line.startswith('#') and 'E' in line and '=' not in line:
                data_lines.append(line)
        # 写入临时数组
        import io
        data = np.loadtxt(io.StringIO('\n'.join(data_lines)))

    omega_ry = data[:, 0]  # Rydberg
    omega_cm = omega_ry * RYDBERG_TO_CM  # cm⁻¹
    a2F_total = data[:, 1]
    # 注意：后面的列是各个模式的贡献
    return omega_cm, a2F_total, data

def calculate_lambda_mode(omega, a2F_mode):
    """
    计算单个模式的λ_ν
    λ_ν = 2 ∫ α²F_ν(ω)/ω dω
    """
    # 排除ω=0的点（避免除零）
    mask = omega > 1e-6
    omega_masked = omega[mask]
    a2F_masked = a2F_mode[mask]

    # 数值积分（梯形法则）
    integrand = a2F_masked / omega_masked
    lambda_mode = 2.0 * np.trapz(integrand, omega_masked)

    return lambda_mode

# 主程序
if __name__ == "__main__":
    # 设置参数
    n_modes = 33  # LaH10: 11原子×3
    base_dir = "./"

    # 存储结果
    lambda_modes = []

    # 读取第一个文件获取频率网格
    omega_cm, a2F_total, data = read_a2F_dos(f"{base_dir}/a2F.dos1")

    # 对每个模式计算λ_ν
    for i in range(1, n_modes + 1):
        filename = f"{base_dir}/a2F.dos{i}"
        try:
            omega, a2F, data = read_a2F_dos(filename)
            # 第3列开始是各模式的贡献，第i+1列是模式i
            a2F_mode = data[:, i + 1]
            lambda_mode = calculate_lambda_mode(omega, a2F_mode)
            lambda_modes.append(lambda_mode)
            print(f"Mode {i:2d}: λ = {lambda_mode:.4f}")
        except Exception as e:
            print(f"Error reading mode {i}: {e}")
            lambda_modes.append(0.0)

    # 统计结果
    lambda_modes = np.array(lambda_modes)
    lambda_total = np.sum(lambda_modes)

    print(f"\n{'='*50}")
    print(f"Total λ from sum: {lambda_total:.4f}")
    print(f"{'='*50}\n")

    # 找出贡献最大的模式
    sorted_indices = np.argsort(lambda_modes)[::-1]  # 降序
    print("Top 10 contributing modes:")
    print(f"{'Rank':<6} {'Mode':<6} {'λ_ν':<10} {'Contribution %':<15}")
    print("-" * 50)
    for rank, idx in enumerate(sorted_indices[:10], 1):
        mode_num = idx + 1
        contribution_pct = 100 * lambda_modes[idx] / lambda_total
        print(f"{rank:<6} {mode_num:<6} {lambda_modes[idx]:<10.4f} {contribution_pct:<15.2f}")

    # 可视化
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # 左图：各模式的λ贡献
    ax1.bar(range(1, n_modes + 1), lambda_modes, alpha=0.7)
    ax1.set_xlabel('Phonon Mode Index', fontsize=12)
    ax1.set_ylabel('$\\lambda_\\nu$', fontsize=12)
    ax1.set_title('Contribution of Each Phonon Mode', fontsize=14)
    ax1.grid(alpha=0.3)

    # 右图：累积贡献
    cumulative = np.cumsum(lambda_modes[sorted_indices])
    ax2.plot(range(1, n_modes + 1), 100 * cumulative / lambda_total, 'o-')
    ax2.axhline(90, color='r', linestyle='--', label='90% threshold')
    ax2.set_xlabel('Number of Top Modes', fontsize=12)
    ax2.set_ylabel('Cumulative Contribution (%)', fontsize=12)
    ax2.set_title('Cumulative $\\lambda$ Contribution', fontsize=14)
    ax2.legend()
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('phonon_mode_analysis.png', dpi=300)
    print(f"\nPlot saved: phonon_mode_analysis.png")
