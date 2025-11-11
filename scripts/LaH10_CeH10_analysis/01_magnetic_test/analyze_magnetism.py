#!/usr/bin/env python3
"""
分析磁性计算结果
从VASP的OUTCAR文件提取能量和磁矩
"""
import os
import re
import sys

def parse_outcar(filename):
    """
    提取OUTCAR中的总能量和原子磁矩

    返回:
        energy: 总能量（eV）
        magnetization: 每个原子的磁矩列表（μB）
    """
    if not os.path.exists(filename):
        print(f"错误：文件不存在 {filename}")
        return None, None

    with open(filename, 'r') as f:
        lines = f.readlines()

    # 提取总能量（取最后一次迭代）
    energy = None
    for line in reversed(lines):
        if 'free  energy   TOTEN' in line:
            energy = float(line.split()[-2])
            break

    if energy is None:
        print(f"警告：未找到能量信息 in {filename}")

    # 提取磁矩（magnetization section）
    magnetization = []
    in_mag_section = False
    for line in lines:
        if 'magnetization (x)' in line:
            in_mag_section = True
            continue
        if in_mag_section:
            if '---' in line or 'tot' in line:
                break
            parts = line.split()
            if len(parts) >= 5:
                try:
                    mag = float(parts[4])
                    magnetization.append(mag)
                except ValueError:
                    pass

    return energy, magnetization


def main():
    """主程序"""
    print("=" * 70)
    print("Ce-based Hydride: 磁性计算结果分析")
    print("=" * 70)

    # 检查文件
    files_to_check = {
        "非磁性": "OUTCAR_NM",
        "铁磁": "OUTCAR_FM"
    }

    results = {}
    for label, filename in files_to_check.items():
        E, mag = parse_outcar(filename)
        if E is not None:
            results[label] = {"energy": E, "magnetization": mag}
        else:
            print(f"错误：无法解析 {filename}")
            return

    # 输出结果
    print(f"\n非磁性计算：")
    print(f"  能量: {results['非磁性']['energy']:.6f} eV")
    print(f"  原子数: {len(results['非磁性']['magnetization'])}")

    print(f"\n铁磁计算：")
    print(f"  能量: {results['铁磁']['energy']:.6f} eV")
    print(f"  原子数: {len(results['铁磁']['magnetization'])}")

    # Ce的磁矩（假设Ce是第一个原子）
    if len(results['铁磁']['magnetization']) > 0:
        Ce_mag = results['铁磁']['magnetization'][0]
        print(f"  Ce磁矩: {Ce_mag:.4f} μB")

        # 总磁矩
        total_mag = sum(results['铁磁']['magnetization'])
        print(f"  总磁矩: {total_mag:.4f} μB")

    # 能量差
    dE = results['铁磁']['energy'] - results['非磁性']['energy']
    print(f"\n能量差 (E_FM - E_NM): {dE:.6f} eV = {dE*1000:.2f} meV")

    # 判断与建议
    print("\n" + "=" * 70)
    print("分析与建议：")
    print("=" * 70)

    Ce_mag = results['铁磁']['magnetization'][0] if len(results['铁磁']['magnetization']) > 0 else 0.0

    if abs(Ce_mag) > 0.5 and dE < -1e-3:
        print(f"\n✓ 结论：Ce有显著磁矩 ({Ce_mag:.2f} μB)，铁磁态能量更低")
        print("  → Ce的4f¹电子局域，表现出磁性")
        print("  → 建议：")
        print("     1. 后续计算使用自旋极化DFT（ISPIN=2）")
        print("     2. 考虑测试DFT+U方法（U=3-7 eV）")
        print("     3. 测试反铁磁构型（可能能量更低）")
        print("     4. 注意磁性与超导的竞争关系")

    elif abs(Ce_mag) < 0.1:
        print(f"\n✓ 结论：Ce磁矩几乎消失 ({Ce_mag:.3f} μB)")
        print("  → Ce的4f电子在高压下被离域，无磁矩")
        print("  → 建议：")
        print("     1. 后续计算可以使用非磁性DFT（ISPIN=1）")
        print("     2. 标准GGA泛函（PBE）应该足够")
        print("     3. 无需使用DFT+U")

    else:
        print(f"\n⚠ 结论：Ce有中等磁矩 ({Ce_mag:.2f} μB)，但能量差很小")
        print("  → 系统处于磁性与非磁性的临界状态")
        print("  → 建议：")
        print("     1. 测试不同Hubbard U值（U=0, 3, 5, 7）")
        print("     2. 测试反铁磁构型")
        print("     3. 检查SOC效应（f电子重元素可能显著）")
        print("     4. 谨慎解释结果（DFT在临界区域不准确）")

    # 额外检查：能量收敛
    if abs(dE) < 1e-5:
        print("\n⚠ 警告：能量差极小（<0.01 meV），可能是：")
        print("   - SCF未完全收敛")
        print("   - K点网格不够密")
        print("   - 初始磁矩设置不当")
        print("   建议检查OUTCAR中的收敛标志")

    # 保存结果到文件
    with open("magnetic_analysis_summary.txt", 'w') as f:
        f.write(f"非磁性能量: {results['非磁性']['energy']:.6f} eV\n")
        f.write(f"铁磁能量: {results['铁磁']['energy']:.6f} eV\n")
        f.write(f"能量差: {dE:.6f} eV\n")
        f.write(f"Ce磁矩: {Ce_mag:.4f} μB\n")
        f.write(f"总磁矩: {total_mag:.4f} μB\n")

    print(f"\n结果已保存到：magnetic_analysis_summary.txt")
    print("=" * 70)


if __name__ == "__main__":
    main()
