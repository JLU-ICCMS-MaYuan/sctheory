# LaH10与CeH10超导温度差异的化学键分析研究计划

---

## 研究背景与动机

### 核心科学问题

**为什么CeH10的超导温度比LaH10低近200K？**

| 材料 | 压力 | Tc (K) | λ | 电子构型 |
|------|------|--------|---|----------|
| LaH10 | 250 GPa | ~250 | 2.3 | La: [Xe]5d¹6s² (4f⁰) |
| CeH10 | 预测 | 50-100? | ? | Ce: [Xe]4f¹5d¹6s² (4f¹) |

**关键差异**：Ce比La多了一个4f电子

**传统解释的不足**：
- 简单的"电子数差异"无法解释如此巨大的Tc差异
- 需要从化学键和电声耦合的微观机制入手

---

### 研究思路来源：YB6/LaB6的启示

**参考文献**：Robert et al., PNAS 2024
*"Chemical bonding determines two seemingly identical superconductors' critical temperature difference"*

**YB6/LaB6对比**：
- **晶体结构**：完全相同（CaB6型，Pm-3m）
- **价电子构型**：Y ([Kr]4d¹5s²) vs La ([Xe]5d¹6s²)
- **Tc差异**：8.4 K vs 0.45 K（相差18倍！）

**核心发现**：
1. **能带反转**：La的低能4f轨道与B的π系统共价，将π带拉低到费米能；Y无可用4f轨道，σ带留在费米能
2. **耦合机制差异**：
   - **YB6**：费米能处是B-B σ键 → 对T2g摆动声子敏感 → 强EPC
   - **LaB6**：费米能处是B-B π键 → 对摆动声子钝化 → 弱EPC
3. **冻结声子证据**：YB6在~2°声子位移时发生能级交叉，LaB6即使5°也不交叉

**应用到LaH10/CeH10的假说**：

**假说1：轨道反转**
- LaH10：La-5d + H-1s形成σ键在费米能
- CeH10：Ce-4f¹参与成键，可能将π系统拉低，发生能带反转

**假说2：耦合机制差异**
- LaH10：σ键 ⊥ H振动 → 强EPC → 高Tc
- CeH10：π键 ∥ H振动 → 弱EPC → 低Tc

**假说3：磁性竞争**
- Ce的4f¹可能有局域磁矩 → 与超导竞争

---

## 研究计划总览

### 五大研究阶段

| 阶段 | 时间 | 目标 | 关键方法 |
|------|------|------|----------|
| **一** | 2-3周 | 确定Ce的4f电子性质 | 磁性、DFT+U、SOC测试 |
| **二** | 3-4周 | 识别前线轨道和能带反转 | 能带、PDOS、轨道可视化、COHP |
| **三** | 4-5周 | 确定主导声子和EPC机制 | 冻结声子、能级交叉分析 |
| **四** | 3-4周 | 超胞与费米面嵌套 | 2×2×2超胞、M点分析 |
| **五** | 2-3周 | 压力和应变效应（可选） | 变压力计算 |

**总时长**：16-19周（约4-5个月）

---

## 阶段一：4f电子性质诊断（2-3周）

### 目标

确定Ce的4f¹电子是**局域磁矩**还是**巡游电子**，这决定后续计算方案的选择。

### 物理背景

**Ce的4f电子问题**：
- **局域情况**（Ce³⁺）：4f¹强关联，可能有磁矩，需要DFT+U或DMFT
- **巡游情况**（Ce⁴⁺）：4f被压力离域，可以用标准DFT
- **中间情况**（价态波动）：最复杂，需要高级方法

**高压下的特殊性**：
- 250 GPa下，原子间距极小（~1.5 Å）
- 4f轨道可能被"压扁"而离域
- 但仍需实验验证

### 计算方案

#### 1.1 非磁性 vs 自旋极化计算

**目的**：检查Ce是否自发形成磁矩

**VASP输入文件**：`01_magnetic_test/INCAR_NM`（非磁性）
```fortran
# 非磁性计算
SYSTEM = CeH10_nonmagnetic
ISTART = 0
ICHARG = 2
ISPIN = 1          # 不考虑自旋
MAGMOM = 33*0.0    # 初始磁矩为0

# 基本设置
ENCUT = 600        # 能量截断（eV）
EDIFF = 1E-6       # 收敛标准
ISMEAR = -5        # 四面体方法+Blöchl修正
LORBIT = 11        # 输出PDOS

# 精确计算
PREC = Accurate
LREAL = .FALSE.
ADDGRID = .TRUE.
```

**VASP输入文件**：`01_magnetic_test/INCAR_FM`（铁磁）
```fortran
# 铁磁计算
SYSTEM = CeH10_ferromagnetic
ISTART = 0
ICHARG = 2
ISPIN = 2          # 考虑自旋
MAGMOM = 1*5.0 32*0.0  # Ce初始磁矩5μB，H为0

# 其他参数同上
...
```

**分析脚本**：`01_magnetic_test/analyze_magnetism.py`
```python
#!/usr/bin/env python3
"""
分析磁性计算结果
"""
import os
import re

def parse_outcar(filename):
    """提取OUTCAR中的磁矩和能量"""
    with open(filename, 'r') as f:
        lines = f.readlines()

    # 提取总能量
    for line in reversed(lines):
        if 'free  energy' in line:
            energy = float(line.split()[-2])
            break

    # 提取磁矩
    magnetization = []
    in_mag_section = False
    for line in lines:
        if 'magnetization (x)' in line:
            in_mag_section = True
            continue
        if in_mag_section and '---' in line:
            break
        if in_mag_section and len(line.split()) >= 5:
            try:
                mag = float(line.split()[4])
                magnetization.append(mag)
            except:
                pass

    return energy, magnetization

# 主程序
print("="*60)
print("磁性计算结果分析")
print("="*60)

# 非磁性结果
E_NM, mag_NM = parse_outcar("OUTCAR_NM")
print(f"\n非磁性计算：")
print(f"  能量: {E_NM:.6f} eV")

# 铁磁结果
E_FM, mag_FM = parse_outcar("OUTCAR_FM")
print(f"\n铁磁计算：")
print(f"  能量: {E_FM:.6f} eV")
print(f"  Ce磁矩: {mag_FM[0]:.3f} μB")
print(f"  总磁矩: {sum(mag_FM):.3f} μB")

# 能量差
dE = E_FM - E_NM
print(f"\n能量差 (E_FM - E_NM): {dE*1000:.2f} meV")

# 判断
if abs(mag_FM[0]) > 0.1 and dE < -1e-3:
    print("\n结论：Ce有显著磁矩且铁磁态能量更低")
    print("      → 需要考虑磁性，使用DFT+U或自旋极化计算")
elif abs(mag_FM[0]) < 0.1:
    print("\n结论：Ce磁矩几乎消失")
    print("      → 4f电子被压力离域，可用标准DFT")
else:
    print("\n结论：磁矩存在但能量差很小")
    print("      → 需要进一步测试（反铁磁、DFT+U等）")
```

#### 1.2 Hubbard U扫描

**目的**：测试4f强关联效应的影响

**VASP输入文件**：`02_hubbard_U/INCAR_U3`
```fortran
# DFT+U计算（U=3eV示例）
SYSTEM = CeH10_U3
LDAU = .TRUE.
LDAUTYPE = 2       # Dudarev方法
LDAUL = 3 -1       # Ce用f轨道，H不用
LDAUU = 3.0 0.0    # U值（eV）
LDAUJ = 0.0 0.0
LMAXMIX = 6        # f轨道需要

# 其他参数同上
...
```

**扫描参数**：U = 0, 3, 5, 7 eV（四个值）

**分析脚本**：`02_hubbard_U/compare_U_values.py`
```python
#!/usr/bin/env python3
"""
对比不同U值的计算结果
"""
import numpy as np
import matplotlib.pyplot as plt

# 手动输入或从文件读取
U_values = [0, 3, 5, 7]
energies = []  # 从各个OUTCAR提取
bandgaps = []  # 费米能附近的能隙（如果有）
dos_at_ef = [] # 费米能处的态密度

# 读取数据（省略具体代码）
...

# 绘图
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

axes[0].plot(U_values, energies, 'o-')
axes[0].set_xlabel('U (eV)')
axes[0].set_ylabel('Total Energy (eV)')
axes[0].set_title('Energy vs U')

axes[1].plot(U_values, dos_at_ef, 'o-')
axes[1].set_xlabel('U (eV)')
axes[1].set_ylabel('N(EF) (states/eV)')
axes[1].set_title('DOS at Fermi Level')

axes[2].plot(U_values, bandgaps, 'o-')
axes[2].set_xlabel('U (eV)')
axes[2].set_ylabel('Band Gap (eV)')
axes[2].set_title('Gap at Fermi Level')

plt.tight_layout()
plt.savefig('U_scan_results.png', dpi=300)
print("结果保存至 U_scan_results.png")
```

#### 1.3 自旋轨道耦合（SOC）测试

**目的**：检查SOC对能带的影响（f电子重元素可能显著）

**VASP输入文件**：`03_SOC_test/INCAR_SOC`
```fortran
# SOC计算
SYSTEM = CeH10_SOC
LSORBIT = .TRUE.   # 打开SOC
ISYM = -1          # 关闭对称性（SOC破坏某些对称）
GGA_COMPAT = .FALSE.

# 需要非共线磁性设置
LNONCOLLINEAR = .TRUE.

# 其他参数同上
...
```

**对比**：绘制有/无SOC的能带结构，检查费米面附近是否有显著差异

### 决策树

```
执行1.1磁性测试
    │
    ├─ Ce磁矩 < 0.1 μB ──→ 4f巡游 ──→ 使用标准DFT（GGA）
    │
    └─ Ce磁矩 > 0.5 μB ──→ 4f局域 ──→ 执行1.2 U扫描
                                      │
                                      ├─ U=0与U=5差异<10% ──→ 可用GGA近似
                                      │
                                      └─ U=0与U=5差异>20% ──→ 必须用DFT+U
                                                             │
                                                             └─ 选择最佳U值（如U=5）
```

**执行1.3 SOC测试**（独立于上述决策）：
- 若SOC改变能带<0.1 eV → 可忽略
- 若SOC改变能带>0.5 eV → 后续计算都要包含SOC

### 预期结果

**情况A：4f巡游**（乐观情况）
- Ce磁矩消失，U扫描影响小
- → 可用标准DFT，计算快速

**情况B：4f局域但弱关联**（中等情况）
- Ce有小磁矩（~1-2 μB），U=3-5可稳定结果
- → 使用DFT+U（U=3或5），可接受

**情况C：4f强关联**（困难情况）
- Ce磁矩大（>3 μB），U扫描差异巨大
- → 需要DMFT或合作（超出本计划范围）

---

## 阶段二：前线轨道与能带反转分析（3-4周）

### 目标

识别LaH10和CeH10在费米能附近的分子轨道性质，检查是否存在σ/π能带反转。

### 2.1 能带结构与PDOS投影

#### VASP计算流程

**步骤1：自洽计算**（`04_band_analysis/01_scf/INCAR`）
```fortran
# SCF计算
SYSTEM = LaH10_scf
ISTART = 0
ICHARG = 2
PREC = Accurate
ENCUT = 600
ISMEAR = 0
SIGMA = 0.05
EDIFF = 1E-7
LORBIT = 11       # 输出投影波函数

# K点网格（KPOINTS文件）
# Gamma-centered 12x12x12
```

**步骤2：能带计算**（`04_band_analysis/02_bands/INCAR`）
```fortran
# Non-SCF bands
SYSTEM = LaH10_bands
ISTART = 1        # 读取WAVECAR
ICHARG = 11       # 固定电荷密度
LORBIT = 11
NEDOS = 3000      # 高分辨DOS

# KPOINTS：高对称路径
# Gamma-X-M-Gamma-R-X-M-R（FCC布里渊区）
```

**KPOINTS文件**（高对称路径）：
```
k-points for band structure (FCC)
40               # 每段的点数
Line-mode
Reciprocal
  0.0 0.0 0.0    ! Gamma
  0.5 0.0 0.5    ! X

  0.5 0.0 0.5    ! X
  0.5 0.5 0.5    ! M

  0.5 0.5 0.5    ! M
  0.0 0.0 0.0    ! Gamma

  0.0 0.0 0.0    ! Gamma
  0.5 0.5 0.0    ! K
```

#### 分析脚本：能带与PDOS

**脚本**：`04_band_analysis/plot_bands_pdos.py`
```python
#!/usr/bin/env python3
"""
绘制能带结构+投影态密度
"""
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter, BSPlotter

def plot_bands_and_pdos(vasprun_file, material_name):
    """绘制能带和PDOS"""
    # 读取VASP输出
    vasprun = Vasprun(vasprun_file, parse_projected_eigen=True)
    bs = vasprun.get_band_structure(line_mode=True)
    dos = vasprun.complete_dos

    # 费米能级
    efermi = vasprun.efermi

    # 创建图形
    fig = plt.figure(figsize=(16, 6))
    gs = fig.add_gridspec(1, 2, width_ratios=[2, 1], wspace=0.05)
    ax_bands = fig.add_subplot(gs[0])
    ax_dos = fig.add_subplot(gs[1])

    # === 左图：能带结构 ===
    # 绘制所有能带
    for spin in bs.bands:
        for band_idx in range(bs.nb_bands):
            energies = bs.bands[spin][band_idx] - efermi
            distances = [d for d in bs.distance]
            ax_bands.plot(distances, energies, 'k-', linewidth=0.5, alpha=0.7)

    # 标记费米能级
    ax_bands.axhline(0, color='r', linestyle='--', linewidth=1.5, label='$E_F$')

    # 高对称点竖线
    for label, dist in zip(bs.labels_dict.keys(), bs.labels_dict.values()):
        ax_bands.axvline(dist, color='gray', linestyle='-', linewidth=0.5)

    # 设置x轴标签
    ax_bands.set_xticks([bs.labels_dict[k] for k in bs.labels_dict])
    ax_bands.set_xticklabels([k for k in bs.labels_dict])

    ax_bands.set_ylabel('Energy - $E_F$ (eV)', fontsize=14)
    ax_bands.set_ylim(-5, 5)
    ax_bands.set_title(f'{material_name} Band Structure', fontsize=16)
    ax_bands.grid(alpha=0.3)
    ax_bands.legend()

    # === 右图：PDOS ===
    # 获取元素投影
    if material_name == 'LaH10':
        metal_symbol = 'La'
    else:
        metal_symbol = 'Ce'

    # 提取PDOS
    energies_dos = dos.energies - efermi
    pdos_metal_d = dos.get_element_spd_dos(vasprun.final_structure[metal_symbol]).get_densities()[2]  # d轨道
    pdos_metal_f = dos.get_element_spd_dos(vasprun.final_structure[metal_symbol]).get_densities()[3]  # f轨道
    pdos_H_s = dos.get_element_spd_dos(vasprun.final_structure['H']).get_densities()[0]  # s轨道

    # 绘制PDOS
    ax_dos.plot(pdos_metal_d, energies_dos, label=f'{metal_symbol}-d', linewidth=2)
    ax_dos.plot(pdos_metal_f, energies_dos, label=f'{metal_symbol}-f', linewidth=2)
    ax_dos.plot(pdos_H_s, energies_dos, label='H-s', linewidth=2)

    # 费米能级
    ax_dos.axhline(0, color='r', linestyle='--', linewidth=1.5)

    ax_dos.set_xlabel('DOS (states/eV)', fontsize=14)
    ax_dos.set_ylim(-5, 5)
    ax_dos.set_yticklabels([])
    ax_dos.set_title('Projected DOS', fontsize=16)
    ax_dos.legend(loc='upper right', fontsize=10)
    ax_dos.grid(alpha=0.3)

    plt.savefig(f'{material_name}_bands_pdos.png', dpi=300, bbox_inches='tight')
    print(f"图像已保存：{material_name}_bands_pdos.png")

    # === 统计费米面交叉次数 ===
    count_crossings = 0
    for spin in bs.bands:
        for band_idx in range(bs.nb_bands):
            energies = bs.bands[spin][band_idx] - efermi
            # 检查能带是否穿过费米能（能量从负到正或从正到负）
            crossings = np.where(np.diff(np.sign(energies)))[0]
            count_crossings += len(crossings)

    print(f"\n{material_name}费米面交叉次数：{count_crossings}")

    # === 计算费米能处的态密度 ===
    idx_fermi = np.argmin(np.abs(energies_dos))
    N_EF_total = dos.densities[0][idx_fermi]  # 总DOS
    N_EF_metal_d = pdos_metal_d[idx_fermi]
    N_EF_metal_f = pdos_metal_f[idx_fermi]
    N_EF_H = pdos_H_s[idx_fermi]

    print(f"\n费米能处的态密度：")
    print(f"  总DOS: {N_EF_total:.3f} states/eV")
    print(f"  {metal_symbol}-d: {N_EF_metal_d:.3f} states/eV ({100*N_EF_metal_d/N_EF_total:.1f}%)")
    print(f"  {metal_symbol}-f: {N_EF_metal_f:.3f} states/eV ({100*N_EF_metal_f/N_EF_total:.1f}%)")
    print(f"  H-s: {N_EF_H:.3f} states/eV ({100*N_EF_H/N_EF_total:.1f}%)")

# 主程序
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("用法：python plot_bands_pdos.py <LaH10|CeH10>")
        sys.exit(1)

    material = sys.argv[1]
    vasprun_file = f"vasprun.xml"  # 假设在当前目录

    plot_bands_and_pdos(vasprun_file, material)
```

**运行**：
```bash
cd 04_band_analysis/LaH10/
python ../plot_bands_pdos.py LaH10

cd ../CeH10/
python ../plot_bands_pdos.py CeH10
```

### 2.2 Kohn-Sham轨道可视化

#### 目标

提取费米能附近的Bloch态，识别H笼子的σ/π键特征。

#### VASP设置

在能带计算中加入：
```fortran
LWAVE = .TRUE.     # 输出WAVECAR
LPARD = .TRUE.     # 输出部分电荷密度
```

#### 提取特定能带的波函数

**脚本**：`05_orbital_visualization/extract_band_charge.sh`
```bash
#!/bin/bash
# 提取特定能带的电荷密度用于可视化

# 参数
BAND_MIN=10      # 费米能附近能带的起始序号
BAND_MAX=15      # 结束序号
KPOINT=1         # Γ点（通常是1）

# 创建INCAR用于部分电荷密度计算
cat > INCAR_parchg << EOF
# Partial charge density
SYSTEM = Extract bands
ISTART = 1
ICHARG = 11
LPARD = .TRUE.
LSEPB = .TRUE.    # 分离每条能带
LSEPK = .TRUE.    # 分离每个k点
NBMOD = -3        # 模式：特定能带范围
EINT = -10 10     # 能量窗口（相对费米能，eV）
EOF

# 运行VASP
echo "提取能带 ${BAND_MIN}-${BAND_MAX} 的电荷密度..."
mpirun -np 16 vasp_std

# 重命名输出文件
for ((i=BAND_MIN; i<=BAND_MAX; i++)); do
    mv PARCHG.${i}.${KPOINT} PARCHG_band${i}_kpt${KPOINT}.vasp
    echo "已保存：PARCHG_band${i}_kpt${KPOINT}.vasp"
done
```

#### 可视化（VESTA）

**步骤**：
1. 打开VESTA软件
2. File → Open → 选择`PARCHG_band10_kpt1.vasp`
3. Edit → Isosurfaces → 添加等值面（推荐值：0.01）
4. 观察：
   - **σ键特征**：电荷密度沿H-H连线方向（pz对齐）
   - **π键特征**：电荷密度在H-H连线侧方（pxy侧向重叠）
   - **金属d/f杂化**：金属原子周围的电荷分布

**对比LaH10和CeH10的轨道差异**：
- 相同能带序号的轨道形状是否不同？
- Ce的4f轨道是否参与H笼子的成键？

### 2.3 COHP成键分析（LOBSTER）

#### 目标

定量分析La-H和Ce-H的成键/反键性质，分离σ/π贡献。

#### LOBSTER输入文件

**步骤1：用VASP生成WAVECAR和PROCAR**
```fortran
# INCAR for LOBSTER
LORBIT = 11
ISYM = -1          # 关闭对称性（LOBSTER要求）
LWAVE = .TRUE.
```

**步骤2：LOBSTER输入**（`lobsterin`）
```
# COHP analysis
COHPstartEnergy  -15.0
COHPendEnergy     10.0
COHPSteps         1000
gaussianSmearingWidth 0.1
useOriginalTetrahedronMethod
basisfunctions La 6s 5d 4f
basisfunctions H  1s
skipCOBI
skipCOOP
```

**步骤3：运行LOBSTER**
```bash
lobster-4.1.0
```

**步骤4：分析COHP**

**脚本**：`05_orbital_visualization/plot_cohp.py`
```python
#!/usr/bin/env python3
"""
绘制COHP（Crystal Orbital Hamilton Population）
"""
import numpy as np
import matplotlib.pyplot as plt

def read_cohp(filename):
    """读取COHPCAR.lobster文件"""
    data = np.loadtxt(filename, comments='#')
    energy = data[:, 0]  # eV
    cohp_total = data[:, 1]
    return energy, cohp_total

def plot_cohp_comparison(cohp_LaH, cohp_CeH, efermi=0.0):
    """对比La-H和Ce-H的COHP"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # La-H COHP
    E_LaH, COHP_LaH = read_cohp(cohp_LaH)
    axes[0].plot(COHP_LaH, E_LaH - efermi, 'b-', linewidth=2, label='La-H')
    axes[0].axhline(0, color='r', linestyle='--', linewidth=1.5, label='$E_F$')
    axes[0].axvline(0, color='gray', linestyle='-', linewidth=0.5)
    axes[0].set_xlabel('-COHP (eV)', fontsize=14)
    axes[0].set_ylabel('Energy - $E_F$ (eV)', fontsize=14)
    axes[0].set_title('La-H Bonding Analysis', fontsize=16)
    axes[0].set_ylim(-5, 5)
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    # 标注成键/反键区域
    axes[0].fill_betweenx(E_LaH[E_LaH < efermi] - efermi,
                          0, COHP_LaH[E_LaH < efermi],
                          where=(COHP_LaH[E_LaH < efermi] > 0),
                          color='green', alpha=0.2, label='Bonding')
    axes[0].fill_betweenx(E_LaH[E_LaH < efermi] - efermi,
                          0, COHP_LaH[E_LaH < efermi],
                          where=(COHP_LaH[E_LaH < efermi] < 0),
                          color='red', alpha=0.2, label='Antibonding')

    # Ce-H COHP（类似）
    E_CeH, COHP_CeH = read_cohp(cohp_CeH)
    axes[1].plot(COHP_CeH, E_CeH - efermi, 'r-', linewidth=2, label='Ce-H')
    axes[1].axhline(0, color='r', linestyle='--', linewidth=1.5)
    axes[1].axvline(0, color='gray', linestyle='-', linewidth=0.5)
    axes[1].set_xlabel('-COHP (eV)', fontsize=14)
    axes[1].set_ylabel('Energy - $E_F$ (eV)', fontsize=14)
    axes[1].set_title('Ce-H Bonding Analysis', fontsize=16)
    axes[1].set_ylim(-5, 5)
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('COHP_comparison.png', dpi=300)
    print("COHP对比图已保存：COHP_comparison.png")

    # 计算积分COHP（ICOHP）
    # 正值：成键，负值：反键
    ICOHP_LaH = np.trapz(COHP_LaH[E_LaH < efermi], E_LaH[E_LaH < efermi])
    ICOHP_CeH = np.trapz(COHP_CeH[E_CeH < efermi], E_CeH[E_CeH < efermi])

    print(f"\n积分COHP（费米能以下）：")
    print(f"  La-H: {ICOHP_LaH:.3f} eV")
    print(f"  Ce-H: {ICOHP_CeH:.3f} eV")

    if ICOHP_LaH > ICOHP_CeH:
        print(f"\n结论：La-H成键更强（差异：{ICOHP_LaH - ICOHP_CeH:.3f} eV）")
    else:
        print(f"\n结论：Ce-H成键更强（差异：{ICOHP_CeH - ICOHP_LaH:.3f} eV）")

# 主程序
if __name__ == "__main__":
    plot_cohp_comparison(
        cohp_LaH="LaH10/COHPCAR.lobster",
        cohp_CeH="CeH10/COHPCAR.lobster",
        efermi=0.0
    )
```

### 预期结果

**能带反转的证据**：
- LaH10：费米能附近主要是La-5d + H-s的σ键，穿过费米能6-8次
- CeH10：费米能附近主要是H-H π键 + Ce-4f杂化，穿过费米能3-5次
- **费米面交叉次数减少 → N(0)下降 → λ下降 → Tc下降**

**COHP分析**：
- 如果Ce-H的成键强度>La-H → Ce的4f参与成键
- 如果COHP谱形不同 → 成键类型改变（σ变π）

---

## 阶段三：关键声子与EPC机制（4-5周）

### 目标

通过冻结声子计算，确定哪些声子模式与电子态强耦合，检查能级交叉现象。

### 3.1 声子模式分解（已有数据）

**使用现有的EPW数据**：
```bash
cd EPW_results/LaH10/
python ../../scripts/analyze_phonon_modes.py

cd ../CeH10/
python ../../scripts/analyze_phonon_modes.py
```

**输出**：识别前5-10个贡献最大的声子模式及其频率

### 3.2 声子模式可视化

**工具**：Phonopy + phononwebsite

**步骤1：生成声子模式动画**
```bash
# 假设有phonopy计算结果（mesh.yaml或band.yaml）
phonopy-load --dim="2 2 2" --pa="0 0.5 0.5  0.5 0 0.5  0.5 0.5 0"
phonopy --anime=1  # 生成anime.ascii

# 上传到 http://henriquemiranda.github.io/phononwebsite/
```

**识别关键模式**：
- **摆动模式**（rocking）：H笼子的整体倾斜/扭转
- **呼吸模式**（breathing）：H笼子的径向膨胀/收缩
- **拉伸模式**（stretching）：H-H键长变化

### 3.3 冻结声子计算（核心！）

#### 原理

**冻结声子方法**：
1. 选择一个声子模式ν，频率ω_ν
2. 沿着该模式的本征矢位移原子：
   ```
   R_i = R_i^0 + A * e_iν
   ```
   其中A是振幅（单位：Å），e_iν是归一化的本征矢
3. 对不同振幅（0°, 0.5°, 1°, 2°, 3°）重算能带
4. **关键观察**：简并的能带是否分裂？分裂点的振幅是否在ZPE范围内？

#### 零点振幅（ZPE）估算

**公式**：
```
A_ZPE = sqrt(ℏ / (2 * M_eff * ω_ν))
```

其中M_eff是有效质量（从本征矢计算）。

**典型值**（H原子，ω=1000 cm⁻¹）：
```
A_ZPE ≈ 0.1 Å ≈ 2-3°的角位移
```

#### 自动化工作流

**脚本**：`06_frozen_phonon/generate_displaced_structures.py`
```python
#!/usr/bin/env python3
"""
生成冻结声子的位移结构
"""
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar

def generate_displaced_structure(
    原始结构,
    本征矢,
    振幅列表=[0.0, 0.05, 0.1, 0.15, 0.2],  # Å
    输出目录="displaced_structures"
):
    """
    生成一系列位移结构

    参数:
        原始结构: pymatgen Structure对象
        本征矢: (N_atoms, 3)数组，归一化的声子本征矢
        振幅列表: 位移振幅（Å）
    """
    import os
    os.makedirs(输出目录, exist_ok=True)

    for amp in 振幅列表:
        # 位移坐标
        structure_disp = 原始结构.copy()
        for i, site in enumerate(structure_disp):
            displacement = amp * 本征矢[i]  # 3D矢量
            site.coords += displacement

        # 保存POSCAR
        poscar = Poscar(structure_disp)
        filename = f"{输出目录}/POSCAR_amp{amp:.3f}"
        poscar.write_file(filename)
        print(f"已生成：{filename}")

    print(f"\n总共生成 {len(振幅列表)} 个位移结构")

# 示例用法
if __name__ == "__main__":
    # 读取原始结构
    structure = Structure.from_file("POSCAR_equilibrium")

    # 读取声子本征矢（需要从phonopy或VASP读取）
    # 这里假设已经提取到numpy数组
    eigenvector = np.load("phonon_mode_25_eigenvector.npy")  # 例如模式25

    # 归一化
    eigenvector /= np.linalg.norm(eigenvector)

    # 生成位移结构
    generate_displaced_structure(
        structure,
        eigenvector,
        振幅列表=[0.0, 0.05, 0.1, 0.15, 0.2, 0.25],
        输出目录="mode25_displaced"
    )
```

#### 批量计算能带

**脚本**：`06_frozen_phonon/run_all_amplitudes.sh`
```bash
#!/bin/bash
# 批量运行冻结声子计算

MODE_NUM=25
AMPLITUDES=(0.000 0.050 0.100 0.150 0.200)

for amp in "${AMPLITUDES[@]}"; do
    DIR="mode${MODE_NUM}_amp${amp}"
    mkdir -p $DIR
    cd $DIR

    # 复制INCAR, KPOINTS, POTCAR
    cp ../INCAR_bands .
    cp ../KPOINTS_bands .
    cp ../POTCAR .

    # 复制位移的POSCAR
    cp ../displaced_structures/POSCAR_amp${amp} POSCAR

    # 提交任务
    echo "提交任务：$DIR"
    sbatch submit_vasp.sh  # 或直接运行

    cd ..
done

echo "所有任务已提交"
```

#### 分析能级交叉

**脚本**：`06_frozen_phonon/analyze_band_splitting.py`
```python
#!/usr/bin/env python3
"""
分析冻结声子计算中的能级交叉
"""
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.vasp import Vasprun

def extract_bands_at_gamma(vasprun_file, efermi):
    """提取Γ点的能带能量"""
    vasprun = Vasprun(vasprun_file)
    bs = vasprun.get_band_structure()

    # Γ点通常是第一个k点
    kpoint_gamma = 0
    energies_at_gamma = []
    for spin in bs.bands:
        for band_idx in range(bs.nb_bands):
            energy = bs.bands[spin][band_idx][kpoint_gamma] - efermi
            energies_at_gamma.append(energy)

    return np.array(sorted(energies_at_gamma))

# 主程序
amplitudes = [0.0, 0.05, 0.1, 0.15, 0.2]
mode_num = 25

# 提取所有振幅的能带能量
all_energies = {}
for amp in amplitudes:
    vasprun_file = f"mode{mode_num}_amp{amp:.3f}/vasprun.xml"
    efermi = 0.0  # 从某个计算中提取
    energies = extract_bands_at_gamma(vasprun_file, efermi)
    all_energies[amp] = energies

# 绘图：能带演化
fig, ax = plt.subplots(figsize=(10, 8))

# 只绘制费米能附近的能带（例如±2 eV）
band_window = (-2.0, 2.0)

for amp in amplitudes:
    energies = all_energies[amp]
    # 筛选能带
    energies_window = energies[(energies > band_window[0]) & (energies < band_window[1])]

    # 绘制
    x = [amp] * len(energies_window)
    ax.scatter(x, energies_window, s=20, alpha=0.7)

    # 连线（追踪同一条能带）
    if amp > amplitudes[0]:
        prev_amp = amplitudes[amplitudes.index(amp) - 1]
        prev_energies = all_energies[prev_amp]
        prev_energies_window = prev_energies[(prev_energies > band_window[0]) & (prev_energies < band_window[1])]

        # 简单匹配：最近邻
        for i, e in enumerate(energies_window):
            if i < len(prev_energies_window):
                ax.plot([prev_amp, amp], [prev_energies_window[i], e], 'b-', alpha=0.3)

# 标记费米能级
ax.axhline(0, color='r', linestyle='--', linewidth=2, label='$E_F$')

ax.set_xlabel('Phonon Amplitude (Å)', fontsize=14)
ax.set_ylabel('Energy - $E_F$ (eV)', fontsize=14)
ax.set_title(f'Band Evolution: Mode {mode_num} Frozen Phonon', fontsize=16)
ax.legend()
ax.grid(alpha=0.3)

plt.savefig(f'mode{mode_num}_band_evolution.png', dpi=300)
print(f"图像已保存：mode{mode_num}_band_evolution.png")

# === 判断能级交叉 ===
print("\n能级交叉分析：")
for i in range(len(amplitudes) - 1):
    amp1, amp2 = amplitudes[i], amplitudes[i+1]
    e1 = all_energies[amp1]
    e2 = all_energies[amp2]

    # 检查能带顺序是否反转（简化判断）
    # 实际中需要更复杂的算法追踪能带
    if len(e1) == len(e2):
        order_change = np.sum(np.argsort(e1) != np.argsort(e2))
        if order_change > 0:
            print(f"  振幅 {amp1:.3f} → {amp2:.3f} Å：检测到 {order_change} 处能级顺序变化")
            print(f"    → 可能的能级交叉！")

# ZPE检查
omega_mode = 1000  # cm^-1（从声子谱读取）
M_eff = 1.008  # amu（H原子）
hbar = 1.054571817e-34  # J·s
amu_to_kg = 1.66053906660e-27
omega_SI = omega_mode * 2 * np.pi * 3e10  # 转换为rad/s

A_ZPE = np.sqrt(hbar / (2 * M_eff * amu_to_kg * omega_SI)) * 1e10  # 转换为Å
print(f"\n零点振幅（ZPE）：{A_ZPE:.3f} Å")
print(f"如果能级交叉发生在振幅<{A_ZPE:.3f} Å → 交叉是物理可达的")
```

### 3.4 耦合机制判据

**判据1：能级交叉**
- **YB6型（强EPC）**：简并带在~0.1 Å（~2°）分裂，且ZPE可达
- **LaB6型（弱EPC）**：简并保持到>0.2 Å，或ZPE不可达

**判据2：轨道类型**
- **σ键 ⊥ 声子**：轨道重叠对位移敏感 → 能级交叉
- **π键 ∥ 声子轴**：轨道重叠稳定 → 无交叉

**预期结果**：
- **LaH10**：主导声子引起σ带交叉 → λ大
- **CeH10**：π带稳定，无交叉 → λ小

---

## 阶段四：超胞与费米面嵌套（3-4周）

### 目标

通过2×2×2超胞研究M点（布里渊区角点）的额外EPC，类似YB6的Peierls效应。

### 4.1 超胞构建

**VASP输入**：
```python
# 生成2×2×2超胞
from pymatgen.core import Structure

structure = Structure.from_file("POSCAR_primitive")
supercell = structure * [2, 2, 2]  # 扩展2×2×2
supercell.to(filename="POSCAR_222supercell")
```

**计算设置**：
- 原子数：11 → 88（11×8）
- K点：原胞12×12×12 → 超胞6×6×6（倒空间缩小）

### 4.2 M点冻结声子

**M点坐标**（原胞布里渊区）：(0.5, 0.5, 0.5)

**物理意义**：
- 在2×2×2超胞中，M点折叠到Γ点
- M点的声子在超胞中变为q=0
- 如果有软化 → 预示结构不稳定性（CDW倾向）

**计算流程**：
1. 在原胞中计算M点的声子频率（phonopy）
2. 提取最软的声子模式（频率最低）
3. 按该模式位移超胞原子
4. 计算能量和能带，检查二聚化效应

### 4.3 费米面可视化与嵌套分析

**工具**：FermiSurfer + VASP

**步骤1：生成费米面数据**
```fortran
# INCAR for Fermi surface
SYSTEM = LaH10_fermisurface
ICHARG = 11
LORBIT = 11
NEDOS = 5000
KPOINTS: 高密度网格（如24×24×24）
```

**步骤2：运行FermiSurfer**
```bash
# 生成BXSF文件
fermi_surface_vasp.py  # 自定义脚本或使用现有工具
fermisurfer LaH10.bxsf
```

**分析**：
- 费米面的形状（球形、椭球、管状？）
- 是否有平行的费米面片段（好嵌套）
- 对比LaH10和CeH10的差异

### 预期结果

**LaH10**：
- 可能有一定的准1维特性（沿某个方向）
- M点声子可能轻微软化
- 类似YB6但效应较弱（3D晶格更稳定）

**CeH10**：
- 费米面形状改变（因为4f改变能带）
- 嵌套减弱 → 额外降低EPC

---

## 阶段五：压力和应变效应（2-3周，可选）

### 目标

理解外界条件如何调控能带顺序和超导性质。

### 5.1 压力依赖

**计算**：
- 150, 175, 200, 225, 250 GPa（5个压力点）
- 每个压力：优化晶格常数 → 能带 → PDOS
- 观察能带反转是否在某个压力发生

**分析**：
- 绘制Tc vs 压力曲线
- 识别最佳压力（如果存在）

### 5.2 拉伸测试

**目的**：验证"4f轨道能量"假说

**方法**：
- 将晶格拉伸5%, 10%, 15%, 20%
- 检查LaH10和CeH10的能带是否趋于相同
- **理论极限**：材料会在~10%应变时破坏（仅理论探索）

---

## 工具脚本总结

所有脚本位于`scripts/LaH10_CeH10_analysis/`目录：

| 脚本名称 | 功能 | 阶段 |
|----------|------|------|
| `analyze_magnetism.py` | 磁性测试分析 | 一 |
| `compare_U_values.py` | Hubbard U扫描对比 | 一 |
| `plot_bands_pdos.py` | 能带+PDOS绘图 | 二 |
| `extract_band_charge.sh` | 提取轨道波函数 | 二 |
| `plot_cohp.py` | COHP成键分析 | 二 |
| `generate_displaced_structures.py` | 生成冻结声子结构 | 三 |
| `run_all_amplitudes.sh` | 批量运行冻结声子 | 三 |
| `analyze_band_splitting.py` | 能级交叉分析 | 三 |
| `build_supercell.py` | 构建2×2×2超胞 | 四 |
| `plot_fermi_surface.py` | 费米面可视化 | 四 |

---

## 时间线与里程碑

| 周数 | 任务 | 可交付物 |
|------|------|----------|
| 1 | 阶段一：磁性测试 | 决定计算方案（GGA/DFT+U） |
| 2 | 阶段一：U扫描+SOC | 确认4f性质报告 |
| 3 | 阶段二：能带+PDOS | LaH10/CeH10能带对比图 |
| 4-5 | 阶段二：轨道可视化 | Kohn-Sham轨道动画 |
| 6 | 阶段二：COHP分析 | 成键差异定量报告 |
| 7-8 | 阶段三：声子分解+可视化 | 主导声子模式列表 |
| 9-11 | 阶段三：冻结声子计算 | 能级交叉证据 |
| 12 | 阶段三：EPC机制总结 | σ/π耦合机制报告 |
| 13-14 | 阶段四：超胞计算 | M点分析结果 |
| 15 | 阶段四：费米面 | 费米面对比图 |
| 16 | 阶段五：压力依赖（可选） | Tc-压力相图 |
| 17-19 | 论文撰写 | 初稿完成 |

---

## 成功标准

### 最低目标（必须达到）
✓ 明确LaH10/CeH10的前线轨道差异
✓ 识别主导声子模式
✓ 提出一个主要的物理机制解释Tc差异

### 理想目标（争取达到）
✓ 通过冻结声子观察到能级交叉差异
✓ 定量解释λ的差异（误差<30%）
✓ 发表高质量论文（PRB或以上）

### 突破性目标（如果顺利）
✓ 发现新的设计原则（预测其他稀土氢化物）
✓ 解决Ce的4f电子角色之谜
✓ 建立"化学键-超导"通用分析框架

---

## 风险管理

### 高风险
1. **Ce需要DMFT**
   - 概率：中等
   - 应对：寻求合作或使用DFT+U近似

2. **氢的非谐性**
   - 概率：高
   - 应对：SSCHA修正（计算量大，可选）

### 中风险
3. **冻结声子计算量巨大**
   - 应对：仅计算top 5模式，每个模式5个振幅

4. **CeH10实验数据缺乏**
   - 应对：以理论预测为主，建议实验验证

### 低风险
5. **能带反转未观测到**
   - 应对：拓展假说，考虑其他机制（磁性、拓扑）

---

## 参考文献

1. **YB6/LaB6研究**：Robert et al., PNAS 2024
2. **LaH10实验**：Drozdov et al., Nature 2019
3. **LaH10理论**：Errea et al., Nature 2020
4. **电声耦合理论**：Allen & Dynes, PRB 1975
5. **冻结声子方法**：Giustino et al., RMP 2017
6. **COHP方法**：Dronskowski & Blöchl, JPC 1993

---

## 附录：计算资源需求

### 硬件需求
- **CPU核数**：16-32核（并行计算）
- **内存**：64-128 GB RAM
- **存储**：~1 TB（WAVECAR和CHGCAR占用大）
- **计算时间**：单个VASP任务~10-50核时

### 软件需求
- **必须**：VASP 5.4+, Quantum ESPRESSO 6.7+
- **推荐**：Phonopy, LOBSTER, VESTA, pymatgen
- **可选**：FermiSurfer, XCrySDen, matplotlib

---

**文件位置**：`/mnt/c/Users/myth6/Desktop/sctheory/docs/04_LaH10_CeH10_research_plan.md`

**最后更新**：2025-11-04
