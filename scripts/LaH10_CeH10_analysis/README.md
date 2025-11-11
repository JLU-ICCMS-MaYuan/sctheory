# LaH10与CeH10超导比较研究：分析工具集

本目录包含用于分析LaH10和CeH10超导温度差异的所有脚本和工具。

---

## 目录结构

```
LaH10_CeH10_analysis/
├── 01_magnetic_test/          # 阶段一：磁性测试
│   └── analyze_magnetism.py
├── 02_hubbard_U/               # 阶段一：Hubbard U扫描
│   └── compare_U_values.py
├── 03_SOC_test/                # 阶段一：自旋轨道耦合
│   └── (待添加)
├── 04_band_analysis/           # 阶段二：能带和态密度
│   └── plot_bands_pdos.py
├── 05_orbital_visualization/   # 阶段二：轨道可视化
│   └── (VESTA脚本或COHP分析)
├── 06_frozen_phonon/           # 阶段三：冻结声子（核心！）
│   ├── generate_displaced_structures.py
│   └── analyze_band_splitting.py
└── README.md                   # 本文件
```

---

## 快速开始

### 前提条件

1. **软件依赖**：
   - Python 3.7+
   - NumPy, Matplotlib
   - VASP 5.4+ 或 Quantum ESPRESSO 6.7+
   - 可选：pymatgen, phonopy, VESTA, LOBSTER

2. **数据准备**：
   - LaH10和CeH10的优化结构（POSCAR）
   - VASP计算完成后的OUTCAR, EIGENVAL, DOSCAR等文件

### 工作流程

```bash
# 1. 磁性测试（确定计算方案）
cd 01_magnetic_test/
python analyze_magnetism.py

# 2. Hubbard U扫描（如果需要）
cd ../02_hubbard_U/
python compare_U_values.py

# 3. 能带和态密度分析
cd ../04_band_analysis/
python plot_bands_pdos.py LaH10 ./LaH10/EIGENVAL ./LaH10/DOSCAR
python plot_bands_pdos.py CeH10 ./CeH10/EIGENVAL ./CeH10/DOSCAR

# 4. 冻结声子分析（最重要！）
cd ../06_frozen_phonon/

# 4a. 生成位移结构
python generate_displaced_structures.py POSCAR_LaH10 \
    --mode 25 \
    --amplitudes 0.0,0.05,0.1,0.15,0.2 \
    --eigenvector mode25_eigenvector.npy \
    --output LaH10_mode25_displaced

# 4b. 提交VASP计算（每个位移结构）
# ... 运行VASP ...

# 4c. 分析能级交叉
python analyze_band_splitting.py 25 0.0,0.05,0.1,0.15,0.2
```

---

## 工具详细说明

### 阶段一：4f电子性质诊断

#### 1. `analyze_magnetism.py` - 磁性分析

**目的**：判断Ce的4f¹电子是局域还是巡游

**输入**：
- `OUTCAR_NM`：非磁性计算的OUTCAR
- `OUTCAR_FM`：铁磁计算的OUTCAR

**输出**：
- 终端输出：能量差、磁矩、建议
- `magnetic_analysis_summary.txt`：结果摘要

**运行**：
```bash
python analyze_magnetism.py
```

**判断标准**：
- Ce磁矩 < 0.1 μB → 4f巡游，用标准DFT
- Ce磁矩 > 0.5 μB → 4f局域，需要DFT+U
- 中间值 → 需要进一步测试

---

#### 2. `compare_U_values.py` - Hubbard U扫描

**目的**：测试不同U值对电子结构的影响

**前提**：
- 目录结构：`U0/`, `U3/`, `U5/`, `U7/`
- 每个目录包含OUTCAR, DOSCAR, EIGENVAL

**输出**：
- `U_scan_results.png`：能量、态密度、能隙 vs U
- 终端输出：分析和建议

**运行**：
```bash
python compare_U_values.py
```

**判断标准**：
- N(EF)变化 < 10% → 可用标准DFT
- N(EF)变化 10-30% → 使用U=3 eV
- N(EF)变化 > 30% → 使用U=5-7 eV

---

### 阶段二：能带与轨道分析

#### 3. `plot_bands_pdos.py` - 能带和态密度

**目的**：绘制能带结构和态密度，统计费米面交叉次数

**输入**：
- `EIGENVAL`：能带能量
- `DOSCAR`：态密度和费米能级

**输出**：
- `<材料名>_band_structure.png`：能带结构图
- `<材料名>_dos.png`：态密度图
- `<材料名>_electronic_summary.txt`：统计信息

**运行**：
```bash
# 单个材料
python plot_bands_pdos.py LaH10 EIGENVAL DOSCAR

# 对比两个材料
python plot_bands_pdos.py LaH10 ./LaH10/EIGENVAL ./LaH10/DOSCAR
python plot_bands_pdos.py CeH10 ./CeH10/EIGENVAL ./CeH10/DOSCAR
```

**关键输出**：
- **费米面交叉次数**：LaH10应该 > CeH10（如果能带反转）
- **N(EF)**：态密度，直接影响λ和Tc

**能带反转的证据**：
- LaH10：6-8次交叉，La-5d + H-s 主导
- CeH10：3-5次交叉，H-H π + Ce-4f 主导

---

### 阶段三：冻结声子分析（核心）

#### 4. `generate_displaced_structures.py` - 生成位移结构

**目的**：沿声子本征矢方向位移原子，生成一系列结构用于能带计算

**输入**：
- `POSCAR`：平衡结构
- `mode<N>_eigenvector.npy`：声子本征矢（numpy数组，shape=(natoms, 3)）
- 振幅列表（Å）

**输出**：
- `POSCAR_mode<N>_amp<A>`：位移后的结构文件

**运行**：
```bash
python generate_displaced_structures.py POSCAR_LaH10 \
    --mode 25 \
    --amplitudes 0.0,0.05,0.1,0.15,0.2,0.25 \
    --eigenvector mode25_eigenvector.npy \
    --output LaH10_mode25_displaced
```

**提示**：
- 本征矢可以从phonopy提取：
  ```python
  from phonopy import load
  phonon = load(phonopy_yaml='phonopy_disp.yaml')
  eigvecs = phonon.get_band_structure().eigenvectors
  np.save('mode25_eigenvector.npy', eigvecs[0, mode_index, :, :])
  ```

- 振幅范围建议：0-0.25 Å（对应0-5°的角位移）

---

#### 5. `analyze_band_splitting.py` - 能级交叉分析

**目的**：分析冻结声子计算中的能带演化，检测能级交叉

**前提**：
- 目录结构：`mode<N>_amp<A>/EIGENVAL`和`DOSCAR`
- 例如：`mode25_amp0.000/`, `mode25_amp0.050/`, ...

**输出**：
- `mode<N>_band_evolution.png`：能带演化图（两个子图）
  * 左图：所有能带（费米能±3 eV）
  * 右图：费米能附近放大（±0.5 eV）
- `mode<N>_crossing_summary.txt`：交叉分析摘要
- 终端输出：详细分析

**运行**：
```bash
# 自动检测振幅
python analyze_band_splitting.py 25

# 手动指定振幅
python analyze_band_splitting.py 25 0.0,0.05,0.1,0.15,0.2
```

**关键输出解读**：

1. **能级交叉检测**：
   ```
   振幅 0.05 → 0.10 Å：
     检测到 4 处能带顺序变化（可能的能级交叉）
       能带 10：-0.1234 eV → 0.0567 eV
       能带 11：0.0456 eV → -0.0987 eV
   ```
   → 说明简并态发生分裂和交叉

2. **简并破缺分析**：
   ```
   组 1：能带 [10, 11]，平均能量 -0.0234 eV
     最大分裂：125.3 meV @ 振幅 0.200 Å
     → 简并显著破缺！强电声耦合的标志
   ```

3. **ZPE判据**（最重要！）：
   ```
   零点振幅 A_ZPE ≈ 0.0892 Å
   第一个能级交叉发生在 ≈ 0.100 Å
     → 0.100 Å > 0.0892 Å：交叉超出ZPE
     → 结论：弱电声耦合（类似LaB6）
   ```

   **判断标准**：
   - 交叉振幅 < A_ZPE → **强EPC**（YB6型，LaH10预期）
   - 交叉振幅 > A_ZPE → **弱EPC**（LaB6型，CeH10预期）
   - 无交叉 → **极弱EPC**

---

## 完整分析流程示例

### 步骤1：Ce的4f性质诊断

```bash
# 1a. 非磁性计算
cd CeH10/magnetic_test/NM/
# ... 编辑INCAR（ISPIN=1）...
mpirun -np 16 vasp_std > vasp.log
cp OUTCAR ../OUTCAR_NM

# 1b. 铁磁计算
cd ../FM/
# ... 编辑INCAR（ISPIN=2, MAGMOM=1*5.0 32*0.0）...
mpirun -np 16 vasp_std > vasp.log
cp OUTCAR ../OUTCAR_FM

# 1c. 分析
cd ..
python ../../scripts/LaH10_CeH10_analysis/01_magnetic_test/analyze_magnetism.py

# 根据结果决定是否需要DFT+U...
```

### 步骤2：能带分析

```bash
cd LaH10/bands/
# ... 运行VASP能带计算 ...
python ../../scripts/LaH10_CeH10_analysis/04_band_analysis/plot_bands_pdos.py LaH10

cd ../../CeH10/bands/
python ../../scripts/LaH10_CeH10_analysis/04_band_analysis/plot_bands_pdos.py CeH10

# 对比费米面交叉次数：
# LaH10: 8次
# CeH10: 4次 ← 能带反转！
```

### 步骤3：声子模式识别

```bash
cd LaH10/phonopy/
# ... 运行phonopy声子计算 ...

# 从EPW数据分析主导模式
python ../../scripts/analyze_phonon_modes.py
# 输出：Mode 25 (λ=0.342, 14.8%), Mode 28 (λ=0.299, 12.9%), ...

# 提取本征矢
python << EOF
from phonopy import load
phonon = load('phonopy_disp.yaml')
bs = phonon.get_band_structure()
import numpy as np
# 假设Mode 25对应索引24（从0开始）
np.save('mode25_eigenvector.npy', bs.eigenvectors[0, 24, :, :])
EOF
```

### 步骤4：冻结声子计算（最关键！）

```bash
# 4a. 生成位移结构
cd frozen_phonon/
python ../../scripts/LaH10_CeH10_analysis/06_frozen_phonon/generate_displaced_structures.py \
    ../POSCAR_optimized \
    --mode 25 \
    --amplitudes 0.0,0.05,0.1,0.15,0.2,0.25 \
    --eigenvector ../phonopy/mode25_eigenvector.npy \
    --output mode25_displaced

# 4b. 批量运行VASP（能带计算）
for amp in 0.000 0.050 0.100 0.150 0.200 0.250; do
    dir="mode25_amp${amp}"
    mkdir -p $dir
    cd $dir
    cp ../INCAR_bands .
    cp ../KPOINTS_bands .
    cp ../POTCAR .
    cp ../mode25_displaced/POSCAR_mode25_amp${amp} POSCAR

    # 提交任务
    sbatch submit_vasp.sh  # 或 mpirun -np 16 vasp_std

    cd ..
done

# 等待所有计算完成...

# 4c. 分析能级交叉
python ../../scripts/LaH10_CeH10_analysis/06_frozen_phonon/analyze_band_splitting.py 25

# 查看结果：
#   - mode25_band_evolution.png
#   - mode25_crossing_summary.txt
```

### 步骤5：LaH10 vs CeH10对比

```bash
# 对LaH10和CeH10的相同模式重复步骤4

# 对比结果：
# LaH10（Mode 25）：
#   - 第一个交叉：0.08 Å < A_ZPE → 强EPC
#   - 简并破缺：150 meV
#
# CeH10（Mode 25）：
#   - 第一个交叉：0.20 Å > A_ZPE → 弱EPC
#   - 简并破缺：30 meV
#
# → 解释了Tc差异！
```

---

## 常见问题

### Q1: 本征矢的符号和归一化

**A**: 脚本会自动归一化本征矢。符号不重要，因为位移振幅覆盖正负（通过0）。

### Q2: 振幅的合适范围

**A**:
- 最小：0.0 Å（平衡位置）
- 最大：~0.3 Å（对应~6°角位移，接近材料破坏极限）
- 步长：0.05 Å（可以更细，但计算量大）
- **关键点**：确保覆盖ZPE振幅（~0.08-0.12 Å for H）

### Q3: 如果没有能级交叉怎么办？

**A**:
1. 检查是否选择了正确的声子模式（应该是主导EPC的模式）
2. 尝试更大的振幅（0.3-0.4 Å）
3. 如果仍无交叉 → 证明弱EPC（这本身就是重要结果！）

### Q4: 如何判断哪些声子模式最重要？

**A**:
1. 从EPW数据：`analyze_phonon_modes.py` → 找λ贡献最大的前5个
2. 可视化声子模式（phonopy website），选择：
   - **摆动模式**（rocking）：H笼子的倾斜/扭转
   - **呼吸模式**（breathing）：径向膨胀
   - 避免纯平移/转动模式（对EPC贡献小）

### Q5: VASP计算设置建议

**能带计算INCAR**：
```fortran
SYSTEM = LaH10_frozen_phonon
ISTART = 0
ICHARG = 2
PREC = Accurate
ENCUT = 600
EDIFF = 1E-7
ISMEAR = 0
SIGMA = 0.05
LORBIT = 11
NEDOS = 3000
# 如果需要PDOS：
LWAVE = .TRUE.
LPARD = .TRUE.
```

**KPOINTS（高对称路径）**：
```
k-points for band structure (FCC)
40
Line-mode
Reciprocal
  0.0 0.0 0.0   ! Gamma
  0.5 0.0 0.5   ! X
  ...
```

---

## 预期结果总结

| 材料 | 费米面交叉 | 主导轨道 | 模式25交叉振幅 | EPC强度 | λ | Tc |
|------|-----------|----------|---------------|---------|---|-----|
| LaH10 | 6-8次 | La-5d + H-s σ键 | ~0.08 Å < A_ZPE | 强 | 2.3 | 250 K |
| CeH10 | 3-5次 | H-H π + Ce-4f | ~0.20 Å > A_ZPE | 弱 | ? | ~50 K? |

**核心机制**：
- LaH10：σ键 ⊥ H振动 → 轨道重叠对位移敏感 → 能级交叉 → 强EPC
- CeH10：π键 ∥ H振动 → 轨道重叠稳定 → 无/弱交叉 → 弱EPC

---

## 进一步分析

完成以上步骤后，可以继续：

1. **COHP分析**（LOBSTER）：量化La-H vs Ce-H的成键差异
2. **轨道可视化**（VESTA）：识别σ/π键的空间分布
3. **2×2×2超胞**：研究M点的Peierls效应
4. **压力依赖**：检查能带反转是否在某个压力发生

---

## 参考文献

1. **方法论**：Robert et al., PNAS 2024（YB6/LaB6）
2. **LaH10实验**：Drozdov et al., Nature 2019
3. **冻结声子**：Giustino et al., Rev. Mod. Phys. 2017

---

## 联系与反馈

如有问题或建议，请参考主目录的`docs/04_LaH10_CeH10_research_plan.md`

**文件位置**：`/mnt/c/Users/myth6/Desktop/sctheory/scripts/LaH10_CeH10_analysis/README.md`

**最后更新**：2025-11-04
