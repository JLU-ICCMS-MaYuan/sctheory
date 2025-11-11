# LaH10/CeH10研究快速启动指南

本文档提供最简洁的操作步骤，帮助你快速开始LaH10和CeH10的比较研究。

---

## 前提准备

### 1. 软件环境

```bash
# 必需
- VASP 5.4+ 或 Quantum ESPRESSO 6.7+
- Python 3.7+ (numpy, matplotlib, scipy)

# 推荐
pip install pymatgen phonopy
```

### 2. 数据准备

你需要已经完成：
- ✅ LaH10和CeH10的结构优化（POSCAR）
- ✅ 声子谱计算（phonopy或QE）
- ✅ 电声耦合计算（EPW），得到α²F(ω)和λ

如果还没有这些数据，参考：[04_LaH10_CeH10_research_plan.md](04_LaH10_CeH10_research_plan.md)

---

## 三步核心流程

### 第一步：确定Ce的4f性质（1-2天）

**目标**：判断是用标准DFT还是DFT+U

```bash
# 1. 运行非磁性和铁磁计算
cd CeH10/magnetic_test/

# 非磁性
cd NM/
# 编辑INCAR：ISPIN=1, MAGMOM=33*0.0
mpirun -np 16 vasp_std
cp OUTCAR ../OUTCAR_NM

# 铁磁
cd ../FM/
# 编辑INCAR：ISPIN=2, MAGMOM=1*5.0 32*0.0
mpirun -np 16 vasp_std
cp OUTCAR ../OUTCAR_FM

# 2. 分析结果
cd ..
python ../../scripts/LaH10_CeH10_analysis/01_magnetic_test/analyze_magnetism.py

# 3. 根据输出决定：
#    - Ce磁矩 < 0.1 μB → 用标准DFT
#    - Ce磁矩 > 0.5 μB → 做DFT+U扫描（U=0,3,5,7）
```

**预期用时**：1-2天（取决于计算资源）

---

### 第二步：能带对比找能带反转（2-3天）

**目标**：证明LaH10和CeH10的费米能附近轨道不同

```bash
# 1. 计算能带结构
cd LaH10/bands/
# 复制INCAR_bands, KPOINTS_bands, POTCAR
# POSCAR使用优化后的结构
mpirun -np 16 vasp_std

cd ../../CeH10/bands/
# 同样操作
mpirun -np 16 vasp_std

# 2. 绘制能带和PDOS
cd ../../LaH10/bands/
python ../../scripts/LaH10_CeH10_analysis/04_band_analysis/plot_bands_pdos.py LaH10

cd ../../CeH10/bands/
python ../../scripts/LaH10_CeH10_analysis/04_band_analysis/plot_bands_pdos.py CeH10

# 3. 对比两张图，检查：
#    - 费米面交叉次数（LaH10应该更多）
#    - N(EF)的值（LaH10应该更大）
#    - PDOS中La-5d vs Ce-4f的占比
```

**关键证据**：
- ✅ LaH10交叉6-8次，CeH10交叉3-5次 → 能带反转！
- ✅ LaH10的PDOS：La-5d主导
- ✅ CeH10的PDOS：Ce-4f显著贡献

**预期用时**：2-3天

---

### 第三步：冻结声子证明EPC差异（1-2周，最关键！）

**目标**：通过能级交叉直接证明LaH10是强EPC，CeH10是弱EPC

#### 3.1 识别主导声子模式

```bash
cd LaH10/phonopy/
python ../../scripts/analyze_phonon_modes.py

# 输出示例：
# Mode 25: λ=0.342 (14.8%)  ← 贡献最大
# Mode 28: λ=0.299 (12.9%)
# ...

# 选择前3个模式进行冻结声子分析
```

#### 3.2 提取声子本征矢

```python
# 在LaH10/phonopy/目录下运行
python << EOF
from phonopy import load
import numpy as np

phonon = load('phonopy_disp.yaml')
bs = phonon.get_band_structure()

# 假设Mode 25对应Γ点第24个模式（从0开始）
kpoint_index = 0  # Γ点
mode_index = 24   # Mode 25

eigvec = bs.eigenvectors[kpoint_index, mode_index, :, :]
np.save('mode25_eigenvector.npy', eigvec)
print(f"本征矢形状：{eigvec.shape}")
EOF
```

#### 3.3 生成位移结构

```bash
cd ../frozen_phonon/

python ../../scripts/LaH10_CeH10_analysis/06_frozen_phonon/generate_displaced_structures.py \
    ../POSCAR_optimized \
    --mode 25 \
    --amplitudes 0.0,0.05,0.1,0.15,0.2,0.25 \
    --eigenvector ../phonopy/mode25_eigenvector.npy \
    --output mode25_displaced

# 输出：
# displaced_structures/POSCAR_mode25_amp0.0000
# displaced_structures/POSCAR_mode25_amp0.0500
# ...
```

#### 3.4 批量运行VASP

```bash
# 方法1：使用批量提交脚本
bash ../../scripts/LaH10_CeH10_analysis/06_frozen_phonon/batch_submit_vasp.sh \
    25 "0.000 0.050 0.100 0.150 0.200 0.250"

# 方法2：手动逐个运行
for amp in 0.000 0.050 0.100 0.150 0.200 0.250; do
    dir="mode25_amp${amp}"
    mkdir -p $dir
    cd $dir
    cp ../templates/INCAR_bands INCAR
    cp ../templates/KPOINTS_bands KPOINTS
    cp ../POTCAR .
    cp ../mode25_displaced/POSCAR_mode25_amp${amp} POSCAR

    # 提交或运行VASP
    sbatch submit_vasp.sh  # 或 mpirun -np 16 vasp_std

    cd ..
done

# 等待所有计算完成（可能需要数小时到数天）
```

#### 3.5 分析能级交叉

```bash
# 所有VASP计算完成后
python ../../scripts/LaH10_CeH10_analysis/06_frozen_phonon/analyze_band_splitting.py 25

# 输出：
# - mode25_band_evolution.png （能带演化图）
# - mode25_crossing_summary.txt （分析报告）
```

**关键判据**：

查看输出中的"零点振幅判据"部分：

```
零点振幅 A_ZPE ≈ 0.0892 Å
第一个能级交叉发生在 ≈ 0.100 Å
  → 0.100 Å > 0.0892 Å：交叉超出ZPE
  → 结论：弱电声耦合（类似LaB6）
```

**预期结果**：
- ✅ **LaH10**：交叉振幅 ~0.08 Å < A_ZPE → **强EPC**
- ✅ **CeH10**：交叉振幅 ~0.20 Å > A_ZPE → **弱EPC**
- ✅ **解释了Tc差异的微观机制！**

**预期用时**：1-2周（主要是VASP计算时间）

---

## 对比CeH10

对CeH10重复第三步的3.2-3.5：

```bash
cd ../../CeH10/phonopy/
# 提取相同模式的本征矢
python << EOF
from phonopy import load
import numpy as np
phonon = load('phonopy_disp.yaml')
bs = phonon.get_band_structure()
np.save('mode25_eigenvector.npy', bs.eigenvectors[0, 24, :, :])
EOF

cd ../frozen_phonon/
# 生成位移结构
python ../../scripts/LaH10_CeH10_analysis/06_frozen_phonon/generate_displaced_structures.py \
    ../POSCAR_optimized \
    --mode 25 \
    --amplitudes 0.0,0.05,0.1,0.15,0.2,0.25 \
    --eigenvector ../phonopy/mode25_eigenvector.npy \
    --output mode25_displaced

# 批量提交VASP
bash ../../scripts/LaH10_CeH10_analysis/06_frozen_phonon/batch_submit_vasp.sh \
    25 "0.000 0.050 0.100 0.150 0.200 0.250"

# 等待完成后分析
python ../../scripts/LaH10_CeH10_analysis/06_frozen_phonon/analyze_band_splitting.py 25
```

---

## 最终对比表

| 材料 | 费米面交叉 | 主导轨道 | 模式25交叉振幅 | vs A_ZPE | EPC | λ | Tc |
|------|-----------|----------|---------------|----------|-----|---|-----|
| LaH10 | 6-8次 | La-5d σ | ~0.08 Å | < A_ZPE | **强** | 2.3 | 250 K |
| CeH10 | 3-5次 | Ce-4f π | ~0.20 Å | > A_ZPE | **弱** | ? | ~50 K |

**物理图像**：
- LaH10：σ键 ⊥ H振动 → 重叠敏感 → 能级交叉在ZPE内 → 强EPC → 高Tc
- CeH10：π键 ∥ H振动 → 重叠稳定 → 能级交叉超出ZPE → 弱EPC → 低Tc

**核心机制**：Ce的4f¹电子将π系统拉低到费米能，改变了电声耦合机制！

---

## 常见问题

### Q1: 计算量太大怎么办？

**A**: 优先级排序：
1. **必须做**：至少分析一个主导声子模式（如Mode 25）
2. **推荐做**：分析前3个主导模式
3. **完整研究**：分析前5-10个模式

每个模式6个振幅 × 2个材料 = 12个VASP计算
如果资源有限，先做LaH10的Mode 25验证方法可行性。

### Q2: 如何知道计算收敛了？

**A**: 检查每个`mode25_amp*/OUTCAR`文件：
```bash
for dir in mode25_amp*/; do
    echo -n "$dir: "
    grep "reached required accuracy" $dir/OUTCAR | tail -1
done
```

应该看到"reached required accuracy"

### Q3: 能级交叉不明显怎么办？

**A**:
1. 增加振幅密度：0.0, 0.025, 0.05, 0.075, ... (间隔更小)
2. 尝试其他声子模式（可能选错了）
3. 放大绘图的能量窗口（analyze_band_splitting.py会自动聚焦±0.5 eV）

### Q4: 两个材料都没交叉？

**A**:
- 可能选错了声子模式，换一个试试
- 或者这个模式确实不耦合（这本身也是重要发现）
- 检查VASP计算的k点网格是否够密（特别是Γ点附近）

---

## 下一步

完成以上三步后，你已经有了核心证据。接下来可以：

1. **扩展研究**（可选）：
   - 分析更多声子模式
   - 2×2×2超胞研究M点效应
   - 压力依赖性
   - COHP成键分析（LOBSTER）
   - 轨道可视化（VESTA）

2. **撰写论文**：
   - 整理所有图表
   - 对比文献（YB6/LaB6, 其他氢化物）
   - 撰写初稿

3. **预测与验证**：
   - 预测其他稀土氢化物（PrH10, NdH10...）的Tc
   - 建议实验测量CeH10的Tc和光谱

---

## 时间估算

| 步骤 | 最快 | 一般 | 完整 |
|------|------|------|------|
| 第一步（4f性质） | 1天 | 2天 | 1周 |
| 第二步（能带） | 1天 | 3天 | 1周 |
| 第三步（冻结声子） | 3天 | 1周 | 2周 |
| **总计** | **5天** | **2周** | **4周** |

"最快"模式：只做LaH10的Mode 25，验证方法
"一般"模式：LaH10和CeH10各做1-2个模式
"完整"模式：按照完整研究计划执行

---

## 成功标志

当你完成以上步骤并且：
- ✅ LaH10的能级交叉在ZPE内（强EPC）
- ✅ CeH10的能级交叉在ZPE外或无交叉（弱EPC）
- ✅ 两者的费米面交叉次数不同（能带反转）

**恭喜！你已经找到了解释Tc差异的微观机制！**

这足以支撑一篇高质量的理论论文（PRB或更好）。

---

**文件位置**：`docs/QUICKSTART.md`
**参考详细计划**：`docs/04_LaH10_CeH10_research_plan.md`
**工具文档**：`scripts/LaH10_CeH10_analysis/README.md`

**最后更新**：2025-11-04
