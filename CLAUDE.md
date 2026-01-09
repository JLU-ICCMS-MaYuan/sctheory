# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目概述

这是一个高压氢化物超导理论研究项目，专注于通过第一性原理计算（DFT）、Eliashberg强耦合理论和群论分析，研究LaH10、CeH10、H3S等高压氢化物超导体的电子结构、电声耦合机制及化学键特性。

**核心科学问题**：为什么结构几乎完全相同的LaH10和CeH10，其超导转变温度（Tc）存在巨大差异（LaH10: ~250K vs CeH10: 50-100K）？

**主要假说**：
- **轨道反转**：Ce的4f轨道下移改变H笼子的成键轨道能级顺序
- **EPC机制差异**：σ键对声子敏感（强EPC），π键对声子钝化（弱EPC）
- **磁性竞争**：Ce的4f局域磁矩可能抑制超导配对

## 代码库架构

### 主要组件

1. **EliashbergEquation/** - Eliashberg方程数值求解器
   - `eliashberg_solver.py`: 从Fortran完整移植的Python实现
   - 求解各向同性Eliashberg方程，计算超导临界温度和能隙
   - 使用Padé解析延拓将Matsubara轴解延拓到实轴

2. **scripts/LaH10_CeH10_analysis/** - LaH10/CeH10对比研究工具集
   - `01_magnetic_test/`: Ce的4f电子磁性诊断
   - `02_hubbard_U/`: Hubbard U参数扫描
   - `04_band_analysis/`: 能带结构和态密度分析
   - `06_frozen_phonon/`: 冻结声子分析（核心方法）
   - 详细使用说明见该目录下的README.md

3. **MH10/, MH6/, MH9/** - 各体系的VASP/QE计算数据
   - 包含结构文件、计算脚本、Mathematica分析笔记本
   - MH10包含LaH10和CeH10的完整电声耦合数据

4. **pythtb/** - PythTB紧束缚模型库（子模块）
   - 用于构建和分析紧束缚模型
   - 核心模块：`tbmodel.py`, `lattice.py`, `wannier.py`

5. **DeePTB/** - 深度学习紧束缚参数化工具（子模块）
   - 基于E3等变神经网络的紧束缚模型训练
   - 支持从第一性原理数据学习哈密顿量

6. **mathematica_GroupTheory1.4/** - 基于GTPack的群论对称性分析
   - 用于分析晶体对称性和轨道投影
   - 包含三个主要包：GroupTheory-V1_4, SimPack_1.0, Book_Tasks
   - 典型笔记本：`MH10/*轨道投影*.nb`, `MH10/*拓扑模型*.nb`

7. **reference/** - 关键参考文献
   - 包含YB6/LaB6化学键分析论文（方法论来源）
   - 其他高压超导理论文献

8. **分子轨道mathematica代码/** - 分子轨道分析工具
   - 用于小分子和笼状结构的轨道计算

**关键架构理解**：
- **数据流**：QE计算 → ALPHA2F.OUT → Eliashberg求解器 → Tc预测
- **分析流**：VASP/QE → 冻结声子 → 能级交叉分析 → EPC机制理解
- **子模块依赖**：pythtb和DeePTB是独立Git仓库，通过子模块管理

## 常用命令

### Eliashberg方程求解

从QE的EPW计算数据求解超导临界温度：

```bash
cd EliashbergEquation/

# 基本用法
python eliashberg_solver.py --input INPUT --alpha2f ALPHA2F.OUT

# 查看帮助
python eliashberg_solver.py --help
```

**输入文件格式**：
- `INPUT`: 两行文本
  - 第1行: `mu*`（库仑赝势，通常0.10-0.13）
  - 第2行: `ntemp`（温度步数，通常100）
- `ALPHA2F.OUT`: 两列数据（频率ω, α²F(ω)），从QE的EPW计算得到

**主要输出**：
- `ELIASHBERG.OUT`: McMillan参数（λ, ω_log, Tc预测）
- `ELIASHBERG_GAP_T.OUT`: 能隙随温度变化
- `ELIASHBERG_GAP_RA.OUT`: 实轴能隙函数（用于Padé延拓验证）

**典型计算时间**：10-30秒（取决于频率点数和温度步数）

### Quantum ESPRESSO计算工作流

LaH10/CeH10的典型QE计算流程（见`MH10/LaH10-Tc-666-paw/`）：

```bash
# 步骤1-3: 结构优化和自洽计算
sbatch s123_prepare.sh
# 依次运行: relax.in -> scffit.in -> scf.in
# 输出: relax.out, scffit.out, scf.out

# 步骤4: 声子计算（无分裂，phonon without splitting）
sbatch s4_PhNoSplit.sh
# 计算动力学矩阵，输出: *.dyn 文件

# 步骤6-8: 声子后处理
sbatch s678_processphono.sh
# 生成: ALPHA2F.OUT, lambda.out, a2F.dos

# 步骤7: 声子能带数据提取（可选）
cd phonoband/
sbatch s7_phonobanddata.sh
```

**关键SLURM参数**（需根据HPC环境调整）：
- `--nodes=1 --ntasks=48`: 单节点48核
- `--partition=public`: 队列名称
- QE路径: `/work/home/mayuan/software/qe-7.4.1/bin/`（需修改为你的路径）
- MPI并行: `mpirun -np 48 pw.x -npool 4`

**计算时间估计**：
- relax: 1-3小时
- scf: 10-30分钟
- phonon: 6-24小时（取决于q点网格）

### 冻结声子分析工作流

这是研究LaH10/CeH10差异的核心方法：

```bash
cd scripts/LaH10_CeH10_analysis/06_frozen_phonon/

# 1. 生成位移结构
python generate_displaced_structures.py POSCAR \
    --mode 25 \
    --amplitudes 0.0,0.05,0.1,0.15,0.2 \
    --eigenvector mode25_eigenvector.npy \
    --output displaced_structures/

# 2. 对每个位移结构运行VASP能带计算
bash batch_submit_vasp.sh

# 3. 分析能级交叉
python analyze_band_splitting.py 25
```

关键输出：
- `mode<N>_band_evolution.png`: 能带随位移演化图
- `mode<N>_crossing_summary.txt`: 能级交叉分析摘要
- 判据：交叉振幅 < 零点振幅(ZPE) → 强电声耦合

### 紧束缚模型构建

以MH6为例：

```bash
cd MH6/

# 生成紧束缚模型
python generate_tb_model.py

# 使用PythTB验证
python verify_with_pythtb.py
```

## 关键物理概念

### Eliashberg谱函数α²F(ω)

**定义**：
```
α²F(ω) = (1/N(0)) Σ_{k,q} |g_{k,k+q}|² δ(ε_k) δ(ε_{k+q}) δ(ω - ω_q)
```

**关键关系**：
- 电声耦合常数：`λ = 2∫[α²F(ω)/ω]dω`
- **重要**：分子中的两个δ函数使得 `α²F ∝ N(0) × ⟨|g|²⟩`
- 因此：**N(0) ↑ → λ ↑ → Tc ↑**（这是常见误解的纠正！）
- N(0)在分母只是归一化因子，不改变λ的N(0)依赖性

**常见误解警告**：
很多人误以为N(0)在分母会导致N(0)越大λ越小，这是错误的！正确理解是N(0)通过两个费米面δ函数出现在分子中，因此费米面态密度越高，电声耦合越强。

### 冻结声子方法

通过沿声子本征矢方向位移原子，观察能带演化：

**物理图像**：
- 声子振动 → 原子位移 → 电子能级改变 → 能级交叉（避免交叉）
- 交叉越容易发生 → 电声耦合越强

**判据**（与零点振幅ZPE比较）：
1. **交叉振幅 < ZPE** → **强电声耦合**（LaH10预期）
   - 能级交叉发生在量子零点振动范围内
   - 系统自发处于强耦合区域

2. **交叉振幅 > ZPE** → **弱电声耦合**（CeH10预期）
   - 需要超出零点振动的大振幅才能交叉
   - 热激发难以达到强耦合区域

3. **无能级交叉** → **极弱电声耦合**
   - 能带结构对声子不敏感

**氢原子的ZPE振幅**：约0.08-0.12 Å（取决于模式频率）

这是理解LaH10（强EPC）vs CeH10（弱EPC）差异的关键物理机制。

### 轨道反转假说

LaH10 vs CeH10的核心假说：
- LaH10: La-5d + H-s σ键主导费米面 → 对声子敏感 → 强EPC
- CeH10: Ce-4f下移导致H-H π键上移 → 对声子钝化 → 弱EPC

## 数据文件格式

### VASP结构文件
- `POSCAR`: 标准VASP结构格式
- `std_*.vasp`: 标准晶胞
- `prim_*.vasp`: 原胞

### QE输入文件
- `relax.in`: 结构优化
- `scf.in`: 自洽计算
- `ph.in`: 声子计算

### 声子数据
- `ALPHA2F.OUT`: α²F(ω)谱函数
- `a2F.dos`: 分模式的α²F贡献
- `lambda.out`: 电声耦合常数

### Mathematica笔记本
- `*轨道投影*.nb`: 轨道投影分析
- `*拓扑模型*.nb`: 紧束缚模型构建
- 使用GTPack进行群论分析

## 子模块管理

项目包含两个Git子模块：

```bash
# 初始化子模块
git submodule update --init --recursive

# 更新子模块
git submodule update --remote

# 子模块位置
# - pythtb/: PythTB紧束缚库
# - DeePTB/: 深度学习紧束缚工具
```

## 研究工作流

完整的LaH10/CeH10对比研究流程：

1. **阶段一：4f电子性质诊断**
   - 磁性测试（NM vs FM）
   - Hubbard U扫描
   - SOC效应测试

2. **阶段二：前线轨道分析**
   - 能带结构和PDOS
   - 轨道可视化（VESTA）
   - COHP成键分析（LOBSTER）

3. **阶段三：关键声子与EPC机制**
   - 识别主导声子模式
   - 冻结声子计算
   - 能级交叉分析

4. **阶段四：超胞与费米面嵌套**
   - 2×2×2超胞计算
   - M点Peierls效应

详细计划见 `docs/04_LaH10_CeH10_research_plan.md`

## 重要文档

- `README.md`: 项目总览和科学问题
- `request.md`: Claude Code问答历史记录
- `docs/01_fortran_to_python_guide.md`: Eliashberg求解器移植指南
- `docs/02_questions_and_answers.md`: Eliashberg理论Q&A
- `docs/03_theory_comprehensive.md`: 电声耦合与超导理论深度解析
- `docs/04_LaH10_CeH10_research_plan.md`: LaH10/CeH10对比研究路线图
- `scripts/LaH10_CeH10_analysis/README.md`: 分析工具详细使用说明

## 计算环境

### 软件依赖
- Python 3.7+ (numpy, scipy, matplotlib)
- VASP 5.4+ 或 Quantum ESPRESSO 6.7+
- Mathematica + GTPack（群论分析）
- 可选：pymatgen, phonopy, VESTA, LOBSTER

### HPC环境
脚本中使用SLURM作业调度系统：
- 典型配置：48核节点，Intel MPI
- QE路径：默认为 `/work/home/mayuan/software/qe-7.4.1/bin/`
  - 如在不同环境运行，需修改脚本中的QE路径
  - 或设置环境变量：`export QE_BIN=/path/to/qe/bin`

### 验证安装

```bash
# 检查Python依赖
python -c "import numpy, scipy, matplotlib; print('Python环境OK')"

# 检查Eliashberg求解器
cd EliashbergEquation/
python eliashberg_solver.py --help

# 检查子模块
git submodule status
```

## 代码风格注意事项

1. **Python脚本**：遵循PEP 8，使用类型注解
2. **物理单位**：
   - 能量：eV（VASP）或Hartree（QE）
   - 温度：Kelvin
   - 长度：Angstrom
3. **文件命名**：
   - 分析脚本：`analyze_*.py`
   - 生成脚本：`generate_*.py`
   - 对比脚本：`compare_*.py`
4. **注释语言**：中文（文档和注释）+ 英文（代码和变量名）

## 常见陷阱与故障排查

### 物理概念陷阱

1. **N(0)的作用**：费米面态密度N(0)在分子中（通过两个δ函数），因此N(0) ↑ → λ ↑，而不是相反
2. **本征矢归一化**：冻结声子脚本会自动归一化，无需手动处理
3. **振幅范围**：氢原子的零点振幅约0.08-0.12 Å，位移振幅应覆盖此范围
4. **能级交叉判据**：必须与ZPE比较，而不是仅看是否有交叉
5. **子模块状态**：修改子模块后需要在主仓库中提交子模块指针的变化

### 常见错误与解决方案

#### Eliashberg求解器问题

```bash
# 错误：找不到输入文件
# 解决：检查文件路径和格式
ls -lh INPUT ALPHA2F.OUT
head -5 INPUT ALPHA2F.OUT

# 错误：收敛失败
# 解决：调整混合参数或增加迭代次数
# 编辑eliashberg_solver.py中的MIXING_BETA和MAX_IT
# MIXING_BETA: 默认0.5，可尝试0.3-0.7
# MAX_IT: 默认1000，可增加到2000
```

#### VASP/QE计算问题

```bash
# 错误：SLURM作业失败
# 解决：检查作业日志
tail -50 slurm-*.out
sacct -j <job_id> --format=JobID,JobName,State,ExitCode

# 错误：内存不足
# 解决：减少KPOINTS密度或增加节点数
# VASP: 调整NCORE参数（INCAR: NCORE = 4）
# QE: 调整-npool参数（mpirun ... pw.x -npool 4）

# 错误：声子计算不收敛
# 解决：
# 1. 增加EDIFF精度（VASP: EDIFF = 1E-8）
# 2. 检查结构是否优化完全（grep "EDIFFG" OUTCAR）
# 3. 增加tr2_ph阈值精度（QE: tr2_ph = 1.0d-14）

# 错误：QE路径找不到
# 解决：修改脚本中的QE_BIN路径
export QE_BIN=/path/to/your/qe/bin
# 或直接编辑.sh脚本中的绝对路径
```

#### 冻结声子分析问题

```bash
# 错误：本征矢维度不匹配
# 解决：检查本征矢shape是否为(natoms, 3)
python -c "import numpy as np; print(np.load('mode25_eigenvector.npy').shape)"
# 应该输出: (33, 3) for MH10 (1 M + 32 H atoms)

# 错误：找不到能级交叉
# 解决：
# 1. 确认选择了正确的声子模式（主导EPC的模式）
#    - 从λ贡献排序选择前5个模式
#    - 优先选择H原子参与的摆动/呼吸模式
# 2. 增大振幅范围（0.3-0.4 Å）
#    python generate_displaced_structures.py ... --amplitudes 0.0,0.1,0.2,0.3,0.4
# 3. 检查能带计算是否包含费米面附近的能级
#    grep "NBANDS" INCAR  # 确保足够多的能带

# 错误：EIGENVAL文件读取失败
# 解决：确认VASP计算完成且EIGENVAL格式正确
grep "NELECT" EIGENVAL
head -20 EIGENVAL  # 检查文件头格式
```

#### 子模块问题

```bash
# 错误：子模块为空（pythtb/ 或 DeePTB/ 目录为空）
# 解决：初始化子模块
git submodule update --init --recursive

# 错误：子模块版本冲突
# 解决：重置子模块到正确版本
cd pythtb/
git checkout master
git pull
cd ..
git add pythtb/
git commit -m "更新pythtb子模块"

# 错误：子模块URL无法访问（SSH密钥问题）
# 解决：检查SSH密钥配置
ssh -T git@github.com
# 或临时切换到HTTPS URL
git config --file .gitmodules submodule.pythtb.url https://github.com/JLU-ICCMS-MaYuan/pythtb.git
git submodule sync
git submodule update --init
```

#### Python环境问题

```bash
# 错误：ImportError: No module named numpy
# 解决：安装必要的Python包
pip install numpy scipy matplotlib
# 或使用conda
conda install numpy scipy matplotlib

# 错误：绘图中文显示乱码
# 解决：检查matplotlib字体配置
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS']  # macOS
# 或 ['SimHei'] for Windows
# 或 ['WenQuanYi Zen Hei'] for Linux

# 查看可用字体
from matplotlib import font_manager
fonts = [f.name for f in font_manager.fontManager.ttflist]
print([f for f in fonts if 'Arial' in f or 'Hei' in f])
```

## Git工作流

主分支：`master`

### 提交信息风格

使用中文，简洁描述变更内容。格式示例：

```bash
# 新功能
feat(MH6): implement and verify TB model for hydrogen cage
feat: 添加LaH10冻结声子分析脚本

# 完成分析
完成了H10晶格的整体分析
增加了纯氢晶格的分析

# 修复
fix(eliashberg): correct Padé extrapolation convergence
fix: 修正声子本征矢归一化问题

# 文档
docs: 更新冻结声子分析README
```

### 提交前检查

```bash
# 查看修改
git status
git diff

# 确保不提交敏感数据
git diff | grep -i "password\|token\|key"

# 检查子模块状态
git submodule status
```

## 数据分析与可视化

### 快速查看计算结果

```bash
# 查看Eliashberg求解结果
tail -20 EliashbergEquation/ELIASHBERG.OUT | grep -i "lambda\|Tc"

# 查看VASP能量收敛
grep "free  energy" OUTCAR | tail -5

# 查看QE声子频率
grep -A 10 "freq" phonon.out

# 查看费米能级
grep "E-fermi" OUTCAR                    # VASP
grep "the Fermi energy" scf.out           # QE

# 查看磁矩（VASP）
grep "magnetization (x)" OUTCAR | tail -5

# 检查SLURM作业状态
squeue -u $USER
tail -50 slurm-*.out
```

### 分析声子模式贡献

识别主导电声耦合的声子模式：

```bash
# 假设有a2F.dos文件（包含分模式α²F数据）
cd MH10/LaH10-Tc-666-paw/

# 查看lambda.out获得总λ
cat lambda.out

# 如果有分模式文件dos1, dos2, ...
# 使用analyze_phonon_modes.py分析各模式贡献
python ../../scripts/analyze_phonon_modes.py

# 输出示例：
# Mode 4: λ=0.5075 (35.22%)
# Mode 5: λ=0.3128 (21.71%)
# ...
```

### 绘图脚本位置

主要分析和绘图工具：

```bash
# 能带结构和态密度
scripts/LaH10_CeH10_analysis/04_band_analysis/plot_bands_pdos.py
# 输入: EIGENVAL, DOSCAR
# 输出: *_band_structure.png, *_dos.png, *_electronic_summary.txt

# 冻结声子演化
scripts/LaH10_CeH10_analysis/06_frozen_phonon/analyze_band_splitting.py
# 输入: mode<N>_amp<A>/EIGENVAL
# 输出: mode<N>_band_evolution.png, mode<N>_crossing_summary.txt

# 磁性分析
scripts/LaH10_CeH10_analysis/01_magnetic_test/analyze_magnetism.py
# 输入: OUTCAR_NM, OUTCAR_FM
# 输出: magnetic_analysis_summary.txt

# Hubbard U扫描
scripts/LaH10_CeH10_analysis/02_hubbard_U/compare_U_values.py
# 输入: U0/, U3/, U5/, U7/ 目录
# 输出: U_scan_results.png
```

**统一绘图风格**：
```python
import matplotlib.pyplot as plt

# 使用统一样式
plt.style.use('seaborn-v0_8')

# 支持中文显示
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS']  # macOS
plt.rcParams['axes.unicode_minus'] = False

# 高分辨率输出
plt.savefig('output.png', dpi=300, bbox_inches='tight')
```

### Mathematica群论分析

在`mathematica_GroupTheory1.4/`中使用GTPack：

```mathematica
(* 加载GTPack *)
<< GroupTheory`

(* 分析晶体点群 *)
pg = GetPointGroup[{Fm-3m}]

(* 轨道投影分析 - 参考MH10/LaH10轨道投影.nb *)
(* 1. 导入VASP能带数据 *)
(* 2. 按对称性分类轨道 *)
(* 3. 识别σ/π键特征 *)
(* 4. 可视化前线轨道 *)

(* 紧束缚模型构建 - 参考MH10/H10晶格-拓扑模型.nb *)
(* 1. 定义晶格结构 *)
(* 2. 根据对称性构建跳跃项 *)
(* 3. 拟合能带结构 *)
```

**重要笔记本文件**：
- `MH10/LaH10轨道投影.nb`: 分析LaH10的Bloch波函数轨道特征
- `MH10/H10晶格-拓扑模型(周期性).nb`: 构建H10笼子的紧束缚模型
- `MH10/M32_cage.nb`: 分析MH10笼状结构的对称性

## 性能优化建议

### VASP计算优化

```bash
# 1. 并行效率测试
# 测试不同核数的性能
for ncores in 24 48 96; do
    echo "Testing $ncores cores..."
    # 修改SLURM脚本中的--ntasks=$ncores
    # 比较wall time
done

# 2. KPOINTS优化
# 对于金属体系，使用Methfessel-Paxton展宽
# INCAR: ISMEAR = 1, SIGMA = 0.2

# 3. 内存优化
# 如果遇到内存问题，减少NCORE或增加节点
# INCAR: NCORE = 4  # 每个orbital在4个核上并行
```

### QE计算优化

```bash
# 声子计算并行
# 将不同q点分配到不同节点
mpirun -np 48 ph.x -nimage 6 -i ph.in > ph.out
# 6个image并行计算6个q点

# EPW计算优化
# 使用更大的电子动量网格加速插值
```

### Python脚本优化

```python
# 使用numpy向量化而非循环
# Bad
result = []
for i in range(len(data)):
    result.append(data[i] ** 2)

# Good
result = data ** 2

# 大数据处理使用生成器
# Bad
data = [process(x) for x in huge_list]

# Good
data = (process(x) for x in huge_list)
```

## 快速参考

### 项目关键文件速查

```
核心计算工具：
├── EliashbergEquation/eliashberg_solver.py    # Tc预测
├── scripts/LaH10_CeH10_analysis/              # 分析工具集
│   ├── 06_frozen_phonon/                      # 最重要！EPC机制分析
│   ├── 04_band_analysis/                      # 能带和费米面
│   └── 01_magnetic_test/                      # 4f电子诊断

计算数据：
├── MH10/LaH10-Tc-666-paw/                     # LaH10完整数据
├── MH10/CeH10-Tc-666-paw/                     # CeH10完整数据
└── MH10/*.nb                                   # 轨道投影和模型

文档：
├── docs/03_theory_comprehensive.md             # 理论深度解析
├── docs/04_LaH10_CeH10_research_plan.md       # 研究路线图
└── scripts/LaH10_CeH10_analysis/README.md     # 工具详细说明
```

### 最常用命令一览

```bash
# 1. 求解Tc
cd EliashbergEquation/ && python eliashberg_solver.py --input INPUT --alpha2f ALPHA2F.OUT

# 2. 提交QE计算
cd MH10/LaH10-Tc-666-paw/ && sbatch s123_prepare.sh

# 3. 冻结声子分析（核心！）
cd scripts/LaH10_CeH10_analysis/06_frozen_phonon/
python generate_displaced_structures.py POSCAR --mode 25 --amplitudes 0.0,0.05,0.1,0.15,0.2 --eigenvector mode25.npy
# ... 运行VASP ...
python analyze_band_splitting.py 25

# 4. 查看计算结果
grep "lambda" lambda.out                        # 电声耦合常数
grep "E-fermi" OUTCAR                          # 费米能级
tail -50 slurm-*.out                           # 作业日志

# 5. 子模块管理
git submodule update --init --recursive         # 初始化
git submodule status                            # 检查状态
```

### 关键物理量典型值

| 物理量 | LaH10 | CeH10（预期）| 单位 | 备注 |
|--------|-------|-------------|------|------|
| Tc | 250 | 50-100 | K | 实验/预测值 |
| λ | 2.3 | ? | - | 电声耦合常数 |
| N(EF) | 5.55 | ? | states/(eV·spin) | 费米面态密度 |
| ω_log | 900-1300 | ? | cm⁻¹ | 特征声子频率 |
| 费米面交叉 | 6-8 | 3-5（预期）| 次 | 能带反转证据 |
| ZPE振幅 | 0.08-0.12 | 0.08-0.12 | Å | H原子零点振幅 |

### 研究进度检查清单

使用这个清单追踪LaH10 vs CeH10研究进展：

- [ ] **阶段一**：4f电子诊断
  - [ ] 磁性测试（NM vs FM）
  - [ ] Hubbard U扫描
  - [ ] SOC效应测试
  - [ ] 确定计算方案（DFT vs DFT+U）

- [ ] **阶段二**：前线轨道分析
  - [ ] 能带结构计算
  - [ ] PDOS和轨道投影
  - [ ] 费米面交叉统计
  - [ ] COHP成键分析（LOBSTER）

- [ ] **阶段三**：关键声子与EPC（核心！）
  - [ ] 识别主导声子模式（从λ贡献）
  - [ ] 提取声子本征矢
  - [ ] 冻结声子计算（至少3个模式）
  - [ ] 能级交叉分析与ZPE比较
  - [ ] LaH10 vs CeH10对比

- [ ] **阶段四**：超胞与费米面嵌套
  - [ ] 2×2×2超胞计算
  - [ ] M点能带分析
  - [ ] Peierls效应检验

- [ ] **论文撰写**
  - [ ] 图表准备（能带、声子、冻结声子演化）
  - [ ] 机制图示（轨道反转示意图）
  - [ ] 数据汇总表格

### 遇到问题时的调试顺序

1. **计算不收敛**：
   - 检查EDIFF精度 → 增加迭代次数 → 检查结构优化

2. **Python脚本报错**：
   - 检查文件路径 → 验证数据格式 → 测试Python环境

3. **找不到能级交叉**：
   - 确认声子模式选择 → 增大振幅范围 → 检查能带数量

4. **SLURM作业失败**：
   - 查看日志文件 → 检查QE路径 → 验证内存/核数配置

5. **子模块问题**：
   - 初始化子模块 → 检查SSH密钥 → 考虑HTTPS URL

---

**最后更新**：2026-01-08
**维护者**：Ma Yuan (mayuan@jlu.edu.cn)
**Git仓库**：https://github.com/JLU-ICCMS-MaYuan/sctheory
