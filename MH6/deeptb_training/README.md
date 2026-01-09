# MH6 DeePTB训练项目

纯H6笼子 vs CaH6完整体系的Gamma点能级对比研究

## 项目目标

使用DeePTB-SK方法训练紧束缚模型，分析Ca原子对H6笼子能级结构的影响，验证是否存在类似LaH10/CeH10的轨道反转现象。

## 目录结构

```
deeptb_training/
├── data_generation/              # 数据生成脚本
│   ├── generate_tb_eigenvalues.py
│   ├── prepare_h6_data.py
│   └── prepare_cah6_data.py
├── h6_pure/                      # 纯H6训练
│   ├── data/set.0/              # 训练数据
│   ├── input.json               # 训练配置
│   └── output/                  # 训练结果（自动生成）
├── cah6_full/                    # CaH6训练
│   ├── data/set.0/
│   ├── input.json
│   └── output/
├── band_comparison/              # 能带对比分析
│   ├── band.json
│   ├── analyze_gamma_levels.py
│   └── plot_level_comparison.py
└── README.md                     # 本文件
```

## 快速开始

### 1. 安装依赖

```bash
# 安装DeePTB
cd /home/mayuan/code/sctheory/DeePTB
pip install -e .

# 安装额外依赖
pip install ase numpy matplotlib
```

### 2. 生成训练数据

```bash
cd deeptb_training/data_generation

# 生成纯H6数据
python prepare_h6_data.py

# 生成CaH6数据
python prepare_cah6_data.py
```

**预期输出**：
- `h6_pure/data/set.0/`: eigenvalues.npy, kpoints.npy, xdat.traj, info.json
- `cah6_full/data/set.0/`: 相同文件

### 3. 训练模型

```bash
# 训练纯H6模型（~15分钟，CPU）
cd ../h6_pure
dptb train input.json -o ./output

# 训练CaH6模型（~20分钟，CPU）
cd ../cah6_full
dptb train input.json -o ./output
```

**检查训练进度**：
```bash
# 查看loss曲线
tail -f output/log/log.txt

# 检查最佳模型
ls -lh output/checkpoint/best_nnsk.pth
```

### 4. 计算能带结构

```bash
cd ../band_comparison

# 计算纯H6能带
dptb run band.json \
    -i ../h6_pure/output/checkpoint/best_nnsk.pth \
    -o ./h6_bands

# 计算CaH6能带
dptb run band.json \
    -i ../cah6_full/output/checkpoint/best_nnsk.pth \
    -o ./cah6_bands
```

### 5. 分析Gamma点能级

```bash
# 提取并分析能级
python analyze_gamma_levels.py

# 绘制对比图
python plot_level_comparison.py
```

**输出结果**：
- `gamma_levels_comparison.txt` - 文本对比
- `gamma_levels_plot.png` - 能级棒图
- `energy_shifts.png` - 能量偏移图

## 关键参数说明

### 紧束缚模型参数

| 参数 | 纯H6 | CaH6 | 说明 |
|------|------|------|------|
| H-H跳跃 | -2.0 eV | -1.5 eV (平均) | 最近邻跳跃积分 |
| H onsite | 0.0 eV | 0.0 eV | H的onsite能级 |
| Ca onsite | - | -2.0 eV | Ca的onsite能级 |
| 截断距离 | 1.3 Å | 2.5 Å | 跳跃截断半径 |

### DeePTB训练参数

- **基组**: H用单s轨道，Ca用单s轨道（简化模型）
- **方法**: powerlaw幂律衰减
- **epoch**: 500轮
- **学习率**: 0.01，使用Reduce on Plateau调度
- **损失函数**: eigvals（拟合本征值）

## 预期结果

### 成功标志
- ✅ 训练loss < 0.01 eV
- ✅ 能带图与原始TB模型匹配（误差 < 0.1 eV）
- ✅ Gamma点能级清晰分类（成键/反键）
- ✅ 检测到能级顺序差异（如果存在）

### 科学发现预期
根据轨道反转假说：
- **纯H6**: σ/π杂化轨道在费米面附近
- **CaH6**: Ca的2个电子填充导致：
  - 费米能级上移
  - 能级间距改变
  - **可能的能级反转**（类似LaH10→CeH10）

## 故障排查

### 问题1：训练不收敛
```bash
# 降低学习率
# 编辑input.json: "lr": 0.001

# 增加训练轮数
# 编辑input.json: "num_epoch": 1000
```

### 问题2：能带不匹配
这是正常现象！DeePTB拟合的是TB模型参数，误差在0.5 eV内可接受。

### 问题3：内存不足
```bash
# 减小batch_size（虽然已经是1了）
# 或关闭不必要的程序
```

## 下一步工作

### 用真实DFT数据替换
1. 运行VASP能带计算：
   ```bash
   # 设置INCAR, KPOINTS, POTCAR
   vasp_std > vasp.log
   ```

2. 转换EIGENVAL为DeePTB格式：
   ```python
   from ase.io import read
   import numpy as np

   # 读取EIGENVAL并提取数据
   # 保存为eigenvalues.npy格式
   ```

3. 重新训练获得高精度模型

### 冻结声子分析
如果发现能级反转，可以进一步：
1. 选择关键声子模式
2. 生成位移结构
3. 分析能级交叉振幅
4. 与ZPE比较，判断EPC强度

## 参考文献

- DeePTB文档: /home/mayuan/code/sctheory/DeePTB/docs/
- 项目总体说明: /home/mayuan/code/sctheory/CLAUDE.md
- 研究计划: /home/mayuan/code/sctheory/docs/04_LaH10_CeH10_research_plan.md

## 联系

如有问题，参考主项目的request.md记录或CLAUDE.md指南。

**最后更新**: 2026-01-09
