# 量子几何论文复现执行计划

本文档是 `Research/quantum-geometry-epc/quantum_geometry/` 目录的执行入口，目标是把以下两篇文献中关于石墨烯和 MgB2 的核心结果，系统地复现到本仓库中：

- `Research/quantum-geometry-epc/quantum_geometry/reference/2024 EPC-量子几何.pdf`
- `Research/quantum-geometry-epc/quantum_geometry/reference/2024 EPC-量子几何-SI.pdf`

这不是一份“泛泛而谈”的阅读提纲，而是一份按图号、按物理量、按软件链、按输出物组织的复现手册。默认工作方式是：

- 你负责外部依赖安装和实际计算
- 我负责定义每一步该算什么、看什么、如何判断是否正确
- 先拿下模型级核心结论，再做第一性原理对齐

---

## 1. 复现目标

第一轮复现不追求一步跑完整篇文章，优先拿下两条核心物理结论：

### 石墨烯

- 正确的 Dirac 锥和低能带结构
- K/K' 附近显著的量子几何结构
- `lambda_geo / lambda` 在小掺杂极限接近 `50%`
- `lambda_topo` 与 `lambda_geo` 趋势一致，并在 Dirac 点附近接近重合

### MgB2

- `pi` 与 `sigma` 两条通道可以独立计算
- `sigma` 通道对总 `lambda` 明显主导
- 主导声子模式是 `E2g(Gamma)`
- 总几何贡献接近文中给出的 `~90%`

---

## 2. 文献图表总览

当前可明确对应到：

- 主文 `3` 张图：`Fig. 1` 到 `Fig. 3`
- SI `10` 张图：`FIG. 1` 到 `FIG. 10`
- 关键数值表 `2` 个：
  - `TABLE I`：石墨烯参数表
  - `TABLE II`：MgB2 的 `lambda` 分解结果表

所有图表的详细任务映射见：

- `Research/quantum-geometry-epc/quantum_geometry/figure_map.md`

本文件只保留执行顺序和关键要求，不重复 `figure_map.md` 中的逐图索引。

---

## 3. 软件链与职责

### 3.1 石墨烯

原文石墨烯的 `ab initio` 路线不是 `QE + EPW`，而是：

- `JDFTx + frozen phonon + Wannier`

职责分工：

- `JDFTx`
  - 负责二维石墨烯的电子结构
  - 负责冻声子位移下的电子态和 EPC 原始信息
  - 负责 Wannier 化前的数据生成
- `Wannier`
  - 把 Bloch 态转到局域 Wannier 表象
  - 便于构造低能电子模型和局域 EPC 矩阵元

石墨烯原文没有使用 `EPW`，原因不是“不能用”，而是“不需要先用”。石墨烯的重点是：

- 两带 Dirac 低能模型
- `E2g(Gamma)` 和 `A1'(K)` 主导的 EPC
- `lambda_E`、`lambda_geo`、`lambda_topo` 的几何分解

因此对石墨烯，第一轮优先考虑：

- 模型级：紧束缚 + 量子几何 + 解析低能公式
- 第一性原理：`JDFTx + frozen phonon + Wannier`

如果不使用 `JDFTx`，可接受的等价路线是：

- `QE + frozen phonon + Wannier90`

但必须明确记录：

- “这是等价复现，不是原文同软件复现”

不建议石墨烯第一轮直接走：

- `QE + ph.x + EPW`

因为这不贴原文石墨烯的主线，也不利于先把 EPC 写成局域轨道基下的几何分解形式。

### 3.2 MgB2

MgB2 原文主路线是：

- `QE + DFPT + Wannier90 + EPW`

职责分工：

- `QE`
  - 做 DFT/DFPT
  - 输出粗网格电子、声子与势导数信息
- `Wannier90`
  - 对近费米面低能带做 Wannier 化
- `EPW`
  - 在致密 `k/q` 网格上插值电子和 EPC
  - 计算 `alpha^2F(omega)`、`lambda`、相关平均量

MgB2 原文还使用了一条直接校验路线：

- `QE/DFPT` 直接计算，不经过 Wannier

MgB2 必须使用 `EPW` 主线的原因是：

- 三维体系
- 多费米面片
- `pi` 和 `sigma` 两条通道并存
- 需要在致密 `k/q` 网格上稳定做 Fermi-surface average

---

## 4. 仓库目录组织

本次复现的所有新增内容只放在 `Research/quantum-geometry-epc/quantum_geometry/` 下，不污染其他研究目录。

### 根目录

- `README.md`：总览
- `figure_map.md`：逐图映射
- `plan.md`：执行计划

### 石墨烯

- `graphene/tb/`
  - 紧束缚带结构
  - quantum metric / Berry curvature / FSM
  - `lambda_E / lambda_geo / lambda_topo` 分解
- `graphene/abinitio/jdftx/`
  - 原文同构路线
- `graphene/abinitio/qe_equiv/`
  - `QE + frozen phonon + Wannier90` 等价路线
- `graphene/plots/`
  - 正式整理后的图
- `graphene/runs/`
  - 每次运行的输入、输出、摘要和日志

### MgB2

- `mgb2/tb_pi/`
  - `pz / pi` 通道
- `mgb2/tb_sigma/`
  - `px,py / sigma` 通道
  - OFSM
  - 几何与拓扑贡献分解
- `mgb2/abinitio/qe_epw/`
  - 原文主路线
- `mgb2/abinitio/qe_direct_check/`
  - 原文直接校验路线
- `mgb2/plots/`
  - 正式整理后的图
- `mgb2/runs/`
  - 每次运行的输入、输出、摘要和日志

---

## 5. 分阶段执行计划

## 阶段 0：建账与对图

目标：

- 把每张图和每个表都明确挂到某个计算任务上
- 避免后面“算了一堆量，但不知道对应哪张图”

当前产物：

- `Research/quantum-geometry-epc/quantum_geometry/README.md`
- `Research/quantum-geometry-epc/quantum_geometry/figure_map.md`

阶段完成标准：

- 主文 `Fig. 1-3`
- SI `FIG. 1-10`
- `TABLE I/II`

都已经有对应的物理量、软件链和输出路径。

---

## 阶段 A：石墨烯模型级复现

目标图：

- 主文 `Fig. 2a`
- 主文 `Fig. 2b`
- 主文 `Fig. 2c`
- 为 `SI FIG. 2-3` 做准备

### A1. 带结构

任务：

- 建立石墨烯最近邻紧束缚模型
- 验证高对称路径 `Gamma-K-M-Gamma` 上的 Dirac 锥

输出：

- `graphene_bands.png`
- 带结构原始数据 `npz/json`

验收：

- `K` 点附近出现正确的 Dirac 交叉
- 能带量级与 `t = -2.751 eV` 一致

### A2. 量子几何

任务：

- 在均匀 `k` 网格上计算 quantum metric、Berry curvature
- 重点看 K/K' 附近的几何增强

输出：

- `graphene_metric_trace.png`
- `graphene_berry_curvature.png`

验收：

- K/K' 附近几何量明显集中
- 与两带 Dirac 物理直觉一致

### A3. `lambda` 分解

任务：

- 输出化学势扫描下的：
  - `D(mu)`
  - `lambda_total(mu)`
  - `lambda_E(mu)`
  - `lambda_geo(mu)`
  - `lambda_topo(mu)`

第一轮策略：

- 现在已经从“解析低能极限占位版”升级到“数值费米线积分版”
- `lambda_total = lambda_E + lambda_geo`
- `lambda_geo` 由费米线上的量子度规积分得到
- `lambda_topo` 单独按绕数下界公式构造，不并入总 `lambda`
- 在 Dirac 点附近额外启用局部高精度 `K/K'` pocket 积分，避免粗网格 contour 破坏极限行为

输出：

- `graphene_lambda_scan.png`
- `graphene_lambda_fractions.png`
- `lambda_scan.json`

验收：

- 小掺杂极限下 `lambda_geo / lambda -> 0.5`
- `lambda_topo / lambda_geo -> 1`

当前数值说明：

- `lambda_topo / lambda_geo -> 1` 现在是数值收敛结果，不是代码中手工设定
- 物理原因是石墨烯低能 Dirac 极限下，拓扑下界对几何项取等
- 数值原因是 Dirac 点附近的费米口袋很小，必须局部构造 `K/K'` pocket，再用闭合弧长积分评估
- 稳定性原因是模型保留了极小 `onsite_delta_ev` 作为正则化；把它压小后，结果更接近论文的无隙极限

当前已验证的代表性点：

- `mu = -0.02 eV`：`lambda_geo / lambda = 0.500023`，`lambda_topo / lambda_geo = 0.999964`
- `mu = -0.003 eV`：`lambda_geo / lambda = 0.500764`，`lambda_topo / lambda_geo = 0.999988`

### A4. 当前状态

当前已经完成：

- `graphene/tb/core.py`
- `graphene/tb/run_workflow.py`

当前脚本能生成：

- 带结构
- quantum metric
- Berry curvature
- `lambda` 扫描结果
- Dirac 点附近的局部高精度 pocket 积分结果

后续待补：

- 更正式的主文 `Fig. 2a/2b/2c` 出图版式
- 与 `ab initio` 矩阵元的对比接口

---

## 阶段 B：石墨烯第一性原理基线

目标图：

- `SI FIG. 1`
- `SI FIG. 2`
- `SI FIG. 3`
- 同时为主文 `Fig. 2a` 提供 `ab initio` 曲线

### B1. 原文同构路线

- `JDFTx + frozen phonon + Wannier`

关键参数必须记录：

- `24x24x1` 电子网格
- `6x6x1` 超胞冻声子
- 7 个 trial Wannier center：
  - 2 个 C `pz`
  - 3 个键中心 `s-like`
  - 2 个六角中心上下 `s-like`

### B2. QE 等价路线

如果不用 `JDFTx`，使用：

- `QE + frozen phonon + Wannier90`

执行要求：

- 使用二维 slab 模型
- 明确记录真空层和二维体系处理方式
- 仍然先走 frozen phonon，不先走 `EPW`

### B3. 阶段产出

至少需要拿到：

- 化学势扫描下的 `lambda` 或相关基础量
- `F_tau_i(K, K+q)` 的模型/ab initio 对比基础数据
- 模型参数和 `ab initio` 参数之间的对接记录

---

## 阶段 C：MgB2 模型级复现

目标图：

- `SI FIG. 4`
- `SI FIG. 5`
- `SI FIG. 6`
- `SI FIG. 7`
- `SI FIG. 8a`
- `SI FIG. 8b`
- `SI FIG. 9`
- 主文 `Fig. 3`

MgB2 必须拆成两条线：

- `pi`
- `sigma`

不能混着做。

### C1. `pi` 通道

任务：

- 建立 `pz` 模型
- 计算带结构、费米面、量子几何量
- 分解：
  - `lambda_pi,E`
  - `lambda_pi,geo`
  - `lambda_pi,topo`

输出：

- `mgb2_pi_*`

验收：

- 能单独给出 `lambda_pi`
- 数值量级接近 `TABLE II` 中的 `~0.16`

### C2. `sigma` 通道

任务：

- 把 `Research/chemical-bonding/reference-notes/reproduce_mgb2_tba.py` 升级为完整 `px,py` 紧束缚
- 不能再用 A 点附近抛物线近似充当最终模型
- 计算 `Gamma-A` 附近有效模型
- 计算 OFSM
- 分解：
  - `lambda_sigma,E`
  - `lambda_sigma,geo`
  - `lambda_sigma,topo`

输出：

- `mgb2_sigma_*`

验收：

- `lambda_sigma` 明显主导
- `lambda_sigma,E` 很小
- `lambda_sigma,geo` 很大

### C3. 合并总量

最后合并：

- `lambda = lambda_pi + lambda_sigma`
- `lambda_geo = lambda_pi,geo + lambda_sigma,geo`
- `lambda_topo = lambda_pi,topo + lambda_sigma,topo`

目标接近：

- `lambda ≈ 0.78`
- `lambda_ab_initio ≈ 0.67`
- `lambda_geo ≈ 0.71`
- `lambda_topo ≈ 0.32`
- 几何贡献 `≈ 91.7%`

---

## 阶段 D：MgB2 第一性原理对齐

目标图：

- `SI FIG. 8a`
- `SI FIG. 8b`
- `SI FIG. 8c`
- 需要时补 `SI FIG. 10`

### D1. 主路线

- `QE + DFPT + Wannier90 + EPW`

关键参数：

- `12x12x12` 电子网格
- `6x6x6` 声子网格
- `60x60x60` 插值 `k/q` 网格
- 5 个 Wannier 初始中心：
  - 2 个 B `pz`
  - 3 个 `s-like` center

目标输出：

- `alpha^2F(omega)`
- 总 `lambda`
- `Gamma_nm` 的模型/ab initio 对比

### D2. 校验路线

- `QE/DFPT` 直接计算，不经过 Wannier

原文参考：

- `lambda_I = 0.67`
- `lambda_II = 0.645`
- 相对差约 `3.7%`

这一步优先级低于主路线，不作为第一轮阻塞项。

---

## 6. 分阶段回传格式

每次只回传以下五类信息：

- 运行命令
- 输入参数文件路径
- 关键数值
- 生成图片路径
- 报错或异常现象

推荐回传顺序：

1. 石墨烯模型级
2. 石墨烯第一性原理基线
3. MgB2 `pi`
4. MgB2 `sigma`
5. MgB2 总 `lambda`
6. MgB2 第一性原理主路线
7. MgB2 双路线一致性检查

---

## 7. 验收标准

### 石墨烯

- 正确 Dirac 锥
- K/K' 附近量子几何量增强
- `lambda_geo / lambda` 在小掺杂极限接近 `50%`
- `lambda_topo` 与 `lambda_geo` 趋势一致

### MgB2

- `pi`、`sigma` 两条通道都可单独运行
- `sigma` 通道主导总 `lambda`
- `E2g(Gamma)` 是主导模式
- 总几何贡献接近 `~90%`

### 工作流

- 每个阶段能通过单条命令重复运行
- 所有结果都在 `runs/` 中留下输入、输出和摘要

---

## 8. 当前进展

截至当前，已经完成：

- `Research/quantum-geometry-epc/quantum_geometry/` 目录骨架
- `README.md`
- `figure_map.md`
- 石墨烯 `tb` 脚本第一版：
  - `graphene/tb/core.py`
  - `graphene/tb/run_workflow.py`

当前未完成：

- 石墨烯 `ab initio` 输入模板与基线结果
- MgB2 `pi` 和 `sigma` 模型脚本
- MgB2 `QE + EPW` 输入模板

---

## 9. 默认决策

- 所有新增内容只放在 `Research/quantum-geometry-epc/quantum_geometry/`
- 石墨烯默认先做模型级，再接 `JDFTx` 或 `QE` 等价路线
- MgB2 默认按 `QE + DFPT + Wannier90 + EPW` 准备主路线
- `SI FIG. 10` 不是第一轮阻塞项
- 第一轮先追结论，不先追论文里每个数值点的高精度重合
