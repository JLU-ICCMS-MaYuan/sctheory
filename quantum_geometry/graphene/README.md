# 石墨烯复现说明

本目录负责复现主文 `Fig. 2`、`SI FIG. 1-3` 和 `TABLE I`。

## 目标

第一轮只追三件事：

- 正确的 Dirac 锥和低能带结构
- `lambda_geo / lambda` 在小掺杂极限接近 50%
- `lambda_topo` 与 `lambda_geo` 趋势一致

## 子目录

- `tb/`
  - 紧束缚基线
  - quantum metric / FSM
  - `lambda_E / lambda_geo / lambda_topo` 分解
- `abinitio/jdftx/`
  - 原文同构路线
- `abinitio/qe_equiv/`
  - 等价路线：`QE + frozen phonon + Wannier90`
- `plots/`
  - 统一存放 `Fig. 2` 和 `SI FIG. 1-3` 对应图片
- `runs/`
  - 每次运行的参数、结果摘要、日志

## 原文参数

需要固定记录：

- `a = 2.46 A`
- `t = -2.751 eV`
- `gamma = -7.308 a^-2`
- `hbar*sqrt(<omega^2>) = 0.1615 eV`
- `hbar*omega_E2g(Gamma) = 0.1935 eV`
- `hbar*omega_A1'(K) = 0.1622 eV`

## 软件链

### 原文同构

- `JDFTx + frozen phonon + Wannier`

作用：

- `JDFTx`：二维石墨烯电子结构、冻声子相关信息
- `Wannier`：把 Bloch 态转到局域基，便于提取低能模型和 EPC 矩阵元

### 等价路线

- `QE + frozen phonon + Wannier90`

适用条件：

- 你希望统一第一性原理工作流到 `QE`
- 接受“物理等价，但不是原文同软件”

记录要求：

- 在 `runs/` 中明确标注“QE 等价复现”

## 不建议的路线

石墨烯不建议一开始直接用：

- `QE + ph.x + EPW`

原因：

- 不贴原文主线
- 不利于先把 EPC 写成局域轨道表象并做几何分解

## 推荐执行顺序

1. 在 `tb/` 中完成石墨烯带结构和量子几何脚本
2. 输出 `D(mu)`、`lambda_total(mu)`、`lambda_E(mu)`、`lambda_geo(mu)`、`lambda_topo(mu)`
3. 在 `abinitio/` 中准备一条你实际采用的输入链
4. 对照 `SI FIG. 2-3` 做模型与 `ab initio` 的 EPC 矩阵元比较
5. 在 `plots/` 中整理主文 `Fig. 2` 的复现图

## 当前 `tb` 数值链路

当前 `tb` 实现已经不是“直接把论文极限写死”，而是按下面的链路数值计算：

1. `build_model()` 建立带极小 `delta` 正则化的两带石墨烯模型
2. `uniform_mesh_quantities()` 在均匀 `k` 网格上输出：
   - 带结构
   - `gap`
   - quantum metric trace
   - Berry curvature
3. `estimate_dirac_winding_numbers()` 围绕 `K/K'` 小回路计算绕数
4. `lambda_decomposition_scan()` 对每个化学势做费米线积分：
   - 大掺杂区间：用全局 contour 抽取费米线
   - Dirac 点附近：改用 `extract_dirac_pocket_segments()` 围绕 `K/K'` 做局部高精度 pocket 求根
5. `integrate_fermi_surface_observables()` 在费米线上逐段积分，得到：
   - `lambda_E`
   - `lambda_geo`
   - 构造 `lambda_topo` 所需的分母积分

和论文附录 F 对应时，应这样理解：

- `lambda_total = lambda_E + lambda_geo`
- `lambda_geo` 来自费米线上的量子度规积分，不是手工把 Berry curvature 加进去
- `lambda_topo` 不是总 `lambda` 的新加项，而是由绕数给出的 `lambda_geo` 下界

## 为什么现在 `lambda_topo / lambda_geo -> 1`

这件事现在是数值结果，不是硬编码假设。原因分成三层：

### 1. 物理层

论文在石墨烯的 Dirac 极限下给出：

- `lambda_E = lambda_geo = lambda_topo = lambda / 2`

这里的意思是：

- 几何项的实际值
- 由绕数控制的拓扑下界

在低能 Dirac 模型里变成同一个值，也就是“不等式取等”。

### 2. 数值层

之前 `lambda_topo / lambda_geo` 偏离 `1`，主要不是公式错，而是费米线在 Dirac 点附近太小，用粗网格全局 contour 容易把局部口袋抽坏。现在改成：

- 围绕 `K/K'` 单独构造局部费米 pocket
- 沿每个角度径向求根，直接解 `E(k) = mu`
- 对 pocket 显式闭合后再做弧长积分

这样 `lambda_geo` 的费米线积分和 `lambda_topo` 的下界分母，都会在同一条高精度局部口袋上评估。

### 3. 正则化层

代码里保留了极小的 `onsite_delta_ev`，只作为数值稳定用的正则化，不代表论文真正想研究的有隙石墨烯。把它从 `1e-3 eV` 降到 `1e-4 eV` 后，Dirac 极限更接近无隙模型，因此：

- `lambda_geo / lambda` 更稳定地逼近 `0.5`
- `lambda_topo / lambda_geo` 更稳定地逼近 `1`

当前已验证的代表性结果：

- `mu = -0.02 eV`：`lambda_geo / lambda = 0.500023`，`lambda_topo / lambda_geo = 0.999964`
- `mu = -0.003 eV`：`lambda_geo / lambda = 0.500764`，`lambda_topo / lambda_geo = 0.999988`

因此当前代码已经具备下面这个解释力：

- `lambda_geo` 是量子度规控制的实际几何贡献
- `lambda_topo` 是由 `K/K'` Dirac 绕数给出的下界
- 在石墨烯低能极限中，这个下界数值上确实逼近实际值

## 每次回传内容

- 运行命令
- 输入参数文件路径
- 关键数值
- 生成图片路径
- 报错或异常现象
