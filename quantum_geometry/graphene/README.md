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

- `a = 2.46 Å`
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
   - 负化学势的 Dirac 锥区间：优先用 `extract_dirac_pocket_segments()` 围绕 `K/K'` 做局部高精度 pocket 求根
   - Dirac 求根失败或不适用时：回退到全局 contour 抽取费米线
5. `integrate_fermi_surface_observables()` 在费米线上逐段积分，得到：
   - `lambda_E`
   - `lambda_geo`
   - 构造 `lambda_topo` 所需的分母积分

和论文附录 F 对应时，应这样理解：

- `lambda_total = lambda_E + lambda_geo`
- `lambda_geo` 来自费米线上的量子度规积分，不是手工把 Berry curvature 加进去
- `lambda_topo` 不是总 `lambda` 的新加项，而是由绕数给出的 `lambda_geo` 下界

## `tb` 代码如何对应论文公式

`tb` 目录的计算可以理解为从“论文公式”到“数值输出”的一条流水线：先建立最近邻石墨烯两带模型，再在 `k` 空间计算能带和量子几何，最后沿费米线积分得到 `lambda_E`、`lambda_geo` 和 `lambda_topo`。

### 1. 固定论文模型参数

入口在 `tb/run_workflow.py` 的 `main()`，首先创建 `GrapheneParameters`。这些参数来自 `SI TABLE I` 和 `SI Eq. (F93)-(F96)`：

- `a = 2.46 A`
- `t = -2.751 eV`
- `gamma = -7.308 a^-2`
- `hbar*sqrt(<omega^2>) = 0.1615 eV`
- `m_C = 12 amu`

代码中的 `unit_cell_area_a2`、`gamma_ev_per_a2`、`mean_omega2_ev2` 和 `carbon_mass_inv_ev_a2` 用来组装 `SI Eq. (F72)` 和 `SI Eq. (F78)` 中的 prefactor。这里 `carbon_mass_inv_ev_a2` 表示 `m_C / hbar^2`，单位是 `1/(eV Å^2)`：

```text
Omega * gamma^2 / [(2*pi)^2 * m_C * <omega^2>]
Omega * gamma^2 / [4 * m_C * <omega^2>]
```

### 2. 构造两带石墨烯哈密顿量

`core.py::build_model()` 调用 `pythtb.models.graphene(delta, t)`，构造蜂窝晶格 A/B 两子晶格最近邻模型，对应 `SI Eq. (F14)-(F15)`。

代码保留了一个很小的 `onsite_delta_ev = 1e-4 eV`。它只是数值正则化，避免 Dirac 点严格简并时量子度规、Berry 曲率和相位绕数出现不稳定；论文的解析结论对应无隙极限。

### 3. 计算能带与量子几何

`core.py::uniform_mesh_quantities()` 在均匀 `k` 网格上计算：

- `evals_ev`：两条能带 `E_1(k), E_2(k)`
- `gap_ev`：`Delta E(k)=E_2(k)-E_1(k)`
- `metric_tensor`：Fubini-Study metric
- `metric_trace`：`Tr g = g_xx + g_yy`
- `berry_curvature`：Berry 曲率
- `velocity_norm`：`|grad_k E_1(k)|`

这里最关键的是量子度规，对应 `SI Eq. (F69)`：

$$ [g_n(k)]_ij = 1/2 Tr[(d_i P_n(k)) (d_j P_n(k))] $$

`metric_trace` 进入 `SI Eq. (F72)` 的 `lambda_geo`。Berry 曲率当前只作为几何/拓扑诊断量输出，不直接进入 `lambda` 分解。

### 4. 抽取费米线

`core.py::lambda_decomposition_scan()` 对每个化学势 `mu` 选择费米线抽取方式：

- 负化学势的 Dirac 锥区间：`extract_dirac_pocket_segments()` 围绕 `K/K'` 沿不同角度径向求根，直接解 `E(k)=mu`，构造两个小 pocket。
- Dirac 求根失败或不适用时：`extract_fermi_surface_segments()` 在全局网格上抽取 `E_1(k)=mu` 的等值线。

Dirac pocket 方式对应 `SI Eq. (F83)-(F86)` 的图像，可以避免 `-0.2 eV` 到 Dirac 点附近混用粗网格 contour 带来的非线性数值误差。

### 5. 沿费米线积分

`core.py::integrate_fermi_surface_observables()` 对每条费米线线段取中点，并调用 `evaluate_band_quantities_on_points()` 重新计算：

- `Delta E(k)`
- `Tr g(k)`
- `|grad_k E(k)|`

然后执行四个线积分：

```text
dos_integral        = integral_FS ds / |grad E|
lambda_e_integral   = integral_FS ds * |grad E|
lambda_geo_integral = integral_FS ds * DeltaE^2 * Tr(g) / |grad E|
fs_inv_gap_integral = integral_FS ds * |grad E| / DeltaE^2
```

对应关系如下：

- `dos_integral`：对应 DOS 定义 `SI Eq. (B72)` 经过 `SI Eq. (B73)` 转成费米线积分后的积分核
- `lambda_e_integral`：对应 `SI Eq. (F64)-(F65)` 和 `SI Eq. (F72)` 的色散项
- `lambda_geo_integral`：对应主文 `Eq. (7)`、`SI Eq. (F68)` 和 `SI Eq. (F72)` 的几何项
- `fs_inv_gap_integral`：对应主文 `Eq. (10)` 和 `SI Eq. (F78)` 中 `lambda_topo` 的分母

### 6. 计算 Dirac winding 与拓扑下界

`core.py::estimate_dirac_winding_numbers()` 围绕 `K` 和 `K'` 各取一个小回路，读取哈密顿量非对角元 `h_AB(k)` 的相位绕数：

```text
W = (1 / 2*pi) * integral d theta(k)
```

这对应 `SI Eq. (F18)`、`SI Eq. (F20)-(F21)`。理想石墨烯中，`K` 和 `K'` 的 winding 分别为 `+1` 和 `-1`，代码取绝对值后得到：

```text
|W_K| + |W_K'| = 2
```

随后 `lambda_decomposition_scan()` 用主文 `Eq. (10)` / `SI Eq. (F78)` 计算：

```text
lambda_topo =
    [Omega * gamma^2 / (4 * m_C * <omega^2>)]
    * (|W_K| + |W_K'|)^2
    / integral_FS ds * |grad E| / DeltaE^2
```

这里的 `lambda_topo` 是 `lambda_geo` 的拓扑下界，不是额外加入总 `lambda` 的新贡献。

### 7. 组装 `lambda` 分解

`core.py::lambda_decomposition_scan()` 最终把费米线积分乘上 prefactor，得到：

```text
lambda_E   = C * integral_FS ds * |grad E|
lambda_geo = C * integral_FS ds * DeltaE^2 * Tr(g) / |grad E|
lambda     = lambda_E + lambda_geo
```

这对应 `SI Eq. (F71)-(F72)`。石墨烯中 `lambda_E-geo` 交叉项因 `C3` 对称性和子晶格手征结构为 0，对应 `SI Eq. (F62)-(F63)`，因此代码没有单独计算该项。

在 Dirac 小掺杂极限，论文给出 `SI Eq. (F89)-(F92)`：

```text
lambda_E = lambda_geo = lambda_topo = lambda / 2
```

当前代码的 `lambda_geo / lambda -> 0.5`、`lambda_topo / lambda -> 0.5` 和 `lambda_topo / lambda_geo -> 1` 就是在数值上复现这个极限。本项目的占比图按当前约定画 `lambda_geo / lambda` 和 `lambda_topo / lambda_geo`。

### 8. 输出文件和图

`tb/run_workflow.py` 保存以下结果：

- `data/graphene_mesh_data.npz`：`k` 网格、能带、gap、`metric_trace`、Berry 曲率、速度
- `data/graphene_band_path.npz`：高对称路径能带
- `data/parameters.json`：`GrapheneParameters`
- `data/winding_numbers.json`：`K/K'` winding number
- `data/lambda_scan.json`：不同 `mu` 下的 `lambda` 分解
- `summary.json`：运行摘要和小掺杂样例点
- `plots/*.png`：能带图、量子度规图、Berry 曲率图、`lambda` 扫描图、占比图

`tb/structure.mmd` 是这条计算链的结构图，节点中同时标注了函数对象、代码位置、论文公式和物理意义。

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
- `lambda_topo / lambda` 更稳定地逼近 `0.5`
- `lambda_topo / lambda_geo` 更稳定地逼近 `1`

当前已验证的代表性结果：

- `mu = -0.02 eV`：`lambda_geo / lambda = 0.500023`，`lambda_topo / lambda = 0.500005`，`lambda_topo / lambda_geo = 0.999964`
- `mu = -0.10 eV`：`lambda_geo / lambda = 0.500148`，`lambda_topo / lambda = 0.499847`，`lambda_topo / lambda_geo = 0.999400`

因此当前代码已经具备下面这个解释力：

- `lambda_geo` 是量子度规控制的实际几何贡献
- `lambda_topo` 是由 `K/K'` Dirac 绕数给出的下界
- `topological_fraction = lambda_topo / lambda_total`
- `topological_fraction_of_geo = lambda_topo / lambda_geo`，当前占比图使用这个量
- 在石墨烯低能极限中，这个下界数值上确实逼近实际值

## 每次回传内容

- 运行命令
- 输入参数文件路径
- 关键数值
- 生成图片路径
- 报错或异常现象
