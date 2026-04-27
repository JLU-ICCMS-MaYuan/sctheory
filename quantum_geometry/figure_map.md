# 图表任务映射

本文档把主文和 SI 中的全部主图、补充图与表格映射到具体复现任务。执行时以“图号 -> 物理量 -> 计算步骤 -> 输出文件”为主线，不按章节泛泛推进。

## 主文

| 图号 | 类型 | 图意 | 需要复现的量 | 主要来源 |
| --- | --- | --- | --- | --- |
| Fig. 1 | 概念图 | 量子几何、拓扑下界与 EPC 的关系 | 无数值复现，需文档解释 `lambda_E / lambda_geo / lambda_topo` 的物理含义 | `README.md` 与后续推导说明 |
| Fig. 2a | 模型/ab initio 对比图 | 石墨烯 `lambda_model(mu)` 与 `lambda_ab_initio(mu)` 的趋势对比 | `D(mu)`、`lambda_total(mu)` | `graphene/tb/` + `graphene/abinitio/` |
| Fig. 2b | 贡献分解图 | 石墨烯几何贡献比例约 50% | `lambda_geo(mu) / lambda(mu)` | `graphene/tb/` |
| Fig. 2c | 拓扑下界图 | `lambda_topo` 接近 `lambda_geo`，Dirac 点附近趋于重合 | `lambda_topo(mu)`、`lambda_geo(mu)` | `graphene/tb/` |
| Fig. 3 | 总结图 | MgB2 的结构、轨道、费米面和 `pi/sigma` 分解结论 | `lambda_pi`、`lambda_sigma`、总 `lambda`、主导模式说明 | `mgb2/tb_pi/`、`mgb2/tb_sigma/`、`mgb2/abinitio/` |

## SI：石墨烯

| 图号/表号 | 类型 | 图意 | 需要复现的量 | 主要来源 |
| --- | --- | --- | --- | --- |
| FIG. 1 | ab initio 基线图 | 石墨烯 `ab initio` 下 `sqrt(<omega^2>)`、`lambda` 或相关扫描基线 | 化学势扫描下的基础量 | `graphene/abinitio/jdftx/` 或 `graphene/abinitio/qe_equiv/` |
| FIG. 2 | EPC 矩阵元对比图 | `F_tau_i(K, K+q)` 的模型与 `ab initio` 对比 | 逐分量 EPC 矩阵元沿路径变化 | `graphene/tb/` + `graphene/abinitio/` |
| FIG. 3 | 误差/贡献图 | 石墨烯模型与 `ab initio` 的误差与贡献对比 | 模型误差、关键点比较 | `graphene/tb/` + `graphene/abinitio/` |
| TABLE I | 参数表 | 石墨烯模型参数 | `a=2.46 A`、`t=-2.751 eV`、`gamma=-7.308 a^-2`、`hbar*sqrt(<omega^2>)=0.1615 eV`、`hbar*omega_E2g(Gamma)=0.1935 eV`、`hbar*omega_A1'(K)=0.1622 eV` | 固化到脚本/配置 |

## SI：MgB2

| 图号/表号 | 类型 | 图意 | 需要复现的量 | 主要来源 |
| --- | --- | --- | --- | --- |
| FIG. 4 | 结构/模式图 | MgB2 结构与 `E2g` 振动模式 | 结构图、模式说明 | `mgb2/abinitio/` 和文档说明 |
| FIG. 5 | 模型图 | 完整与简化 MgB2 紧束缚能带、费米面比较 | band structure、Fermi surface、`DOS(E_F)` | `mgb2/tb_pi/` + `mgb2/tb_sigma/` |
| FIG. 6 | 拓扑证明图 | `spxpy` 子空间、Wilson loop、Euler/拓扑结构 | `kz=0` 子空间、Wilson loop、拓扑指标 | `mgb2/tb_sigma/` |
| FIG. 7 | 子空间分析图 | `E1u@1a` 相关成分分析 | 轨道投影、子空间分解 | `mgb2/tb_sigma/` |
| FIG. 8a | 模型/ab initio 对比图 | `sigma` 通道 `Gamma_nm(Gamma, Gamma+(0,0,kz))` 对比 | `sigma` 通道 EPC 线宽相关量 | `mgb2/tb_sigma/` + `mgb2/abinitio/qe_epw/` |
| FIG. 8b | 模型/ab initio 对比图 | `pi` 通道 `Gamma_nm(K, K+(0,0,kz))` 对比 | `pi` 通道 EPC 线宽相关量 | `mgb2/tb_pi/` + `mgb2/abinitio/qe_epw/` |
| FIG. 8c | Eliashberg 图 | `alpha^2F(omega)` 峰值靠近 `E2g(Gamma)` | `alpha^2F(omega)`、`omega_E2g(Gamma)` | `mgb2/abinitio/qe_epw/` |
| FIG. 9 | frozen phonon 证据图 | `sigma` 通道主要来自键长变化，不是 onsite EPC | frozen phonon 畸变下 bonding/antibonding 能级移动 | `mgb2/tb_sigma/` 或 `mgb2/abinitio/` |
| FIG. 10 | 双路线一致性图 | MgB2 两条 `ab initio` 方法得到的 `Gamma_nm` 接近 | `QE+EPW` 与 `QE direct` 的数值比较 | `mgb2/abinitio/qe_epw/` + `mgb2/abinitio/qe_direct_check/` |
| TABLE II | 结果表 | MgB2 各项 `lambda` 的分解 | `lambda≈0.78`、`lambda_ab_initio≈0.67`、`lambda_pi≈0.16`、`lambda_sigma≈0.62`、`lambda_geo≈0.71`、`lambda_topo≈0.32`、几何贡献 `≈91.7%` | `mgb2/tb_pi/` + `mgb2/tb_sigma/` + `mgb2/abinitio/` |

## 阶段与图号对应

### 阶段 A：石墨烯

- 先完成主文 `Fig. 2a`、`Fig. 2b`、`Fig. 2c`
- 同步准备 `SI FIG. 1`、`SI FIG. 2`、`SI FIG. 3`
- 固化 `TABLE I` 参数

### 阶段 B：MgB2

- 先完成 `SI FIG. 4`、`SI FIG. 5`
- 再拆成 `pi` 和 `sigma` 两条线推进 `SI FIG. 6-9`
- 最后汇总成主文 `Fig. 3`
- 需要时补 `SI FIG. 10`

## 当前默认软件链

- 石墨烯原文同构：`JDFTx + frozen phonon + Wannier`
- 石墨烯等价路线：`QE + frozen phonon + Wannier90`
- MgB2 主路线：`QE + DFPT + Wannier90 + EPW`
- MgB2 校验路线：`QE/DFPT` 直接计算
