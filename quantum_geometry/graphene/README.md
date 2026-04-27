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

## 每次回传内容

- 运行命令
- 输入参数文件路径
- 关键数值
- 生成图片路径
- 报错或异常现象
