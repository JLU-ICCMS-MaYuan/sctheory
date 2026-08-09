# MgB2 复现说明

本目录负责复现主文 `Fig. 3`、`SI FIG. 4-10` 和 `TABLE II`。

## 目标

第一轮只追四件事：

- `pi` 和 `sigma` 两条通道都能单独计算
- `sigma` 通道主导总 `lambda`
- `E2g(Gamma)` 是主导模式
- 总几何贡献接近文中约 90%

## 子目录

- `tb_pi/`
  - `pz / pi` 通道模型、量子几何和 `lambda_pi` 分解
- `tb_sigma/`
  - 完整 `px,py / sigma` 紧束缚
  - OFSM
  - `lambda_sigma` 分解
- `abinitio/qe_epw/`
  - 原文主路线：`QE + DFPT + Wannier90 + EPW`
- `abinitio/qe_direct_check/`
  - 原文校验路线：`QE/DFPT` 直接计算
- `plots/`
  - 主文 `Fig. 3`、`SI FIG. 4-10`
- `runs/`
  - 每次运行的参数、结果摘要、日志

## 模型级硬约束

- `pi` 和 `sigma` 必须分开推进
- `sigma` 通道禁止继续使用 A 点附近抛物线近似作为最终结果
- `Research/chemical-bonding/reference-notes/reproduce_mgb2_tba.py` 只能作为历史参考，不能直接当最终 `sigma` 模型

## 原文关键数值目标

- `lambda ≈ 0.78`
- `lambda_ab_initio ≈ 0.67`
- `lambda_pi ≈ 0.16`
- `lambda_sigma ≈ 0.62`
- `lambda_geo ≈ 0.71`
- `lambda_topo ≈ 0.32`
- 几何贡献约 `91.7%`

## 第一性原理主路线

- `QE + DFPT + Wannier90 + EPW`

关键参数：

- `12x12x12` 电子网格
- `6x6x6` 声子网格
- `60x60x60` 插值 `k/q` 网格
- 5 个 Wannier 初始中心
  - 2 个 B `pz`
  - 3 个 `s-like` center

## 第一性原理校验路线

- `QE/DFPT` 直接计算，不经过 Wannier

原文一致性检查目标：

- `lambda_I = 0.67`
- `lambda_II = 0.645`
- 相对差约 `3.7%`

## 推荐执行顺序

1. 在 `tb_pi/` 中建立 `pz` 模型，完成 `pi` 通道
2. 在 `tb_sigma/` 中建立完整 `px,py` 紧束缚，完成 `sigma` 通道
3. 合并总 `lambda`、`lambda_geo`、`lambda_topo`
4. 在 `abinitio/qe_epw/` 中准备原文主路线输入
5. 需要时在 `abinitio/qe_direct_check/` 中补校验路线

## 图号对应

- `SI FIG. 4`：结构与 `E2g` 模式
- `SI FIG. 5`：完整与简化模型 band structure、Fermi surface、`DOS(E_F)`
- `SI FIG. 6-7`：`sigma` 通道拓扑与子空间
- `SI FIG. 8a`：`sigma` 通道 `Gamma_nm` 模型/ab initio 对比
- `SI FIG. 8b`：`pi` 通道 `Gamma_nm` 模型/ab initio 对比
- `SI FIG. 8c`：`alpha^2F(omega)` 与 `E2g(Gamma)` 峰值
- `SI FIG. 9`：frozen phonon 证明 two-center approximation 合理
- `SI FIG. 10`：两条 ab initio 路线一致性

## 每次回传内容

- 运行命令
- 输入参数文件路径
- 关键数值
- 生成图片路径
- 报错或异常现象
