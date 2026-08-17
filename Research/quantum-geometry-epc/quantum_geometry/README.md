# 量子几何论文复现总览

本目录用于复现以下两篇文献中关于石墨烯和 MgB2 的核心结果：

- `Research/quantum-geometry-epc/quantum_geometry/reference/2024 EPC-量子几何.pdf`
- `Research/quantum-geometry-epc/quantum_geometry/reference/2024 EPC-量子几何-SI.pdf`

本复现按“先模型级，再第一性原理对齐”的顺序推进。第一轮目标不是一次跑完整篇文章，而是先拿下两条核心结论：

- 石墨烯：`lambda_geo / lambda` 在小掺杂极限接近 50%
- MgB2：`sigma` 通道主导，总几何贡献接近 90%

## 目录结构

- `figure_map.md`：主文与 SI 全部图表的任务映射
- `plan.md`：历史规划记录，保留，不作为执行入口
- `graphene/`：石墨烯复现
- `mgb2/`：MgB2 复现

### `graphene/`

- `tb/`：紧束缚、量子几何、`lambda_E / lambda_geo / lambda_topo` 分解脚本
- `abinitio/jdftx/`：原文同构路线，`JDFTx + frozen phonon + Wannier`
- `abinitio/qe_equiv/`：等价路线，`QE + frozen phonon + Wannier90`
- `plots/`：主文 `Fig. 2` 和对应 SI 图
- `runs/`：每次运行的参数、结果摘要、日志

### `mgb2/`

- `tb_pi/`：`pz / pi` 通道模型与图
- `tb_sigma/`：`px,py / sigma` 通道模型、OFSM 与 EPC 分解
- `abinitio/qe_epw/`：原文主路线，`QE + DFPT + Wannier90 + EPW`
- `abinitio/qe_direct_check/`：原文校验路线，`QE/DFPT` 直接算
- `plots/`：主文 `Fig. 3` 和对应 SI 图
- `runs/`：每次运行的参数、结果摘要、日志

## 软件链选择

### 石墨烯

原文不是 `QE + EPW`，而是：

- `JDFTx + frozen phonon + Wannier`

原因是石墨烯这里重点是两带 Dirac 低能模型、局域 EPC 矩阵元和几何分解，不需要先上三维多带体系那种致密 `k/q` 插值流程。

如果你为了统一工作流不想用 `JDFTx`，可以改为：

- `QE + frozen phonon + Wannier90`

但必须在记录中注明：

- “这是等价复现，不是原文同软件复现”

### MgB2

原文主路线是：

- `QE + DFPT + Wannier90 + EPW`

原因是 MgB2 是三维多费米面体系，需要在致密 `k/q` 网格上稳定计算 `alpha^2F(omega)`、`lambda` 和 Fermi-surface average。原文还用了：

- `QE/DFPT` 直接计算

作为不经 Wannier 的数值校验路线。

## 推荐执行顺序

1. 完成 `figure_map.md` 对照，明确每张图的输出量
2. 跑石墨烯模型级结果
3. 跑石墨烯第一性原理基线
4. 跑 MgB2 `pi` 通道
5. 跑 MgB2 `sigma` 通道
6. 合并 MgB2 总 `lambda`
7. 跑 MgB2 第一性原理主路线
8. 需要时补 MgB2 双路线一致性检查

## 每次回传内容

每次只回传以下五类内容：

- 运行命令
- 输入参数文件路径
- 关键数值
- 生成图片路径
- 报错或异常现象
