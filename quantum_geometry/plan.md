# 量子几何论文在本仓库中的分阶段复现计划

## 摘要

目标是在仓库根目录下的空目录 quantum_geometry/ 中，建立一套“你计算、我分析、再推进下一步”的渐进式复现流程.

覆盖 reference/2024 EPC-量子几何.pdf 与  
reference/2024 EPC-量子几何-SI.pdf 里关于石墨烯和 MgB2 的核心结果。

默认采用“先模型级复现，后 ab initio 对齐”的顺序，因为仓库当前已有 pythtb、HamEPC、reference/reproduce_mgb2_tba.py 等基础，但没有现成的 MgB2/石墨烯 DFPT  
原始输入。

## 实施方案

### 1. 目录与产物组织(这个文章里面有多少个图片,每个图片什么意思, 我应该做什么复现出来, 这是你需要在目标阶段, 明确告诉我每一个阶段你要做的事情如何和每一个图片对应起来)

在 quantum_geometry/ 下建立两条主线：

- quantum_geometry/graphene/

  - tb/：紧束缚与量子几何脚本
  - plots/：带结构、Berry curvature、量子度量、EPC 分解图
  - runs/：每次运行的参数与原始输出
  - README.md：每一步要跑什么、预期结果是什么

- quantum_geometry/mgb2/

  - tb_pi/：MgB2 的 pz/π 部分
  - tb_sigma/：MgB2 的 px,py/σ 部分
  - abinitio/：QE/Wannier90/EPW 输入模板与运行记录(这一部分,我来完成我会在

    quantum_geometry/graphene, quantum_geometry/mgb2这两个路径里面分别做这两个第一性原里计算的, 分别用qe和wannier90来进行)

  - plots/、runs/、README.md 同上


### 2. 复现顺序与每一步要算什么

#### 阶段 A：石墨烯模型级复现，先跑最稳的结果

先做石墨烯，因为仓库里 external/pythtb 已有石墨烯模型和量子几何 API，可最快建立“量子几何影响 EPC”的数值直觉。

这一阶段要实现并运行：

1. 石墨烯紧束缚基线

   - 计算带结构，确认 Dirac 点位置与能带形状
   - 输出 graphene_bands.*

2. 石墨烯量子几何

   - 计算 Berry curvature、quantum metric、必要时 Fubini-Study metric
   - 重点看 K/K' 附近几何量增强
   - 输出 graphene_metric_map.*、graphene_curvature_map.*

3. 石墨烯 EPC 的 energetic / geometric 分解

   - 按主文和 SI 的公式，把 EPC 常数 λ 分成普通能量项与几何项(我应该用什么计算EPC, qe+EPW吗?文章中是这么搞的吗)
   - 先做文中的解析/紧束缚近似，不直接追求 ab initio 数值
   - 输出 lambda_total、lambda_geo、lambda_energetic、lambda_topo

4. 石墨烯检查点

   - 你把数值结果和图发给我
   - 我判断是否已经达到“几何贡献显著、趋势与文中一致”的标准
   - 通过后再进入 MgB2


#### 阶段 B：MgB2 模型级复现，分 π 和 σ 两部分

MgB2 不能直接沿用当前 reference/reproduce_mgb2_tba.py 的简化版本，因为那份笔记已经明确记录：σ 带目前还是抛物线近似，会导致 DOS 与耦合常数失真。这里要改  
成更接近文中 SI 的完整紧束缚工作流。

这一阶段要实现并运行：

1. π 通道

   - 建立 pz 两带模型
   - 计算带结构、费米面附近 DOS、量子几何量
   - 做 lambda_pi 的 energetic / geometric / topological 分解
   - 输出 mgb2_pi_*

2. σ 通道

   - 将现有 MgB2 脚本改为完整的 px,py 紧束缚哈密顿量，不再用 A 点附近抛物线替代
   - 计算 Γ-A 附近主导 σ 态、对应几何量和 EPC 分解
   - 输出 mgb2_sigma_*

3. 总 λ 合并

   - 汇总 lambda_pi、lambda_sigma
   - 检查几何贡献是否接近主文“MgB2 几何贡献占总 λ 的大头，约 90%”这一结论的量级与趋势
   - 不要求第一轮就完全数值对齐，但要先保证分解逻辑和主导通道正确


#### 阶段 C：ab initio 对齐准备，只搭可运行骨架

你已说明外部依赖由你安装，因此这一阶段只负责把工作流设计完整，不假设当前机器已经全配好。

这一阶段要补：

1. 石墨烯 ab initio 输入模板

   - 依据 SI：JDFTx 风格参数为参考，但在本仓库里优先准备等价的可执行说明和输入模板
   - 记录 24×24×1 电子网格、6×6×1 超胞冻声子、Wannier 试探中心等关键参数

2. MgB2 ab initio 输入模板

   - 依据 SI：QE + Wannier90 + EPW 路线
   - 记录 12×12×12 电子网格、6×6×6 声子网格、60×60×60 插值网格、五个 Wannier 初始中心等关键参数

3. 对齐策略

   - 先用模型级结果确认“几何项为何大”
   - 再用 ab initio 检查“数量级是否接近文中”


### 3. 你每一步需要回传给我的内容

每跑完一个阶段，只回传这四类信息，我就能判断下一步怎么改：

- 运行命令
- 关键数值
- 报错日志或不收敛现象

1. 先完成石墨烯模型级
2. 我分析
3. 再推进 MgB2 π
4. 我分析
5. 再推进 MgB2 σ
6. 我分析
7. 最后再搭 ab initio 对齐

## 测试与验收

- 石墨烯

  - 能带应出现正确 Dirac 锥
  - 量子几何量在 K/K' 附近出现显著结构
  - λ 分解后几何项占比应与文中趋势一致，至少定性正确

- MgB2

  - π、σ 两类通道要分别可计算、可出图、可汇总
  - σ 通道必须基于完整紧束缚，而不是当前脚本中的抛物线近似
  - 总 λ 的主导贡献应来自文中强调的几何相关部分，至少先达到定性复现

- 工作流

  - 每个阶段都应能通过单条命令重复运行
  - 每次运行都在 runs/ 下留下参数和结果快照，避免后续分析时丢上下文


## 关键接口与默认假设

- 默认新代码落在 quantum_geometry/，不污染 reference/ 和已有材料目录
- 默认先复现模型级结果，再做 ab initio
- 默认你负责安装和运行外部依赖，我负责定义每一步要算什么、怎样判断结果是否对
- 默认第一轮不追求一步到位重现论文所有图，而是优先拿下这几个核心量：

  - 石墨烯：lambda_total、lambda_geo、量子度量分布

- MgB2：lambda_pi、lambda_sigma、总 λ 的几何贡献比例



















# 量子几何论文复现计划 2.0：按图对应的分阶段路线

## 摘要

这篇文章当前可明确对应到：

- 主文 3 张主图：Fig. 1 到 Fig. 3
- SI 10 张补充图：FIG. 1 到 FIG. 10
- 关键数值表至少 2 个：

  - 石墨烯参数表 TABLE I
  - MgB2 结果表 TABLE II


复现目标不是“泛泛做石墨烯/MgB2”，而是把每一张图背后的物理对象、数值流程、以及仓库中的实现位置逐一对应起来。  
默认代码与结果都放在 quantum_geometry/graphene/ 和 quantum_geometry/mgb2/ 下。

关于你问的 EPC 路线，结论先说清楚：

- 论文里石墨烯的 ab initio 不是 QE+EPW
- 石墨烯原文路线是：JDFTx + frozen phonon + Wannier
- MgB2 原文主路线是：QE + Wannier90 + EPW
- MgB2 还有一条原文校验路线：QE/DFPT 直接算，不经过 Wannier

如果你的目标是“尽量忠实复现原文软件链”，石墨烯应优先按 JDFTx；如果你的目标是“在你现有环境下尽量复现物理结果”，石墨烯也可以用 QE + frozen phonon +  
Wannier90 做等价替代，但这属于“等价实现”，不是“原文同构实现”。

## 图谱与阶段对应

### 阶段 0：读图建账，不做数值拟合

这一阶段只建立图和任务的映射表，输出一个总索引文档。

对应图片：

- 主文 Fig. 1

  - 含义：量子几何、拓扑下界、EPC 的概念图
  - 不需要数值复现
  - 需要在文档里解释它对应后面所有公式分解的物理意义

- 主文 Fig. 2

  - 含义：石墨烯中 λ、λ_geo、λ_topo 随化学势变化，以及模型和 ab initio 的对比
  - 这是石墨烯阶段的核心验收图

- 主文 Fig. 3

  - 含义：MgB2 的结构、主导轨道/费米面、以及 π/σ 分解后的 EPC 结论
  - 这是 MgB2 阶段的核心验收图

- SI FIG. 1 到 FIG. 10

  - 作为主文三张图的展开支撑
  - 每一张都要被挂到某个复现步骤上，不能悬空


阶段 0 的产物：

- quantum_geometry/README.md
- quantum_geometry/figure_map.md

figure_map.md 必须逐条写清楚：

- 图号
- 图在讲什么
- 是“概念图 / 模型图 / ab initio 对比图 / 拟合支撑图 / 拓扑证明图”
- 你要计算什么量才能复现它
- 该量由哪个脚本或哪条外部流程产生

### 阶段 A：石墨烯模型级复现

目标是先把主文 Fig. 2 和石墨烯相关 SI 图做出来。

对应图片与任务：

- 主文 Fig. 2a

  - 目标：复现 λ_model(μ) 与 λ_ab_initio(μ) 的趋势对比
  - 你要算：

    - 石墨烯能带
    - 化学势扫描下的 D(μ)
    - λ_total(μ)
    - λ_E(μ)
    - λ_geo(μ)
    - λ_topo(μ)


- 主文 Fig. 2b

  - 目标：证明石墨烯 λ_geo / λ ≈ 50%
  - 你要算：

    - 几何贡献比例随 μ 的变化


- 主文 Fig. 2c

  - 目标：证明 λ_topo 接近 λ_geo，尤其在 Dirac 点附近趋于重合
  - 你要算：

    - 拓扑下界与几何项的对比曲线


- SI FIG. 1

  - 含义：石墨烯 ab initio 下 ⟨ω²⟩、λ 或相关量随化学势的基础数值
  - 作用：给主文 Fig. 2 的参数输入和趋势基线
  - 任务：你做 ab initio 时要把这类基础扫描先跑出来

- SI FIG. 2

  - 含义：Fτi(K, K+q) 的模型与 ab initio 对比
  - 作用：验证石墨烯 EPC 近似模型是可信的
  - 任务：要输出 EPC 矩阵元沿高对称路径的逐分量对比

- SI FIG. 3

  - 含义：石墨烯更进一步的 EPC/贡献对比与模型误差说明
  - 作用：说明为什么主文 Fig. 2 的 λ 拟合成立
  - 任务：保留模型与 ab initio 的误差曲线或关键点比较

- SI TABLE I

  - 含义：石墨烯模型参数表
  - 任务：在你的脚本里固定记录

    - a = 2.46 Å
    - t = -2.751 eV
    - γ = -7.308 a^{-2}
    - ℏ√⟨ω²⟩ = 0.1615 eV
    - ℏω_E2g(Γ) = 0.1935 eV
    - ℏω_A1'(K) = 0.1622 eV



这一阶段在仓库中要有：

- quantum_geometry/graphene/tb/

  - 带结构
  - quantum metric / FSM
  - λ_E / λ_geo / λ_topo 分解

- quantum_geometry/graphene/abinitio/

  - 你来放 JDFTx 或 QE+frozen-phonon+Wannier90 输入

- quantum_geometry/graphene/plots/

  - 对应 Fig. 2 和 SI FIG. 1-3

- quantum_geometry/graphene/runs/

  - 每次跑的参数和原始输出


石墨烯阶段验收标准：

- 出现正确 Dirac 锥
- K/K' 附近量子几何量增强
- λ_geo / λ 在小掺杂极限接近 50%
- λ_topo 与 λ_geo 接近，趋势一致
- EPC 矩阵元模型对比图可以解释为何主文近似有效

### 阶段 B：MgB2 模型级复现

目标是复现主文 Fig. 3，并用 SI FIG. 4-9 把 MgB2 的结构、能带、拓扑、EPC 和主导模式全部拆清楚。

对应图片与任务：

- SI FIG. 4

  - 含义：MgB2 晶体结构与 E2g 振动模式
  - 任务：建立结构输入和声子模式解释图
  - 这是 MgB2 一切后续建模的出发点

- SI FIG. 5

  - 含义：完整与简化 MgB2 紧束缚能带、费米面比较
  - 任务：先复现完整低能紧束缚，再复现去掉 Mg 后的简化模型
  - 你要输出：

    - band structure
    - Fermi surface
    - DOS at E_F


- SI FIG. 6

  - 含义：spxpy 子空间、Wilson loop、Euler/拓扑结构
  - 任务：这是 σ 通道几何项的拓扑根源
  - 你要算：

    - kz=0 面上的目标子空间
    - Wilson loop
    - 相关拓扑指标


- SI FIG. 7

  - 含义：E1u@1a 与相关子空间的成分分析
  - 任务：补足 σ 有效模型中的轨道选择规则
  - 你要输出：

    - 轨道投影着色图
    - 子空间分解结果


- SI FIG. 8a

  - 含义：MgB2 σ 通道 Γnm(Γ, Γ+(0,0,kz)) 模型与 ab initio 对比
  - 任务：这是 σ 通道 EPC 模型有效性的核心数值图

- SI FIG. 8b

  - 含义：MgB2 π 通道 Γnm(K, K+(0,0,kz)) 模型与 ab initio 对比
  - 任务：这是 π 通道 EPC 模型有效性的核心数值图

- SI FIG. 8c

  - 含义：MgB2 的 Eliashberg 函数 α²F(ω)，并显示峰值靠近 E2g(Γ)
  - 任务：验证主导声子模式确实是 E2g

- SI TABLE II

  - 含义：MgB2 各项 λ 的分解结果
  - 你最终要复现的关键数值目标是：

    - λ ≈ 0.78，ab initio ≈ 0.67
    - λπ ≈ 0.16
    - λσ ≈ 0.62
    - λ_geo ≈ 0.71
    - λ_topo ≈ 0.32
    - 几何贡献约 91.7%


- SI FIG. 9

  - 含义：冻结 E2g 声子后，证明 σ 通道主要来自键长变化，不是 onsite EPC
  - 任务：这是“为什么 two-center approximation 合理”的证据图
  - 你要做：

    - frozen phonon 畸变
    - 比较 bonding/antibonding 的能级移动方式



这一阶段必须拆成 π 和 σ 两条线，不能混在一起。

#### B1：MgB2 π 通道

要做：

- 建立 pz 模型
- 计算带结构和费米线/费米面
- 计算量子几何量
- 计算 λπ,E、λπ,geo、λπ,topo
- 与 SI FIG. 8b 对应验证

#### B2：MgB2 σ 通道

要做：

- 把 reference/reproduce_mgb2_tba.py 升级成完整 px,py 紧束缚，不允许继续用抛物线近似充当最终结果
- 计算 Γ-A 附近有效模型
- 计算 OFSM
- 计算 λσ,E、λσ,geo、λσ,topo
- 与 SI FIG. 6-9 对应验证

#### B3：MgB2 总结图

最后合并：

- λ = λπ + λσ
- λ_geo = λπ,geo + λσ,geo
- λ_topo = λπ,topo + λσ,topo

并生成一张“对照主文 Fig. 3 的总结图”。

MgB2 阶段验收标准：

- σ、π 两条通道分别可算
- σ 通道明显主导总 λ
- λσ,E 很小，λσ,geo 很大
- E2g(Γ) 是主导模式
- 总几何贡献接近文中 ~90%

### 阶段 C：第一性原理对齐

这一阶段的目标不是先做新物理，而是给阶段 A/B 的模型结论做数值背书。

#### C1：石墨烯 ab initio

如果你坚持严格贴论文，建议：

- JDFTx + frozen phonon + Wannier

如果你希望工作流统一到你的仓库习惯，可以接受：

- QE + frozen phonon + Wannier90

但要在文档里明确标注：

- “这是等价复现，不是原文同软件复现”

石墨烯必须记录的原文参数：

- 24×24×1 电子网格
- 6×6×1 超胞冻声子
- 7 个 trial Wannier center

  - 2 个 C pz
  - 3 个键中心 s-like
  - 2 个六角中心上下 s-like


石墨烯这一阶段的目标图：

- SI FIG. 1-3
- 主文 Fig. 2a

#### C2：MgB2 ab initio

按原文主路线：

- QE + DFPT + Wannier90 + EPW

原文关键参数：

- 12×12×12 电子网格
- 6×6×6 声子网格
- 60×60×60 插值 k/q 网格
- 5 个 Wannier 初始中心

  - 2 个 B pz
  - 3 个 s-like center


原文还有校验路线：

- QE/DFPT 直接算，不经过 Wannier
- 结果：

  - λ_I = 0.67
  - λ_II = 0.645
  - 相对差约 3.7%


对应图片：

- SI FIG. 10

  - 含义：MgB2 两条 ab initio 方法得到的 Γnm 非常接近
  - 任务：如果你愿意做双路线，就拿它作为数值可靠性检验
  - 如果你只做一条路线，这张图可以暂时不作为第一阶段强制目标


## 每一步你需要回传给我的内容

每次只回传这五类，不要多也不要少：

- 关键数值

我会按这个顺序接你的结果：

1. 石墨烯模型级
2. 石墨烯 ab initio 基线
3. MgB2 π
4. MgB2 σ
5. MgB2 总 λ 汇总
6. MgB2 ab initio
7. 双路线一致性检查

## 默认实现约束

- 所有新实现只放在 quantum_geometry/
- 不修改 reference/ 中原始文献文件
- 允许你在

  - quantum_geometry/graphene/abinitio/
  - quantum_geometry/mgb2/abinitio/  
    里分别放外部第一性原理输入

- 石墨烯阶段先追主文 Fig. 2 和 SI FIG. 1-3
- MgB2 阶段先追主文 Fig. 3 和 SI FIG. 4-9
- SI FIG. 10 作为 MgB2 数值可靠性增强项，不是第一轮阻塞项

## 关键假设

- 你负责安装并运行外部依赖
- 我负责逐阶段告诉你“现在该算什么、输出什么、如何判断是否对”
- 第一轮不要求你一步跑完整篇文章
- 第一轮的成功标准是：

  - 石墨烯先拿下 50% 几何贡献
  - MgB2 先拿下 σ 主导和 ~90% 几何贡献的物理结论
