# GEMINI 项目：高压氢化物超导研究环境

本项目是一个专注于高压氢化物超导体（如 $LaH_{10}$, $CeH_{10}$, $H_3S$）理论分析的研究平台。它结合了密度泛函理论 (DFT)、Eliashberg 强耦合超导理论以及群论（用于对称性分析和化学键研究）。

## 项目概览

主要研究目标是探索氢化物超导的微观机制，重点关注：
- **对比研究**：通过轨道反转 (Orbital Inversion) 和电声耦合 (EPC) 机制，解释 $LaH_{10}$ 与 $CeH_{10}$ 之间巨大的临界温度 ($T_c$) 差异。
- **光谱分析**：求解 Eliashberg 方程，计算 $T_c$、超导能隙及质量重整化因子。
- **对称性与成键**：利用 `GTPack` (Mathematica) 和晶体轨道哈密顿布居 (COHP) 分析成键环境。

## 核心物理洞察 (QA 汇总)

为解决研究过程中的常见困惑，本项目强调以下关键物理理解：

1.  **$N(0)$ 与超导强度的关系**：
    - **纠正误区**：$N(0)$ 并不在 $\alpha^2F(\omega)$ 的分母导致反比关系。
    - **正确结论**：$\alpha^2F \propto N(0) \times \langle|g|^2\rangle$。因此，费米面态密度 $N(0)$ 越大，电声耦合常数 $\lambda$ 越大，超导越强。
2.  **质量重整化因子 $Z$**：
    - $Z(0) \approx m^*/m$。$Z=1$ 对应自由费米气体。在强耦合系统中，$Z$ 可达 2-3，表示电子受声子云包裹导致有效质量显著增加。
3.  **能隙函数 $\Delta$ 的实部与虚部**：
    - 在 **Matsubara 虚频轴**上，$\Delta(i\omega_n)$ 是实数，适合数值自洽求解。
    - 在 **实频轴**上，$\Delta(\omega+i0^+)$ 是复数。实部对应能隙大小，虚部对应准粒子寿命。
4.  **轨道反转假说**：
    - 在 $CeH_{10}$ 中，Ce 的 4f 轨道可能改变 H 笼子的 $\sigma/\pi$ 键能级顺序。$\sigma$ 键对声子敏感（强 EPC），而 $\pi$ 键对声子钝化（弱 EPC）。

## 技术栈

- **模拟计算**：VASP (5.4+), Quantum ESPRESSO (6.7+), EPW, LOBSTER。
- **数据分析**：Python (3.7+) - `numpy`, `scipy`, `matplotlib`, `pymatgen`。
- **对称性分析**：Mathematica - `GTPack` (v1.4)。
- **可视化**：VESTA, Phonopy, FermiSurfer。

## 关键工作流

### 1. Eliashberg 方程求解器
位于 `EliashbergEquation/`，用于从谱函数 $\alpha^2F(\omega)$ 计算超导性质。
- **输入**：`ALPHA2F.OUT` (谱函数), `INPUT` (包含 $\mu^*$ 和温度步长)。
- **运行**：`python eliashberg_solver.py --input INPUT --alpha2f ALPHA2F.OUT`
- **输出**：$T_c$、McMillan 参数、随频率变化的能隙 $\Delta(\omega)$。

### 2. $LaH_{10}/CeH_{10}$ 对比分析
遵循 `docs/04_LaH10_CeH10_research_plan.md` 中的研究计划：
1.  **磁性测试**：确定 Ce 4f 电子的局域性/磁性 (`01_magnetic_test`)。
2.  **能带分析**：识别费米面附近的轨道反转及态密度特征 (`04_band_analysis`)。
3.  **冻结声子**：通过监测原子位移导致的能带分裂，直接探测 EPC 强度 (`06_frozen_phonon`)。

### 3. 对称性分析 (Mathematica)
使用 `GTPack` 进行：
- 构建紧束缚 (Tight-Binding) 模型。
- 轨道不可约表示分解（如 $D_{6h}$ 群下 f 轨道的分解）。
- 对称性适配波函数 (SAWPs) 的构造。

## 目录结构

- `docs/`：理论指南与研究计划。
  - `QUICKSTART.md`：$LaH_{10}/CeH_{10}$ 工作流简要步骤。
  - `03_theory_comprehensive.md`：Eliashberg 理论与 EPC 深层解析。
  - `04_LaH10_CeH10_research_plan.md`：当前研究的路线图（核心文档）。
- `EliashbergEquation/`：核心 Python 求解器及示例 I/O。
- `scripts/`：特定任务的分析脚本。
- `reference/`：支持研究假设的关键文献（PDF 及提取文本）。
- `MH10/`, `H3S/`, `MH9/`：特定体系的计算数据与 Mathematica 笔记本。

## 开发规范

- **数据追踪**：分析脚本应将日志和绘图输出到各自的子目录中。
- **代码风格**：Python 脚本需包含类型提示 (Type Hints) 和文档字符串。
- **解析延拓**：使用 Padé 近似从虚频数据获取实轴物理量，确保使用递归算法以保证数值稳定性。

## 核心文件
- `EliashbergEquation/eliashberg_solver.py`：主要的数值计算工具。
- `docs/04_LaH10_CeH10_research_plan.md`：当前研究方向的“真理来源”。
- `request.md`：记录 AI 交互历史与关键问题的解决记录。