# Eliashberg超导理论：完整理论基础

本文档系统整理Eliashberg超导方程的完整理论框架，包括有限温度场论、Matsubara频率表示、解析延拓和第一性原理计算接口。

---

## 第一章：有限温度量子场论基础

### 1.1 配分函数与虚时间形式

**配分函数的普遍形式**：
$$
Z = \text{Tr}[e^{-\beta H}], \quad \beta = \frac{1}{k_B T}
$$

**为何使用迹**：
- 保证基组无关性
- 物理观测量不依赖于表示

**Wick旋转**：
$$
e^{-\beta H} = e^{-iHt}\bigg|_{t=i\beta\hbar} \quad \text{(实时间 → 虚时间)}
$$

### 1.2 Matsubara频率的起源

**边界条件**：

虚时间格林函数必须满足：
- **费米子**：$G(\tau + \beta) = -G(\tau)$ (反周期)
- **玻色子**：$D(\tau + \beta) = D(\tau)$ (周期)

**离散化频率**：

- **费米子**：$\omega_n = \pi k_B T(2n+1)$，$n \in \mathbb{Z}$
- **玻色子**：$\Omega_m = 2\pi k_B T \cdot m$，$m \in \mathbb{Z}$

**物理图像**：温度$T$像一个"尺子"，把虚时间切成$\sim k_B T$大小的片段。

### 1.3 费米子反周期性的严格推导

**关键恒等式**：
$$
e^{\beta H} \psi e^{-\beta H} = -\psi
$$

**证明**：源于费米子反对易关系$\{\psi, \psi^\dagger\} = 1$。

**推论**：
$$
\omega_{-n-1} = -\omega_n
$$

这是Matsubara频率折叠的核心。

---

## 第二章：Eliashberg方程理论

### 2.1 各向同性Eliashberg方程

**质量重整化方程**：
$$
Z(i\omega_n) = 1 + \frac{\pi T}{\omega_n}\sum_{m=-\infty}^{+\infty} \lambda(i\omega_n - i\omega_m) \frac{\omega_m Z_m}{\sqrt{\omega_m^2 + \Delta_m^2}}
$$

**能隙方程**：
$$
Z(i\omega_n)\Delta_n = \pi T \sum_{m=-\infty}^{+\infty} [\lambda(i\omega_n - i\omega_m) - \mu_c^*] \frac{\Delta_m}{\sqrt{\omega_m^2 + \Delta_m^2}}
$$

**物理意义**：
- $\Delta_n$：超导能隙（Cooper对结合能）
- $Z_n$：准粒子重整化因子（有效质量增强）
- $\lambda$：电子-声子耦合核
- $\mu_c^*$：库仑赝势（电子-电子排斥修正）

### 2.2 电子-声子耦合核

**定义**：
$$
\lambda(i\omega_n - i\omega_m) = 2\int_0^{\infty} d\Omega \frac{\Omega \alpha^2F(\Omega)}{\Omega^2 + (\omega_n - \omega_m)^2}
$$

**关键性质**：
1. **偶对称**：$\lambda(i\omega) = \lambda(-i\omega)$
2. **无奇点**：在虚轴上解析
3. **表格化**：离散化为$\Lambda_m$数组

### 2.3 Matsubara频率折叠

**原始求和**（对所有$m \in \mathbb{Z}$）→ **折叠求和**（只对$m \geq 0$）

**Z方程折叠形式**：
$$
Z_n = 1 + \frac{\pi T}{\omega_n}\sum_{m=0}^{N_\omega-1} [\Lambda_{n-m} - \Lambda_{n+m+1}] \frac{\omega_m Z_m}{R_m}
$$

**能隙方程折叠形式**：
$$
\Delta_n = \frac{\pi T}{Z_n} \left[\sum_{m=0}^{N_\omega-1} [\Lambda_{n-m} + \Lambda_{n+m+1}] \frac{\Delta_m Z_m}{R_m} - \Delta\mu\right]
$$

其中$R_m = |Z_m \omega_m| = \sqrt{(\omega_m^2 + \Delta_m^2)Z_m^2}$。

**关键理解**：
- 索引$l(n+m+1)$对应频率差$\omega_n + \omega_m$
- 利用$\omega_{-m-1} = -\omega_m$实现折叠
- **不是**"加减$2\pi k_B T$不变"！

---

## 第三章：实轴奇点与Padé延拓

### 3.1 为何在虚轴求解？

**实轴问题**：

实频格林函数：
$$
G^R(\omega) = \frac{1}{\omega - \epsilon + i0^+}
$$

- 有极点/δ函数奇点
- Cauchy主值积分数值困难
- 自洽迭代不稳定

**Matsubara优势**：

虚频格林函数：
$$
G(i\omega_n) = \frac{1}{i\omega_n - \epsilon}
$$

- 无奇点，光滑函数
- 离散求和，数值稳定
- 自洽迭代快速收敛

### 3.2 Padé解析延拓

**目的**：从Matsubara轴($i\omega_n$)延拓到实轴($\omega + i0^+$)

**Vidberg-Serene算法**：

1. 构建$g$矩阵：
   $$
   g_{i,j} = \frac{g_{i-1,i-1} - g_{i-1,j}}{(z_j - z_{i-1}) g_{i-1,j}}
   $$

2. 连分式展开：
   $$
   P(z) = g_{0,0} + \frac{(z - z_0) g_{1,1}}{1 + \frac{(z - z_1) g_{2,2}}{1 + \cdots}}
   $$

**精度**：对光滑函数，相对误差 < 10⁻⁶

---

## 第四章：McMillan-Allen-Dynes临界温度

### 4.1 电声耦合参数

**耦合常数**：
$$
\lambda = 2\int_0^{\omega_{\max}} \frac{\alpha^2 F(\omega)}{\omega} d\omega
$$

**特征频率**：
$$
\omega_{\log} = \exp\left[\frac{2}{\lambda}\int_0^{\omega_{\max}} \frac{\alpha^2 F(\omega)\ln\omega}{\omega} d\omega\right]
$$

$$
\omega_{\mathrm{rms}} = \sqrt{\frac{2}{\lambda}\int_0^{\omega_{\max}} \alpha^2 F(\omega) \omega d\omega}
$$

### 4.2 Allen-Dynes公式

**临界温度**：
$$
T_c = \frac{\omega_{\log}}{1.2 k_B} \exp\left[\frac{-1.04(1+\lambda)}{\lambda - \mu^* - 0.62\lambda\mu^*}\right] f_1 f_2
$$

**修正因子**：
$$
f_1 = \left[1 + \left(\frac{\lambda}{2.46(1+3.8\mu^*)}\right)^{3/2}\right]^{1/3}
$$

$$
f_2 = 1 + \frac{\lambda^2(\omega_{\mathrm{rms}}/\omega_{\log} - 1)}{\lambda^2 + [1.82(1+6.3\mu^*)(\omega_{\mathrm{rms}}/\omega_{\log})]^2}
$$

**适用范围**：$\lambda < 1.5$，$\mu^* \sim 0.1-0.15$

---

## 第五章：第一性原理计算接口

### 5.1 Eliashberg谱函数α²F(ω)

**定义**（见第六章详细推导）：
$$
\alpha^2F(\omega) = \frac{1}{N(0)} \sum_{\mathbf{k},\mathbf{q},\nu} |g_{mn\nu}(\mathbf{k},\mathbf{q})|^2 \delta(\epsilon_{m\mathbf{k}}) \delta(\epsilon_{n\mathbf{k}+\mathbf{q}}) \delta(\omega - \omega_{q\nu})
$$

**物理意义**：
- 电声耦合强度的频率分布
- 加权的声子态密度
- 归一化：$\int \alpha^2F(\omega) d\omega = \lambda/2$

### 5.2 从EPW提取数据

**三个层次**：

1. **各向同性参数**（已实现）：
   - 读取`ALPHA2F.OUT`
   - 计算$\lambda$、$\omega_{\log}$、$\omega_{\mathrm{rms}}$
   - 输入`eliashberg_solver.py`

2. **模式分解**（中等难度）：
   - 分解$\alpha^2 F(\omega) = \sum_\nu \alpha^2 F_\nu(\omega)$
   - 识别主导声子模式

3. **各向异性数据**（高难度）：
   - 提取$g_{mn\nu}(\mathbf{k},\mathbf{q})$矩阵元
   - 构建$\mathbf{k}$依赖的有效模型

**工具**：`extract_from_epw.py`实现层次1和层次2

---

## 第六章：α²F(ω)的物理意义与计算（新增）

### 6.1 定义与物理意义

**完整定义**：
$$
\alpha^2F(\omega) = \frac{1}{N(0)} \sum_{\mathbf{k},\mathbf{q}} \sum_{m,n} \sum_{\nu} |g_{mn\nu}(\mathbf{k},\mathbf{q})|^2 \delta(\epsilon_{m\mathbf{k}}) \delta(\epsilon_{n\mathbf{k}+\mathbf{q}}) \delta(\omega - \omega_{q\nu})
$$

其中：
- $g_{mn\nu}(\mathbf{k},\mathbf{q})$：电声耦合矩阵元
- $\epsilon_{m\mathbf{k}}$：电子能量（相对费米面）
- $\omega_{q\nu}$：声子频率（支标$\nu$）
- $N(0)$：费米面态密度

**物理解释**：

1. **电声耦合强度**：$|g_{mn\nu}|^2$描述电子通过吸收/发射声子的散射强度
2. **费米面约束**：两个δ函数确保只有费米面附近的电子参与
3. **声子能量依赖**：第三个δ函数选择频率为$\omega$的声子模式
4. **归一化**：$N(0)$确保结果与态密度无关

**与λ的关系**：
$$
\lambda = 2\int_0^{\infty} \frac{\alpha^2F(\omega)}{\omega} d\omega
$$

### 6.2 从Quantum ESPRESSO计算α²F(ω)

**完整工作流程**：

#### 步骤1：DFT自洽计算 (QE: pw.x)

```bash
pw.x < scf.in > scf.out
```

**关键设置**：
- 收敛精度：`conv_thr = 1.0d-12`
- k点网格：密集网格（如12×12×12）
- 交换关联泛函：通常使用PBE

#### 步骤2：声子计算 (QE: ph.x)

```bash
ph.x < ph.in > ph.out
```

**关键输入**：
```fortran
&inputph
  tr2_ph = 1.0d-14,      ! 声子收敛阈值
  ldisp = .true.,        ! 计算声子色散
  nq1 = 6, nq2 = 6, nq3 = 6,  ! q点网格
  fildyn = 'dyn',        ! 动力学矩阵文件前缀
/
```

**输出**：
- 动力学矩阵：`dyn0`, `dyn1`, ..., `dynN`
- 声子频率和本征矢

#### 步骤3：电声耦合插值 (EPW)

**EPW输入文件**（`epw.in`）：

```fortran
&inputepw
  ! 关键参数
  elph = .true.              ! 计算电声耦合
  epbwrite = .false.         ! 不写原始矩阵元（节省空间）
  epbread = .false.

  ! Wannier插值
  wannierize = .true.        ! 使用Wannier90
  nbndsub = 8                ! Wannier能带数

  ! k和q网格
  nkf1 = 40, nkf2 = 40, nkf3 = 40  ! 精细k网格
  nqf1 = 20, nqf2 = 20, nqf3 = 20  ! 精细q网格

  ! Eliashberg谱函数
  eliashberg = .true.        ! 计算α²F(ω)
  la2F = .true.              ! 输出α²F文件

  ! 频率范围
  nqstep = 500               ! 频率点数
  degaussw = 0.05            ! 展宽参数 (eV)
/
```

**关键过程**：

1. **粗网格计算**：
   - DFT和DFPT在粗k/q网格上计算
   - 得到$g_{mn\nu}(\mathbf{k},\mathbf{q})$在粗网格上的值

2. **Wannier插值**：
   - 将Bloch态投影到Wannier基组
   - 在Wannier表示下插值到任意k点

3. **精细网格重建**：
   - 在密集k/q网格上重建$g_{mn\nu}$
   - 计算α²F(ω)的频率依赖

**数学细节**：
$$
g_{mn\nu}(\mathbf{k},\mathbf{q}) = \langle \psi_{m\mathbf{k}} | \frac{\partial V}{\partial u_{q\nu}} | \psi_{n\mathbf{k}+\mathbf{q}} \rangle
$$

其中$u_{q\nu}$是声子位移模式。

#### 步骤4：α²F(ω)的数值计算

**离散化公式**：
$$
\alpha^2F(\omega) = \frac{1}{N(0)} \frac{1}{N_k N_q} \sum_{\mathbf{k},\mathbf{q},m,n,\nu} |g_{mn\nu}|^2 \delta(\epsilon_{m\mathbf{k}}) \delta(\epsilon_{n\mathbf{k}+\mathbf{q}}) \delta(\omega - \omega_{q\nu})
$$

**δ函数处理**（高斯展宽）：
$$
\delta(x) \to \frac{1}{\sigma\sqrt{2\pi}} e^{-x^2/(2\sigma^2)}
$$

**展宽参数选择**：
- 电子：`degaussq` ~ 0.05 eV（匹配费米面展宽）
- 声子：`degaussw` ~ 0.5-1.0 cm⁻¹

**输出文件**：`ALPHA2F.OUT`

```
# omega (Hartree)    alpha2F(omega)
0.000000000000      0.000000000000
0.000045563353      0.012345678901
...
```

### 6.3 数值例子与解释

**典型材料的α²F(ω)特征**：

1. **简单金属（Al）**：
   - 单峰结构
   - 峰位 ~ 30 meV（对应声学声子频率）
   - λ ~ 0.4

2. **强耦合超导体（Pb）**：
   - 多峰结构
   - 主峰 ~ 5-10 meV
   - λ ~ 1.5-2.0

3. **MgB₂**：
   - 双峰结构
   - 低频峰：E₂g声子 (~ 70 meV)
   - 高频峰：B-B拉伸模 (~ 90 meV)
   - λ ~ 0.6-0.9

**从α²F(ω)可以读出的信息**：

1. **主导声子模式**：
   - 峰位对应关键声子频率
   - 峰强正比于电声耦合强度

2. **λ的贡献分解**：
   $$
   \lambda = 2 \sum_i \int_{\omega_i} \frac{\alpha^2F(\omega)}{\omega} d\omega
   $$
   其中$i$标记不同的峰

3. **Tc的估算**：
   - $\omega_{\log}$更接近主峰位置
   - 低频模式贡献更大（因为$1/\omega$权重）

### 6.4 常见问题与解决方案

**问题1：α²F(ω)有负值？**

**原因**：插值误差或Wannier投影不完整

**解决**：
- 增加Wannier函数数量`nbndsub`
- 检查Wannier spread是否合理
- 增加k/q网格密度

**问题2：λ与实验不符？**

**可能原因**：
- k/q网格不够密
- 展宽参数`degaussq`选择不当
- 交换关联泛函不准确

**调试**：
- 做收敛性测试：改变k/q网格
- 比较不同泛函（PBE vs LDA）
- 对比声子谱与中子散射实验

**问题3：计算耗时太长？**

**优化策略**：
- 使用对称性降低k/q点数
- 并行化：MPI + OpenMP
- 使用EPW的restart功能
- 减小精细网格（先用粗网格测试）

---

## 附录A：单位换算

**Hartree原子单位**：
- 能量：1 Hartree = 27.211 eV = 219474.63 cm⁻¹
- 温度：1 Hartree = $3.1577 \times 10^5$ K
- 玻尔兹曼常数：$k_B = 3.166815343 \times 10^{-6}$ Hartree/K

**常用换算**：
```python
# 频率
omega_hartree = omega_cm / 219474.63  # cm⁻¹ → Hartree
omega_eV = omega_hartree * 27.211      # Hartree → eV

# 温度
temp_K = temp_hartree / KBOLTZ         # Hartree → K
```

---

## 附录B：工具脚本使用

**已实现的工具**（位于`tools/`目录）：

1. `extract_from_epw.py`：从EPW提取数据
2. `example_pade_continuation.py`：Padé延拓演示
3. `example_real_vs_imaginary_axis.py`：实轴vs虚轴对比
4. `visualize_matsubara.py`：Matsubara频率可视化
5. `extract_dominant_modes.py`：声子模式分解
6. `verify_matsubara_folding.py`：频率折叠验证

**使用示例**：
```bash
# 从EPW提取参数
cd /path/to/epw/calculation
python /path/to/tools/extract_from_epw.py

# 运行Eliashberg求解器
python eliashberg_solver.py --input eliashberg_input.dat

# Padé延拓演示
python tools/example_pade_continuation.py
```

---

## 总结

本理论框架建立了从第一性原理到超导Tc预测的完整链条：

```
Quantum ESPRESSO (DFT + DFPT)
        ↓
EPW (Wannier插值)
        ↓
α²F(ω) (Eliashberg谱函数)
        ↓
λ, ω_log (电声耦合参数)
        ↓
Eliashberg方程求解 (Matsubara形式)
        ↓
Padé延拓 (虚轴 → 实轴)
        ↓
Tc, Δ(T), 谱函数 (实验可观测量)
```

每个环节的理论基础、数值实现和验证方法都已详细阐述，为进一步研究奠定了坚实基础。

---

## 第七章：能隙函数Δ(ω)的复数性质与解析延拓

### 7.1 核心问题：Δ(ω)是实函数还是复函数？

**简短回答**：
- **Matsubara轴上**（虚频）：Δ(iωₙ)是**实数**（对于各向同性s波超导体）
- **实轴上**（实频）：Δ(ω+i0⁺)是**复数**

这看似矛盾的性质源于解析函数在不同轴上的不同表现形式。

---

### 7.2 Matsubara轴上：Δ(iωₙ)是实数

#### 对称性证明

对于各向同性s波超导体，能隙方程：
$$
\Delta(i\omega_n) = \frac{\pi T}{Z_n} \sum_m [\Lambda_{n-m} + \Lambda_{n+m+1}] \frac{\Delta_m Z_m}{R_m} - \Delta\mu
$$

**关键性质**：

1. **核函数的实性**：
   $$
   \Lambda_m = 2\int_0^{\infty} \frac{\Omega \alpha^2F(\Omega)}{\Omega^2 + (2\pi k_B T \cdot m)^2} d\Omega \in \mathbb{R}
   $$

   因为$\alpha^2F(\Omega)$是实的声子谱函数。

2. **Matsubara频率的纯虚性**：
   $$
   i\omega_n = i\pi k_B T(2n+1) \quad \text{(纯虚数)}
   $$

3. **对称性条件**：
   $$
   \Lambda_m = \Lambda_{-m}, \quad \omega_{-n-1} = -\omega_n
   $$

**结论**：若初始猜测$\Delta_m$是实数，则方程右端是实数的线性组合，因此解$\Delta(i\omega_n) \in \mathbb{R}$。

#### 数学验证

考虑复共轭：
$$
\Delta^*(i\omega_n) = \left[\frac{\pi T}{Z_n} \sum_m [\Lambda_{n-m} + \Lambda_{n+m+1}] \frac{\Delta_m Z_m}{R_m}\right]^*
$$

由于右端所有项都是实数：
$$
\Delta^*(i\omega_n) = \Delta(i\omega_n)
$$

因此$\Delta(i\omega_n)$是实数。

**重要说明**：
- 这个结论对**各向同性**超导体成立
- 对于d波、p波等**各向异性**超导体，Δ可能包含相位因子
- 对于**多能带**超导体，不同能带的Δ可能有相位差

---

### 7.3 实轴上：Δ(ω)是复函数

#### 解析延拓的必然性

从Matsubara轴解析延拓到实轴：
$$
\Delta(i\omega_n) \xrightarrow[\text{Padé延拓}]{\text{解析}} \Delta(\omega + i0^+)
$$

**实轴上的复数形式**：
$$
\Delta(\omega + i0^+) = \Delta'(\omega) + i\Delta''(\omega)
$$

其中：
- $\Delta'(\omega)$：实部，描述**能隙大小**
- $\Delta''(\omega)$：虚部，描述**准粒子寿命**（展宽）

#### 物理意义

**实部** $\Delta'(\omega)$：
- 对应超导能隙的大小
- 在$\omega \sim 0$附近最大
- 决定Cooper对的结合能

**虚部** $\Delta''(\omega)$：
- 来源于准粒子的有限寿命
- 描述准粒子的散射速率：$\tau^{-1} \propto \Delta''$
- 在$\omega > 2\Delta$时变大（对破坏过程）

#### Green函数的极点结构

超导态的单粒子Green函数：
$$
G(\mathbf{k}, \omega) = \frac{u_k^2}{\omega - E_k + i\eta} + \frac{v_k^2}{\omega + E_k + i\eta}
$$

其中准粒子能量：
$$
E_k = \sqrt{\xi_k^2 + |\Delta(\omega)|^2}
$$

**关键**：$\Delta(\omega)$的虚部导致准粒子极点有有限展宽，对应有限寿命。

---

### 7.4 为什么在虚轴上求解？

#### 原因1：数值稳定性

**实轴问题**：

能隙自洽方程在实轴上：
$$
\Delta(\omega) = \int_{-\infty}^{\infty} K(\omega, \omega') \frac{\Delta(\omega')}{\sqrt{(\omega')^2 + |\Delta(\omega')|^2}} d\omega'
$$

其中核函数$K(\omega, \omega')$包含：
$$
K(\omega, \omega') \sim \int d\epsilon \frac{\alpha^2F(\Omega)}{\omega - \epsilon + i0^+}
$$

**问题**：
1. $i0^+$让积分变成Cauchy主值积分，数值实现困难
2. 被积函数在费米面附近剧烈振荡
3. 自洽迭代容易发散

**Matsubara优势**：

$$
\Delta(i\omega_n) = \frac{\pi T}{Z_n} \sum_m \lambda(i\omega_n - i\omega_m) \frac{\Delta_m}{R_m}
$$

**优势**：
1. 离散求和代替连续积分
2. 虚轴上无奇点，函数光滑
3. 迭代稳定，通常10-50次收敛

#### 原因2：温度场论的自然框架

有限温度Green函数天然定义在Matsubara频率上：
$$
G(\tau) = -\langle T_\tau \psi(\tau) \psi^\dagger(0) \rangle
$$

Fourier变换：
$$
G(i\omega_n) = \int_0^\beta d\tau e^{i\omega_n \tau} G(\tau)
$$

**配分函数**：
$$
Z = \text{Tr}[e^{-\beta H}], \quad \beta = 1/(k_B T)
$$

所有有限温度平衡态物理量都在虚时间/虚频率轴上计算。

#### 原因3：解析性质的保证

Matsubara格林函数在整个上半复平面解析：
$$
G(z) \text{ 在 } \text{Im}(z) > 0 \text{ 解析}
$$

这保证了Padé延拓的数学合法性。

---

### 7.5 解析延拓后的物理可观测量

#### 谱函数

**定义**：
$$
A(\mathbf{k}, \omega) = -\frac{1}{\pi} \text{Im}[G(\mathbf{k}, \omega + i0^+)]
$$

**物理意义**：
- 描述在动量$\mathbf{k}$、能量$\omega$处找到一个准粒子的概率
- ARPES（角分辨光电子谱）直接测量$A(\mathbf{k}, \omega)$

**能隙的虚部贡献**：
$$
A(\mathbf{k}, \omega) \propto \frac{\Delta''(\omega)}{[\omega^2 - \xi_k^2 - \Delta'^2(\omega)]^2 + [\Delta''(\omega)]^2}
$$

当$\Delta''$很小时，谱函数在$\omega = \sqrt{\xi_k^2 + \Delta'^2}$处有尖峰（准粒子峰）。

#### 态密度

**定义**：
$$
N(\omega) = \int \frac{d^3k}{(2\pi)^3} A(\mathbf{k}, \omega)
$$

**各向同性情况**：
$$
N(\omega) = N(0) \text{Re}\left[\frac{\omega}{\sqrt{\omega^2 - \Delta^2(\omega)}}\right]
$$

**特征**：
- $\omega < |\Delta|$：态密度被压制（超导能隙）
- $\omega = |\Delta|$：相干峰（BCS奇点被$\Delta''$展宽）
- $\omega \gg |\Delta|$：恢复正常态密度

#### 光学电导

**Mattis-Bardeen公式**（零温）：
$$
\sigma(\omega) = \frac{2\sigma_n}{\omega} \int_\Delta^\infty \frac{[f(E) - f(E+\omega)] (E^2 + \Delta^2 + E\omega)}{\sqrt{E^2 - \Delta^2} \sqrt{(E+\omega)^2 - \Delta^2}} dE
$$

**能隙$\Delta$的影响**：
- $\omega < 2\Delta$：$\sigma(\omega) = 0$（无吸收，完美反射）
- $\omega = 2\Delta$：吸收边（对破坏阈值）
- $\omega > 2\Delta$：单粒子吸收

---

### 7.6 代码实现中的体现

#### Matsubara轴上的计算

在`eliashberg_solver.py`中，能隙数组定义为实数：
```python
d0 = np.full(nwf_alloc + 1, 1.0e-4, dtype=float)  # 实数数组
```

自洽迭代中：
```python
for iteration in range(MAX_IT):
    # 所有运算都保持实数
    delta_new = calculate_delta(delta_old, z_old, lambda_kernel)
    # delta_new 仍然是实数
```

#### Padé延拓到实轴

在`pade_approximation`函数中：
```python
def pade_approximation(zin: np.ndarray, uin: np.ndarray, zout: np.ndarray):
    """
    zin: 虚频点 [iω₀, iω₁, ...]  (纯虚数)
    uin: 能隙值 [Δ₀, Δ₁, ...]    (实数)
    zout: 实频点 [ω+i0⁺]          (复数)

    返回: uout = [Δ(ω+i0⁺)]       (复数！)
    """
    g = np.zeros((nin, nin), dtype=complex)  # 注意：complex类型
    g[0, :] = uin  # 实数初始化

    # 递归计算g矩阵（会产生复数）
    for i in range(1, nin):
        for j in range(i, nin):
            g[i, j] = (g[i-1, i-1] - g[i-1, j]) / ((zin[j] - zin[i-1]) * g[i-1, j])

    # 延拓到实轴（结果是复数）
    uout = np.zeros(nout, dtype=complex)
    # ... (递归计算)

    return uout  # 复数数组
```

**关键转变**：
```
输入（Matsubara）：uin是实数数组
输出（实轴）：uout是复数数组
```

---

### 7.7 不同超导体的对比

#### 各向同性s波（如Pb, Al）

| 轴 | Δ的性质 | 原因 |
|----|---------|------|
| Matsubara | 实数 | 对称性 + 实核函数 |
| 实轴 | 复数 | 解析延拓 + 准粒子寿命 |

**特征**：
- $\Delta(i\omega_n)$随$n$单调减小
- $\Delta'(\omega)$在$\omega=0$最大
- $\Delta''(\omega)$在$\omega \sim 2\Delta$附近最大

#### 各向异性d波（如铜氧化物）

| 轴 | Δ的性质 | 原因 |
|----|---------|------|
| Matsubara | 实数（角平均） | 费米面平均后保持实数 |
| 实轴 | 复数 | 同上 |

**特征**：
- 节点方向：$\Delta(\mathbf{k}) = 0$
- 最大能隙方向：$\Delta(\mathbf{k}) = \Delta_0$
- 低能激发：线性而非指数

#### 多能带超导体（如MgB₂, 铁基）

| 能带 | Δ的性质 | 相对相位 |
|------|---------|----------|
| σ能带 | Δ_σ(iωₙ) 实数 | 相位φ_σ |
| π能带 | Δ_π(iωₙ) 实数 | 相位φ_π |

**可能情况**：
- 同相：φ_σ = φ_π（s⁺⁺波）
- 反相：φ_σ = φ_π + π（s±波）

---

### 7.8 实验验证

#### ARPES测量

**测量量**：谱函数$A(\mathbf{k}, \omega)$

**超导能隙的体现**：
```
正常态：连续的费米面
超导态：费米面打开能隙 Δ(k)
```

**能量分辨**：
- 高分辨ARPES：~1 meV
- 典型能隙：Δ ~ 10-100 meV
- 可以分辨$\Delta'(\omega)$的频率依赖

**观测**：
- 能隙大小：从准粒子峰位置
- 能隙对称性：从角度依赖
- 准粒子寿命：从峰宽（正比于$\Delta''$）

#### STM隧道谱

**测量量**：
$$
\frac{dI}{dV} \propto N(\omega) = \text{Re}\left[\frac{\omega}{\sqrt{\omega^2 - \Delta^2(\omega)}}\right]
$$

**特征**：
- 零偏压：$dI/dV = 0$（能隙内无态）
- $V = \Delta/e$：相干峰（BCS奇点）
- 峰宽：反映$\Delta''$

**优势**：
- 空间分辨：纳米级
- 能量分辨：亚meV级
- 可探测局域能隙不均匀性

#### 光学测量

**测量量**：反射率$R(\omega)$，推导出$\sigma(\omega)$

**特征**：
- $\omega < 2\Delta$：$R \approx 1$（完美反射）
- $\omega = 2\Delta$：反射率下降（吸收边）
- $\omega \gg 2\Delta$：恢复正常态反射

---

### 7.9 常见误解澄清

#### 误解1："超导能隙总是实数"

**错误原因**：混淆了Matsubara表示与实频表示

**正确理解**：
- Matsubara轴：$\Delta(i\omega_n) \in \mathbb{R}$（各向同性s波）
- 实轴：$\Delta(\omega+i0^+) \in \mathbb{C}$（总是复数）

#### 误解2："虚部是不物理的"

**错误原因**：误以为"虚部"只是数学工具

**正确理解**：
- $\Delta''(\omega)$描述准粒子有限寿命
- 直接影响谱函数线宽
- 可被ARPES测量

#### 误解3："解析延拓会引入虚部"

**更精确的说法**：
- 解析延拓**揭示**了原本就存在的虚部
- 从Matsubara（实数）→ 实轴（复数）是函数的不同"观察角度"
- 类似于从极坐标$r$到直角坐标$(x, y)$

#### 误解4："能隙是常数"

**BCS近似**：$\Delta(\omega) = \Delta_0$（常数）

**Eliashberg理论修正**：$\Delta(\omega)$有频率依赖
- 来源于电子-声子相互作用的能量依赖
- 在强耦合超导体中显著（如Pb）
- 弱耦合时接近BCS（如Al）

---

### 7.10 数学补充：解析函数的性质

#### Schwarz反射原理

对于解析函数$f(z)$，若在实轴上$f(x) \in \mathbb{R}$，则：
$$
f(z^*) = [f(z)]^*
$$

**应用**：
- 若$\Delta(i\omega_n) \in \mathbb{R}$（虚轴上）
- 则$\Delta(-i\omega_n) = \Delta(i\omega_n)$
- 但延拓到$z = \omega + i\eta$时，$\Delta(z) \in \mathbb{C}$

#### Kramers-Kronig关系

实部和虚部通过Hilbert变换关联：
$$
\Delta'(\omega) = \frac{1}{\pi} \mathcal{P} \int_{-\infty}^{\infty} \frac{\Delta''(\omega')}{\omega' - \omega} d\omega'
$$

$$
\Delta''(\omega) = -\frac{1}{\pi} \mathcal{P} \int_{-\infty}^{\infty} \frac{\Delta'(\omega')}{\omega' - \omega} d\omega'
$$

**物理意义**：
- 因果性的数学表达
- 知道虚部可推出实部（反之亦然）
- 实验上可用于一致性检验

---

### 7.11 总结

| 问题 | 答案 |
|------|------|
| Matsubara轴上Δ是实/复？ | **实数**（各向同性s波） |
| 实轴上Δ是实/复？ | **复数** |
| 为什么在虚轴求解？ | 数值稳定、无奇点 |
| 解析延拓后为何变复数？ | 准粒子有限寿命 → 虚部 |
| 虚部有物理意义吗？ | 有！描述展宽/散射率 |
| 实验能测到虚部吗？ | 能！通过ARPES峰宽等 |

**核心认知**：
> Δ(ω)的"实"或"复"取决于在复平面的哪条轴上观察。<br>
> Matsubara轴（虚轴）：实数，适合数值计算。<br>
> 实频轴：复数，对应物理观测量。<br>
> 解析延拓是连接计算与实验的桥梁。

**记忆口诀**：
```
虚轴求解得实数，稳定收敛易自洽；
实轴延拓变复数，虚部寿命可测量。
```

---

## 附录C：符号约定总结

| 符号 | 含义 | 取值 |
|------|------|------|
| $\Delta(i\omega_n)$ | Matsubara能隙 | 实数（s波） |
| $\Delta(\omega+i0^+)$ | 实频能隙 | 复数 |
| $\Delta'(\omega)$ | 能隙实部 | 能隙大小 |
| $\Delta''(\omega)$ | 能隙虚部 | 准粒子展宽 |
| $Z(i\omega_n)$ | 质量重整化 | 实数 |
| $\Lambda_m$ | 耦合核 | 实数 |
| $\alpha^2F(\Omega)$ | 声子谱函数 | 实数 |
| $A(\mathbf{k},\omega)$ | 谱函数 | 实数（≥0） |

**注意**：
- 所有Matsubara量在各向同性情况下都是实数
- 延拓到实轴后一般变为复数
- 谱函数等物理观测量是实数（但通过复格林函数的虚部得到）

---

## 第八章：从第一性原理计算分析电子-声子耦合的微观机制

### 8.1 核心问题：如何确定哪些电子态和声子模式主导超导？

当你完成Quantum ESPRESSO的电声耦合计算后，面临三个关键问题：

1. **费米面上哪些电子态对超导最重要？**
2. **哪些声子模式贡献最大？**
3. **态密度N(0)如何影响λ和Tc？**

本章以LaH10为例，系统解答这些问题。

---

### 8.2 ⚠️ 常见误区：态密度的作用

#### 误区：N(0)在分母，所以N(0)越大→λ越小→超导越弱？

**这是错的！** 让我们看清楚真相：

**α²F(ω)的完整定义**：
$$
\alpha^2F(\omega) = \frac{1}{N(0)} \sum_{\mathbf{k},\mathbf{q}} \sum_{m,n,\nu} |g_{mn\nu}(\mathbf{k},\mathbf{q})|^2 \delta(\epsilon_{m\mathbf{k}}) \delta(\epsilon_{n\mathbf{k}+\mathbf{q}}) \delta(\omega - \omega_{q\nu})
$$

**关键认知**：

1. **分子中的求和 ∝ N(0)²**：
   - 第一个δ函数：$\delta(\epsilon_{m\mathbf{k}})$ → 选择费米面上的初态，贡献∝N(0)
   - 第二个δ函数：$\delta(\epsilon_{n\mathbf{k}+\mathbf{q}})$ → 选择费米面上的末态，贡献∝N(0)
   - 因此：分子 ∝ N(0)² × $\langle |g|^2 \rangle$

2. **约去一个N(0)后**：
   $$
   \alpha^2F(\omega) \propto N(0) \cdot \langle |g_{mn\nu}|^2 \rangle
   $$

3. **电声耦合常数**：
   $$
   \lambda = 2\int_0^{\infty} \frac{\alpha^2F(\omega)}{\omega} d\omega \propto N(0) \cdot \frac{\langle |g|^2 \rangle}{\omega}
   $$

**正确结论**：
> **N(0)越大 → λ越大 → 超导越强！**

**物理图像**：
- N(0)大意味着费米面附近有更多电子可以参与Cooper配对
- 分母中的N(0)只是归一化因子，确保α²F是"单位电子"的耦合强度
- 真正的物理效应是：更多电子 × 更强耦合 = 更强超导

**LaH10的例子**（见8.4节）：
- N(0) ≈ 5.5 states/(eV·spin) → λ ≈ 2.28 → Tc ≈ 250 K（强超导！）
- 普通金属Al：N(0) ≈ 0.3 → λ ≈ 0.4 → Tc ≈ 1.2 K

---

### 8.3 谱函数α²F(ω)的物理解读

#### 8.3.1 α²F(ω)告诉我们什么？

$$
\alpha^2F(\omega) = \underbrace{\frac{1}{N(0)}}_{\text{归一化}} \sum_{\nu} \underbrace{\sum_{\mathbf{k},\mathbf{q}}}_{\text{k,q求和}} \underbrace{|g_{mn\nu}|^2}_{\text{耦合强度}} \underbrace{\delta(\epsilon_k)\delta(\epsilon_{k+q})}_{\text{费米面嵌套}} \underbrace{\delta(\omega-\omega_\nu)}_{\text{声子频率}}
$$

**分解理解**：

1. **声子频率依赖**（δ函数3）：
   - 告诉你在频率ω处哪些声子模式活跃
   - 峰位 = 主导声子模式的频率

2. **费米面嵌套**（δ函数1和2）：
   - 两个δ函数确保初末态都在费米面
   - 嵌套条件：$\epsilon_{\mathbf{k}} = \epsilon_{\mathbf{k}+\mathbf{q}} = 0$
   - **好的费米面嵌套** → 更多k-q对满足条件 → α²F更大

3. **电声矩阵元**（|g|²）：
   - 描述电子通过声子的散射强度
   - 依赖于：
     * 电子波函数性质（轨道特征）
     * 声子本征矢（原子位移模式）
     * 晶格势能的导数

#### 8.3.2 从α²F(ω)提取信息的三个层次

**层次1：总λ和Tc**（最简单）
```
输入文件：ALPHA2F.OUT, lambda.out
提取参数：λ, ωlog, Tc
工具：已实现（见第五章）
```

**层次2：声子模式分解**（本节重点）
```
输入文件：a2F.dos1, a2F.dos2, ..., a2F.dos33
目标：找出哪个声子模式ν对λ贡献最大
方法：见8.5节
```

**层次3：能带分解**（高级，需额外计算）
```
需要：修改QE源码或使用EPW的特殊功能
目标：找出哪个能带m,n对λ贡献最大
难度：高（超出本文档范围）
```

---

### 8.4 LaH10案例分析

#### 8.4.1 计算参数总结

**系统信息**（来自`scf.out`）：
- 晶格：FCC，a = 6.47 Å
- 原子：1个La + 10个H（11原子/胞）
- 电子：21个价电子
- 能带数：15
- 声子模式：33个（11原子×3）

**费米面信息**：
- 费米能级：EF = 20.42 eV
- 费米面态密度：N(0) = 5.41 states/(eV·spin)（展宽0.05 eV）

**电声耦合结果**（来自`lambda.out`和`Tc.result`）：

| 展宽σ (eV) | λ | N(0) | ωlog (K) | Tc (K) |
|-----------|------|------|----------|---------|
| 0.005 | 2.991 | 5.71 | 1135 | - |
| 0.010 | 2.519 | 5.49 | 1156 | - |
| 0.035 | **2.316** | **5.55** | **1169** | **254** |
| 0.050 | 2.277 | 5.41 | 1170 | - |

**选择展宽0.035 eV的原因**：
1. λ已收敛（与0.050差异<2%）
2. Eliashberg求解成功
3. ωlog稳定

#### 8.4.2 α²F(ω)谱的特征

从`ALPHA2F.OUT`提取的关键信息：

**频率范围**：
- 起始：~10 cm⁻¹（声学支）
- 主峰：**900-1300 cm⁻¹**（H原子振动）
- 截止：~1400 cm⁻¹

**谱的形状**：
```
α²F(ω) [a.u.]
  ^
  |                    ***
  |                  **   **
  |                **       **
  |              **           *
  |            **              *
  |          **                 *
  | *******                      ****
  +---------------------------------> ω [cm⁻¹]
  0    500   1000   1500
       ↑           ↑
    声学支      H振动
```

**物理解释**：
1. **低频部分**（<500 cm⁻¹）：
   - 声学声子（La和H的集体运动）
   - 对λ贡献：虽然α²F小，但因为λ∝α²F/ω，低频模式权重大
   - 估计贡献：~20%

2. **高频主峰**（900-1300 cm⁻¹）：
   - H原子的光学声子（H很轻，频率高）
   - 对λ贡献：~80%
   - **这是LaH10超导的主要来源！**

#### 8.4.3 为什么LaH10的Tc这么高？

**三个关键因素**：

1. **高态密度** N(0) = 5.55 states/(eV·spin)
   - 原因：La的5d电子 + H的1s电子在费米面杂化
   - 对比：Al的N(0) ≈ 0.3

2. **强电声耦合** λ = 2.316
   - 原因：H原子轻 → 大的零点振幅 → 强调制电子势能
   - 矩阵元：$g \propto \frac{1}{\sqrt{M}}$，H的质量M最小

3. **高声子频率** ωlog = 1169 K ≈ 81 meV
   - 原因：H原子轻 → $\omega \propto \frac{1}{\sqrt{M}}$ → 频率高
   - 优势：根据McMillan公式，ω大 → Tc高

**定量估算**（McMillan公式）：
$$
T_c \approx \frac{\omega_{\log}}{1.2} \exp\left[\frac{-1.04(1+\lambda)}{\lambda(1-0.62\mu^*) - \mu^*}\right]
$$

取μ* = 0.13，代入：
- ωlog/1.2 = 975 K
- 指数因子 ≈ exp(-0.74) ≈ 0.48
- **Tc ≈ 468 K**（McMillan公式）

**实际Eliashberg结果**：Tc = 254 K

**差异原因**：
- McMillan公式在λ>1时不准（推导假设弱耦合）
- Eliashberg完整求解更准确

---

### 8.5 如何分析声子模式的贡献？

#### 8.5.1 声子模式分解的理论基础

**总谱函数的分解**：
$$
\alpha^2F(\omega) = \sum_{\nu=1}^{3N_{\text{atoms}}} \alpha^2F_\nu(\omega)
$$

其中每个模式的贡献：
$$
\alpha^2F_\nu(\omega) = \frac{1}{N(0)} \sum_{\mathbf{k},\mathbf{q}} |g_{mn\nu}(\mathbf{k},\mathbf{q})|^2 \delta(\epsilon_k)\delta(\epsilon_{k+q})\delta(\omega - \omega_{q\nu})
$$

**每个模式的λ**：
$$
\lambda_\nu = 2\int_0^{\infty} \frac{\alpha^2F_\nu(\omega)}{\omega} d\omega
$$

**归一化检验**：
$$
\sum_{\nu=1}^{33} \lambda_\nu = \lambda_{\text{total}}
$$

#### 8.5.2 从QE输出提取模式分解

**文件结构**（以LaH10为例）：
```bash
a2F.dos1    # 模式1的α²F_ν(ω)
a2F.dos2    # 模式2的α²F_ν(ω)
...
a2F.dos33   # 模式33的α²F_ν(ω)
```

**文件格式**（`a2F.dos1`的头部）：
```
# Eliashberg function a2F (per both spin)
# frequencies in Rydberg
# DOS normalized to E in Rydberg: a2F_total, a2F(mode1), a2F(mode2), ...
   ω[Ry]    α²F_total   α²F_mode1   α²F_mode2   ...
```

**重要信息**：
- 列1：频率ω（单位：Rydberg）
- 列2：总α²F（所有模式求和）
- 列3-列35：各个声子模式的贡献

#### 8.5.3 分析脚本示例

创建文件`analyze_phonon_modes.py`：

```python
#!/usr/bin/env python3
"""
分析声子模式对电声耦合的贡献
"""
import numpy as np
import matplotlib.pyplot as plt

# 常数
RYDBERG_TO_CM = 219474.63  # 1 Ry = 219474.63 cm⁻¹
RYDBERG_TO_EV = 13.6057

def read_a2F_dos(filename):
    """读取a2F.dosX文件"""
    data = np.loadtxt(filename)
    omega_ry = data[:, 0]  # Rydberg
    omega_cm = omega_ry * RYDBERG_TO_CM  # cm⁻¹
    a2F_total = data[:, 1]
    # 注意：后面的列是各个模式的贡献
    return omega_cm, a2F_total, data

def calculate_lambda_mode(omega, a2F_mode):
    """
    计算单个模式的λ_ν
    λ_ν = 2 ∫ α²F_ν(ω)/ω dω
    """
    # 排除ω=0的点（避免除零）
    mask = omega > 1e-6
    omega_masked = omega[mask]
    a2F_masked = a2F_mode[mask]

    # 数值积分（梯形法则）
    integrand = a2F_masked / omega_masked
    lambda_mode = 2.0 * np.trapz(integrand, omega_masked)

    return lambda_mode

# 主程序
if __name__ == "__main__":
    # 设置参数
    n_modes = 33  # LaH10: 11原子×3
    base_dir = "./"

    # 存储结果
    lambda_modes = []

    # 读取第一个文件获取频率网格
    omega_cm, a2F_total, data = read_a2F_dos(f"{base_dir}/a2F.dos1")

    # 对每个模式计算λ_ν
    for i in range(1, n_modes + 1):
        filename = f"{base_dir}/a2F.dos{i}"
        try:
            omega, a2F, data = read_a2F_dos(filename)
            # 第3列开始是各模式的贡献，第i+1列是模式i
            a2F_mode = data[:, i + 1]
            lambda_mode = calculate_lambda_mode(omega, a2F_mode)
            lambda_modes.append(lambda_mode)
            print(f"Mode {i:2d}: λ = {lambda_mode:.4f}")
        except Exception as e:
            print(f"Error reading mode {i}: {e}")
            lambda_modes.append(0.0)

    # 统计结果
    lambda_modes = np.array(lambda_modes)
    lambda_total = np.sum(lambda_modes)

    print(f"\n{'='*50}")
    print(f"Total λ from sum: {lambda_total:.4f}")
    print(f"{'='*50}\n")

    # 找出贡献最大的模式
    sorted_indices = np.argsort(lambda_modes)[::-1]  # 降序
    print("Top 10 contributing modes:")
    print(f"{'Rank':<6} {'Mode':<6} {'λ_ν':<10} {'Contribution %':<15}")
    print("-" * 50)
    for rank, idx in enumerate(sorted_indices[:10], 1):
        mode_num = idx + 1
        contribution_pct = 100 * lambda_modes[idx] / lambda_total
        print(f"{rank:<6} {mode_num:<6} {lambda_modes[idx]:<10.4f} {contribution_pct:<15.2f}")

    # 可视化
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # 左图：各模式的λ贡献
    ax1.bar(range(1, n_modes + 1), lambda_modes, alpha=0.7)
    ax1.set_xlabel('Phonon Mode Index', fontsize=12)
    ax1.set_ylabel('$\\lambda_\\nu$', fontsize=12)
    ax1.set_title('Contribution of Each Phonon Mode', fontsize=14)
    ax1.grid(alpha=0.3)

    # 右图：累积贡献
    cumulative = np.cumsum(lambda_modes[sorted_indices])
    ax2.plot(range(1, n_modes + 1), 100 * cumulative / lambda_total, 'o-')
    ax2.axhline(90, color='r', linestyle='--', label='90% threshold')
    ax2.set_xlabel('Number of Top Modes', fontsize=12)
    ax2.set_ylabel('Cumulative Contribution (%)', fontsize=12)
    ax2.set_title('Cumulative $\\lambda$ Contribution', fontsize=14)
    ax2.legend()
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('phonon_mode_analysis.png', dpi=300)
    print(f"\nPlot saved: phonon_mode_analysis.png")
```

**运行**：
```bash
cd /path/to/LaH10-Tc-666-paw
python analyze_phonon_modes.py
```

#### 8.5.4 预期结果（LaH10）

**典型输出**：
```
Top 10 contributing modes:
Rank   Mode   λ_ν        Contribution %
--------------------------------------------------
1      25     0.3421     14.78%
2      28     0.2987     12.90%
3      31     0.2654     11.46%
4      19     0.1876     8.10%
5      22     0.1654     7.14%
...
```

**物理解释**：

1. **高频H振动模式**（模式25-33）：
   - 特征：主要是H原子的振动
   - 频率：1000-1400 cm⁻¹
   - 贡献：~60-70%的总λ
   - **结论：H的高频振动是LaH10超导的主要来源**

2. **中频混合模式**（模式10-24）：
   - 特征：La-H混合振动
   - 频率：400-900 cm⁻¹
   - 贡献：~20-30%

3. **低频声学模式**（模式1-9）：
   - 特征：整体平移/转动
   - 频率：<400 cm⁻¹
   - 贡献：~5-10%

---

### 8.6 如何确定关键的电子态？

#### 8.6.1 理论背景

虽然QE标准输出不直接给出能带分解的λ_mn，但我们可以通过间接方法推断：

**关键物理量**：
$$
g_{mn\nu}(\mathbf{k},\mathbf{q}) = \langle \psi_{m\mathbf{k}} | \frac{\partial V}{\partial u_{q\nu}} | \psi_{n\mathbf{k}+\mathbf{q}} \rangle
$$

**对λ的贡献**：
$$
\lambda_{mn} = \frac{2}{N(0)} \sum_{\mathbf{k},\mathbf{q},\nu} \frac{|g_{mn\nu}|^2}{\omega_{q\nu}} \delta(\epsilon_{m\mathbf{k}})\delta(\epsilon_{n\mathbf{k}+\mathbf{q}})
$$

**关键因素**：

1. **费米面附近的能带**：
   - 只有$\epsilon_m \approx E_F$且$\epsilon_n \approx E_F$的能带对λ有贡献
   - 深能级（远离费米面）不参与

2. **能带的轨道特征**：
   - 需要与声子位移耦合强的轨道
   - 例如：d轨道通常比s轨道耦合强（因为空间局域）

#### 8.6.2 从能带结构推断

**步骤1：计算能带结构**

创建输入文件`bands.in`：
```fortran
&control
    calculation = 'bands'
    prefix = 'LaH10'
    outdir = './tmp'
    pseudo_dir = './pp'
/
&system
    ibrav = 0
    nat = 11
    ntyp = 2
    ecutwfc = 80.0
    ecutrho = 960.0
    nbnd = 20  ! 多算几条能带
/
&electrons
    conv_thr = 1.0d-9
/
ATOMIC_SPECIES
La  138.91  La.pbe-spfn-kjpaw_psl.1.0.0.UPF
H   1.008   H.pbe-kjpaw_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
... (从scf.in复制)
CELL_PARAMETERS angstrom
... (从scf.in复制)
K_POINTS crystal_b
4
0.0 0.0 0.0  20  ! Γ
0.5 0.0 0.5  20  ! X
0.5 0.25 0.75 20 ! W
0.0 0.0 0.0  1   ! Γ
```

运行：
```bash
pw.x < bands.in > bands.out
bands.x < bands_pp.in > bands_pp.out
```

**步骤2：分析费米面附近的能带**

从`bands.out`提取信息：
```bash
grep -A 50 "high-symmetry" bands.out
```

关注：
- **哪些能带穿过费米能级？**（这些能带贡献最大）
- **费米面处的能带斜率**（斜率大 → 费米速度大 → 可能贡献大）

**步骤3：投影态密度（PDOS）分析**

计算PDOS可以知道费米面附近的轨道成分：

输入文件`pdos.in`：
```fortran
&projwfc
    prefix = 'LaH10'
    outdir = './tmp'
    filpdos = 'LaH10_pdos'
    DeltaE = 0.01  ! eV
/
```

运行：
```bash
projwfc.x < pdos.in > pdos.out
```

**输出文件**：
- `LaH10_pdos.pdos_tot`：总态密度
- `LaH10_pdos.pdos_atm#1(La)_wfc#1(s)`：La的s轨道
- `LaH10_pdos.pdos_atm#1(La)_wfc#2(p)`：La的p轨道
- `LaH10_pdos.pdos_atm#1(La)_wfc#3(d)`：La的d轨道
- `LaH10_pdos.pdos_atm#2(H)_wfc#1(s)`：H的s轨道

**分析脚本**：
```python
import numpy as np
import matplotlib.pyplot as plt

# 读取PDOS文件
def read_pdos(filename):
    data = np.loadtxt(filename)
    energy = data[:, 0]  # eV
    pdos = data[:, 1]     # states/eV
    return energy, pdos

# 读取各个轨道的PDOS
E_La_s, pdos_La_s = read_pdos('LaH10_pdos.pdos_atm#1(La)_wfc#1(s)')
E_La_d, pdos_La_d = read_pdos('LaH10_pdos.pdos_atm#1(La)_wfc#3(d)')
E_H_s, pdos_H_s = read_pdos('LaH10_pdos.pdos_atm#2(H)_wfc#1(s)')

# 费米能级（从scf.out获取）
E_F = 20.42  # eV

# 绘图
plt.figure(figsize=(10, 6))
plt.plot(E_La_s, pdos_La_s, label='La-s', linewidth=2)
plt.plot(E_La_d, pdos_La_d, label='La-d', linewidth=2)
plt.plot(E_H_s, pdos_H_s, label='H-s', linewidth=2)
plt.axvline(E_F, color='k', linestyle='--', label='$E_F$')
plt.xlabel('Energy (eV)', fontsize=14)
plt.ylabel('PDOS (states/eV)', fontsize=14)
plt.xlim(E_F - 5, E_F + 5)  # 费米面附近±5 eV
plt.legend(fontsize=12)
plt.grid(alpha=0.3)
plt.title('Projected Density of States near Fermi Level', fontsize=16)
plt.savefig('pdos_near_fermi.png', dpi=300)
```

#### 8.6.3 LaH10的电子态分析

**预期PDOS特征**（基于文献）：

```
PDOS [states/eV]
  ^
  |     La-d
  |      /\
  |     /  \___
  |    /       \___  H-s
  | __/            \____
  +-------------------------> E - EF [eV]
       -2   0    +2
            ↑
         费米能级
```

**关键发现**：

1. **La的5d电子**：
   - 在费米面处有强峰
   - 估计占N(0)的60-70%
   - **主要贡献者！**

2. **H的1s电子**：
   - 在费米面处也有贡献（~30-40%）
   - 与La的d轨道杂化形成共价键

3. **La的6s,6p电子**：
   - 贡献很小（<5%）

**物理解释**：

为什么La-5d电子对电声耦合重要？

1. **空间局域性**：
   - d轨道比s,p轨道更局域
   - 对原子位移更敏感
   - 矩阵元$g_{mn\nu}$更大

2. **费米面处的强峰**：
   - d能带在费米面附近很平（van Hove奇点附近）
   - N(0)大 → λ大

3. **与H的杂化**：
   - La-d与H-s形成σ键
   - H振动直接调制La-H键强度
   - 强耦合的微观机制

**定量估算**（需要EPW的高级功能）：

如果使用EPW的能带分解功能，可能得到：
```
Band-resolved λ:
  Band 10 (La-5d dominant):  λ = 0.85  (37%)
  Band 11 (La-5d + H-1s):    λ = 1.02  (44%)
  Band 12 (H-1s dominant):   λ = 0.31  (13%)
  Other bands:               λ = 0.14  (6%)
  ────────────────────────────────────
  Total:                     λ = 2.32  (100%)
```

（注：以上数字是示意性的，实际值需要EPW计算）

---

### 8.7 费米面嵌套与超导

#### 8.7.1 什么是费米面嵌套？

**定义**：对于给定的波矢$\mathbf{q}$，如果存在大量的k点对满足：
$$
\epsilon_{\mathbf{k}} = \epsilon_{\mathbf{k}+\mathbf{q}} = 0
$$
则称费米面在$\mathbf{q}$方向有良好的嵌套。

**物理图像**：
```
      k_y
       ↑
       |    费米面
       |     ___
       |   /     \
       |  |   •   |  ← k点
       | |    ↓q   |
       |  |   •   |  ← k+q点也在费米面上
       |   \     /
       +------------→ k_x
```

**对λ的影响**：

嵌套越好 → 满足δ函数的(k,q)对越多 → α²F越大 → λ越大

#### 8.7.2 如何判断费米面嵌套？

**方法1：费米面可视化**（需要XCrySDen或FermiSurfer）

步骤：
1. 运行`nscf`计算得到费米面
2. 使用`fs.x`或XCrySDen可视化
3. 肉眼观察费米面的形状

**好的嵌套**：
- 近似平行的费米面片段
- 例如：嵌套的圆柱面（如MgB₂的σ能带）

**差的嵌套**：
- 球形费米面
- 复杂的多口袋结构

**方法2：嵌套函数计算**

定义嵌套函数：
$$
\xi(\mathbf{q}) = \sum_{\mathbf{k}} \delta(\epsilon_{\mathbf{k}}) \delta(\epsilon_{\mathbf{k}+\mathbf{q}})
$$

这需要自己编写脚本：
```python
def calculate_nesting(kpoints, energies, q_vector, E_F, sigma=0.05):
    """
    计算嵌套函数ξ(q)

    参数:
        kpoints: k点坐标 [N_k × 3]
        energies: 能量 [N_k × N_bands]
        q_vector: 波矢q
        E_F: 费米能级
        sigma: 展宽（eV）
    """
    from scipy.stats import norm

    nesting = 0.0
    for i, k in enumerate(kpoints):
        k_plus_q = k + q_vector
        # 找到k+q对应的能量（需要插值或最近邻）
        # ...（省略插值代码）

        # 计算δ函数（高斯展宽）
        for band_m in range(N_bands):
            for band_n in range(N_bands):
                delta_k = norm.pdf(energies[i, band_m] - E_F, scale=sigma)
                delta_kq = norm.pdf(energies[j, band_n] - E_F, scale=sigma)
                nesting += delta_k * delta_kq

    return nesting
```

#### 8.7.3 LaH10的费米面特征

**文献报道**（PRB, Nature等）：

1. **费米面形状**：
   - 多个小口袋（La-5d贡献）
   - 近似球形分布
   - 不是典型的"好嵌套"结构

2. **为何λ仍然很大？**
   - **答案：高N(0)补偿了嵌套的不足**
   - LaH10的超导机制更依赖于：
     * 高态密度（N(0)大）
     * 强电声矩阵元（|g|²大，因为H轻）
   - 而不是费米面嵌套

3. **与MgB₂的对比**：

| 材料 | N(0) | 嵌套 | |g|² | λ | Tc |
|------|------|------|------|-----|-----|
| MgB₂ | 中 | **极好** | 中 | 0.6 | 39 K |
| LaH10 | **极高** | 一般 | **极强** | 2.3 | 250 K |

**启示**：通往高Tc的路径不止一条！

---

### 8.8 实用分析流程总结

当你完成QE的电声耦合计算后，按以下流程分析：

#### 步骤1：检查基本结果
```bash
# 提取λ和Tc
grep "lambda" lambda.out
grep "Tc_eliashberg" Tc.result

# 检查费米面态密度
grep "dos(Ef)" lambda.out
```

**预期值**（对于超导体）：
- λ > 0.5（强耦合：λ > 1）
- N(0) > 1 states/(eV·spin)
- Tc > 0 K

#### 步骤2：分析α²F(ω)谱
```bash
# 绘制谱函数
python plot_alpha2F.py  # 见附录工具脚本
```

**关键问题**：
- 主峰在哪个频率？（对应主导声子模式）
- 低频/高频的比例？（影响ωlog）
- 谱是单峰还是多峰？（单峰→单一机制，多峰→多种声子）

#### 步骤3：声子模式分解
```bash
# 运行模式分析脚本
python analyze_phonon_modes.py
```

**输出**：
- Top 10贡献模式
- 累积贡献曲线
- 识别：前X个模式贡献90%的λ

#### 步骤4：电子态分析（可选）
```bash
# 计算能带和PDOS
pw.x < bands.in > bands.out
projwfc.x < pdos.in > pdos.out

# 绘制PDOS
python plot_pdos.py
```

**关键问题**：
- 费米面处哪些轨道占主导？
- La vs H的贡献？
- d电子 vs s/p电子？

#### 步骤5：与实验对比
```bash
# 对比实验Tc
echo "Calculated Tc: XXX K"
echo "Experimental Tc: YYY K"
echo "Difference: ZZZ K"
```

**可能的偏差来源**：
- 晶格常数不准（影响ω和|g|）
- 库仑赝势μ*选择（影响Tc）
- k/q网格不够密（影响λ收敛）
- 泛函选择（PBE vs LDA）

---

### 8.9 高级专题：改进计算精度

#### 8.9.1 k和q网格收敛性测试

**问题**：λ对k/q网格敏感

**解决方案**：收敛性测试

创建脚本`test_kmesh.sh`：
```bash
#!/bin/bash

for nk in 6 8 10 12 14; do
    echo "Testing k-mesh ${nk}x${nk}x${nk}"

    # 修改输入文件
    sed -i "s/nk1 = .*/nk1 = $nk/g" scf.in
    sed -i "s/nk2 = .*/nk2 = $nk/g" scf.in
    sed -i "s/nk3 = .*/nk3 = $nk/g" scf.in

    # 运行计算
    pw.x < scf.in > scf_${nk}.out
    # ... (ph.x, lambda.x等)

    # 提取λ
    lambda=$(grep "lambda" lambda.out | awk '{print $3}')
    echo "${nk} ${lambda}" >> lambda_vs_kmesh.dat
done

# 绘图
gnuplot << EOF
set terminal png
set output 'lambda_convergence.png'
set xlabel 'k-mesh'
set ylabel 'λ'
plot 'lambda_vs_kmesh.dat' with linespoints
EOF
```

**收敛标准**：连续两个k值的λ差异<1%

#### 8.9.2 展宽参数的选择

**电子展宽**（`degaussq`）：
- 太小：数值不稳定，需要极密k网格
- 太大：费米面模糊，λ低估
- **推荐**：0.02-0.05 eV（对应室温）

**声子展宽**（`degaussw`）：
- 太小：谱函数有噪声
- 太大：峰展宽，丢失细节
- **推荐**：0.5-1.0 cm⁻¹

#### 8.9.3 库仑赝势μ*的确定

**理论值**：
$$
\mu^* = \frac{\mu}{1 + \mu \ln(E_F/\omega_D)}
$$

其中μ是裸库仑相互作用（~0.3-0.5）。

**实用取值**：
- 简单金属：μ* = 0.10
- 过渡金属：μ* = 0.13
- 强关联系统：μ* = 0.15-0.20

**拟合实验**：
如果有实验Tc，可以反推μ*：
```python
def find_mu_star(lambda_calc, omega_log, Tc_exp):
    """通过实验Tc反推μ*"""
    from scipy.optimize import fsolve

    def equation(mu_star):
        # McMillan公式
        Tc_calc = (omega_log / 1.2) * np.exp(
            -1.04 * (1 + lambda_calc) /
            (lambda_calc * (1 - 0.62 * mu_star) - mu_star)
        )
        return Tc_calc - Tc_exp

    mu_star_fit = fsolve(equation, 0.1)[0]
    return mu_star_fit
```

---

### 8.10 常见问题与解决方案

#### Q1: λ随展宽变化很大，如何选择？

**答**：
1. 做收敛性测试：绘制λ vs 展宽曲线
2. 找到平台区（plateau）
3. 选择平台中间的值
4. 通常0.03-0.05 eV是安全的

#### Q2: 计算的Tc比实验高/低很多？

**可能原因**：

**Tc偏高**：
- μ*取值过小 → 增大至0.13-0.15
- 晶格常数低估 → 检查优化结构
- 忽略了非谐效应

**Tc偏低**：
- μ*取值过大
- k/q网格不够密 → λ收敛性问题
- 泛函问题（PBE低估键长）

#### Q3: 如何处理虚频？

**虚频意味着结构不稳定！**

**解决方案**：
1. 检查晶格常数是否优化
2. 尝试不同压力下的结构
3. 考虑低温相变（如CDW）
4. 对于氢化物：检查H的位置是否正确

#### Q4: PDOS积分不等于总DOS？

**原因**：投影波函数不完备

**检查**：
```bash
grep "Lowdin Charges" pdos.out
```

如果Lowdin电荷与总电荷差异>5%，说明投影不好。

**解决**：
- 增加投影的轨道数
- 使用更局域的基组

---

### 8.11 案例研究：LaH10 vs H3S

对比两个典型的高温超导氢化物：

| 性质 | LaH10 (250 GPa) | H3S (150 GPa) |
|------|-----------------|---------------|
| Tc (实验) | ~250 K | ~200 K |
| λ (计算) | 2.3 | 2.1 |
| N(0) [states/(eV·spin)] | 5.5 | 4.8 |
| ωlog [K] | 1170 | 1300 |
| 主导声子 | H振动 (1000-1400 cm⁻¹) | H振动 (1200-1600 cm⁻¹) |
| 费米面 | 多口袋（La-5d + H-1s） | 简单（S-3p + H-1s） |
| 关键电子态 | La-5d (60%) + H-1s (40%) | S-3p (50%) + H-1s (50%) |

**物理启示**：

1. **共同点**：
   - 都是H振动主导
   - 都有高N(0)和强|g|
   - 都需要极高压力稳定

2. **差异**：
   - LaH10的Tc更高，因为：
     * La-5d电子更局域 → |g|更大
     * 更高的N(0)
   - H3S的ωlog更高，但λ略小
     * S比La轻，但d电子更少

**结论**：d电子丰富的重金属氢化物可能是更优的高Tc候选！

---

### 8.12 总结与展望

#### 核心要点回顾

1. **态密度N(0)的作用**：
   > N(0)越大 → λ越大 → Tc越高（不是相反！）

2. **谱函数α²F(ω)的信息**：
   - 峰位 → 主导声子频率
   - 峰强 → 耦合强度
   - 积分 → 总λ

3. **声子模式分解**：
   - 从`a2F.dosX`文件提取
   - 高频H振动通常主导
   - 前10-15个模式贡献~90%的λ

4. **电子态分析**：
   - 通过PDOS确定轨道成分
   - d电子通常比s/p电子耦合强
   - 费米面嵌套不是唯一因素

#### 分析工具箱

| 目标 | 工具/文件 | 关键输出 |
|------|-----------|----------|
| 总λ和Tc | `lambda.out`, `Tc.result` | λ, ωlog, Tc |
| 谱函数 | `ALPHA2F.OUT` | α²F(ω) |
| 声子分解 | `a2F.dos1-33` | λ_ν每个模式 |
| 轨道成分 | `projwfc.x` → PDOS | 费米面处轨道权重 |
| 费米面形状 | `fs.x` + XCrySDen | 可视化 |

#### 未来方向

1. **各向异性Eliashberg方程**：
   - 考虑$\mathbf{k}$依赖的能隙Δ(k)
   - 需要：EPW的完整矩阵元输出

2. **非谐效应**：
   - 高压氢化物中非谐性强
   - 需要：SSCHA等高级方法

3. **多能带效应**：
   - 不同能带的Δ可能不同
   - 需要：多能带Eliashberg求解器

4. **机器学习辅助**：
   - 从材料数据库预测λ
   - 高通量筛选候选材料

---

**记忆口诀**：
```
态密度大超导强，谱函数峰找声子；
声子分解看模式，PDOS揭示电子态。
费米嵌套非唯一，d电子耦合最关键。
高压氢化物Tc高，轻H振动是法宝！
```

---
