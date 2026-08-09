# deeptb使用帮助手册

## 第1步: 从结构中提取sk参数
需要准备提取sk参数的输入文件. 使用到的命令是`dptb esk sk_in.json -m poly4 `, 需要用到的输入文件的模板为：比如：
1. si
```json
{
    "common_options": {
        "basis": {
            "Si": ["3s","3p","d*"]
        }
    }
}
```

2. GaAs
```json
{
    "common_options": {
        "basis": {
            "Ga": ["s","p","d"],
            "As": ["s","p","d"]
        }
    }
}
```

## 第2步: 绘制sk参数的能带和DFT能带对比
需要用到的命令是：`dptb run band.json -i sktb.json -o  band -stu ../data/silicon.vasp`

## 第3步: 产生deeptb的训练的输入文件
需要用到的命令是：`dptb config ./ -tr -sk -m ./sktb.json`. 运行完成后修改input_templete.json里面关于路径的参数
特别注意：
```shell
- r_max: 
    取 hopping.rs + 3*w（若方法是 powerlaw/varTang96 则用 rs + 8*w）。你设的 rs=1.2, w=0.35，所以 1.2 + 3*0.35 = 2.25，对应 {'H-H':
2.25}。
- oer_max: 
    取 onsite.rs + 3*w（strain/NRL + powerlaw/varTang96 用 rs + 8*w）。你设的 rs=1.2, w=0.3，所以 1.2 + 3*0.3 ≈ 2.1。
- er_max: 
    只有用 env 校正/embedding 时才会有；当前 nnsk 配置没给 env cutoff，所以是 None。
```
## 第4步: 训练钱准备好info.json

## 第5步：训练

## 附录: model_options.nnsk 方法详解

在 `model_options.nnsk` 中，`onsite` 和 `hopping` 的 `method` 决定了模型如何描述在位能校正和跳跃积分。

### 1. model_options.nnsk.onsite.method (在位能校正)
用于定义原子在位能（矩阵对角块）的修正方式：

- **`none`**: (默认) 不进行额外校正。直接使用 DeePTB 数据库中预设的隔离原子能量。
- **`uniform`**: 在数据库隔离原子能量的基础上，为每种原子的每种轨道（s, p, d...）增加一个可训练的统一偏移量 $\epsilon_l^\prime$。
- **`uniform_noref`**: 与 `uniform` 类似，但不加数据库的参考值，完全由模型从头拟合。
- **`strain`**: **(推荐)** 考虑周围原子对当前原子产生的场（应变效应）。使用类似于 Slater-Koster 的积分形式来修正 onsite 矩阵，考虑 cutoff 范围内的原子贡献。需要设置 `rs` 和 `w`。
- **`NRL`**: 使用 NRL-TB (Naval Research Laboratory) 的在位能公式。

### 2. model_options.nnsk.hopping.method (跳跃积分公式)
用于定义 Slater-Koster 跳跃积分随键长 $r$ 的衰减函数 $f(r)$：

- **`powerlaw`**: 幂律衰减 $1/r^n$。公式简洁，物理意义明确。
- **`poly1pow` ~ `poly4pow`**: **(常用)** 结合多项式和幂律的形式。例如 `poly4pow` 通常能提供更精细的拟合效果。
- **`poly2exp` ~ `poly4exp`**: 结合多项式和指数衰减的形式。
- **`varTang96`**: 基于 Tang 等人 1996 年提出的公式的变体。
- **`NRL0` / `NRL1`**: 使用 NRL-TB 提供的跳跃项拟合公式。
- **`custom`**: 用户自定义函数。

### 3. 参数 rs, w 与 Cutoff 的计算逻辑
对于绝大多数方法（除 NRL 系列外）：
- **`rs`**: 截断平滑函数 $f_c(r)$ 取值为 0.5 的位置。
- **`w`**: 衰减宽度。
- **实际截断距离**:
    - **`r_max` (Hopping)**: 
        - 若 method 为 `powerlaw` 或 `varTang96`：`r_max = rs + 8 * w`
        - 其他方法：`r_max = rs + 3 * w`
    - **`oer_max` (Onsite)**:
        - 仅在 `strain` 或 `NRL` 模式下生效。
        - 规则同上：`onsite.rs + 3*w` (或 8*w)。

> **提示**: 在 `input.json` 中设置 `rs` 和 `w` 时，务必确保 `data_options` 中的 `r_max` 和 `oer_max` 覆盖了上述计算出的理论截断距离，否则会导致超出范围的原子对被遗漏。


