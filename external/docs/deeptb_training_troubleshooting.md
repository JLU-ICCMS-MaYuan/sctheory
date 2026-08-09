# DeePTB 训练常见问题与解决方案

本文档记录了在运行 `dptb train` 命令时遇到的一系列配置错误及其解决方案。通过逐一解决这些问题，最终成功使训练任务得以运行。

---

### 问题 1: `AssertionError` - 损失函数方法错误

- **错误信息**: `AssertionError` on `... in ["eigvals"]`
- **根本原因**: 在 `input.json` 文件中，`train_options.loss_options.train` 下的 `method` 参数被设置为一个不被当前训练模式所支持的值（例如 `"skints"`）。代码断言此处的值必须是 `["eigvals"]` 之一。
- **解决方案**: 修改 `input.json`，将 `loss_options.train` 块的内容替换为符合要求的配置。
  ```json
  // "train_options": {
  //   ...
      "loss_options": {
          "train": {
              "method": "eigvals",
              "diff_on": false,
              "eout_weight": 0.001,
              "diff_weight": 0.01
          }
      },
  //   ...
  // }
  ```

---

### 问题 2: `AssertionError` - 未找到数据路径

- **错误信息**: `AssertionError: No trajectory folders are found. Please check the prefix.`
- **根本原因**: `input.json` 中 `data_options.train.root` 指定的路径是相对于**命令执行的当前目录**的。如果路径设置不正确，程序将无法找到数据。例如，如果数据在 `h6_pure/data` 中，而 `root` 仅设置为 `"./data"`，则会导致路径错误。
- **解决方案**: 修正 `input.json` 中的 `root` 路径，确保它能从命令执行目录正确地指向数据所在的文件夹。在此案例中，我们将路径从 `"./data"` 修改为 `"h6_pure/data"`。

---

### 问题 3: `ArgumentKeyError` - 数据集 `info.json` 包含非法键

- **错误信息**: 
    - `ArgumentKeyError: undefined key 'description' is not allowed in strict mode`
    - `ArgumentKeyError: undefined key 'parameters' is not allowed in strict mode`
- **根本原因**: `DeePTB` 在严格模式下解析数据集的元信息文件（例如 `data/set.0/info.json`）。该文件包含了程序 schema 中未定义的键（如 `description`, `parameters`），从而导致校验失败。
- **解决方案**: 编辑 `info.json` 文件，直接移除 `description` 和 `parameters` 这两个不被允许的键值对。

---

### 问题 4: `KeyError: 's'` - 基组名称不匹配

- **错误信息**: `KeyError: 's'`
- **根本原因**: `input.json` 的 `common_options.basis` 中为氢原子（'H'）指定的基组是 `["s"]`。然而，`DeePTB` 的一个内部数据库 `dptb/nn/sktb/onsiteDB.py` 中，为氢原子 s 轨道在位能定义的键名是 `'1s'`。由于键名不匹配，导致查找失败。
- **解决方案**: 修改 `input.json`，将氢原子的基组从 `["s"]` 改为 `["1s"]`，以匹配内部数据库的期望格式。
  ```json
  // "common_options": {
  //   ...
      "basis": {
          "H": ["1s"]
      },
  //   ...
  // }
  ```
