# quantum_geometry PPT

该目录存放论文 `Non-trivial quantum geometry and the strength of electron-phonon coupling` 的教学型 `LaTeX Beamer` 幻灯片源码。

## 文件

- `main.tex`：主幻灯片源码，使用 `ctexbeamer`，当前版本为 36 页推导型教学讲稿。

## 编译

优先使用 `XeLaTeX`：

```bash
xelatex "main.tex"
```

如果本机安装了 `latexmk`，也可以使用：

```bash
latexmk -xelatex "main.tex"
```

## 说明

- 当前版本默认使用“图示占位框”，没有直接嵌入论文原图。
- 每页都包含简短的讲者备注，默认设置为 `hide notes`，不会显示在最终 PDF 中。
- 内容基于主文与 SI 的教学性重组，重点突出：
  - 背景
  - 研究现状
  - 方法中的连续推导
  - 结论
- 当前版本将 SI 中的关键推导主干直接放入正文，而不是只做概述。

## 已知限制

- 当前工作环境中未检测到 `xelatex`/`latexmk`，因此本次只完成源码生成，未在本机做编译验收。
- 如需进一步提升可展示性，下一步可补：
  - 论文原图抽取与裁切
  - 更细致的公式排版微调
  - 标题页个人署名信息
