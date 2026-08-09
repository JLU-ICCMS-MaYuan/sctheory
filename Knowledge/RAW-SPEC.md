# Raw 证据层规范（Raw Evidence Specification）

Raw 层保存事实来源和可再生成的拆解物。它不是 Wiki，也不追求顺畅叙事；它追求忠实、定位准确和可重复生成。

## 目录模型

```text
Knowledge/
├── Sources/PDF/               # 不可变原始 PDF，本地保存，Git 忽略
└── Extracted/<source-id>/      # PDF 的结构化拆解结果
    ├── source.md
    ├── figures.md
    ├── tables.md
    ├── equations.md
    ├── references.md
    ├── extraction-report.md
    ├── manifest.json
    └── assets/                 # 截图，本地保存，Git 忽略
```

这里的 Raw 是逻辑层，由 `Sources/PDF` 与 `Extracted` 共同组成；不再额外创建含义重复的根目录 `raw/`。

## 不变量

1. 原始 PDF 只读，禁止覆盖、重排或清理。
2. `manifest.json` 必须记录 PDF SHA-256、页数、工具版本和每个截图的哈希及裁剪坐标。
3. 提取文本不得伪装成作者排版原文；多栏顺序、公式和表格存在不确定性时必须标注。
4. 原文陈述与 Agent 解释分开保存。
5. 找不到页码或图表编号的内容不能进入稳定 Wiki 主张。
6. `Journal/**` 不是来源，禁止进入 Raw 层。

## 主文与 SI

主文和 Supplementary Information 使用不同 `source-id`，通过以下字段配对：

```yaml
document_role: main | supplementary
work_id: doi-10.1038-s41567-024-02486-0
companions:
  - doi-10.1038-s41567-024-02486-0-si
```

主文负责论点与主要结果；SI 负责推导、参数、收敛、补充图表和可复现性。Wiki 引用必须指向真正承载证据的那一份文档，不能只引用主文 DOI 代替 SI 页码。

## 提取状态

```text
queued → extracted → needs-review → approved
                    ↘ rejected
```

- `extracted`：工具运行完成，不代表正确。
- `needs-review`：存在需检查的公式、表格、阅读顺序或裁剪。
- `approved`：人类允许进入 Wiki 编译流程。
- `rejected`：提取质量不足，必须重跑或人工修订。

## 人类审核 Raw

审核者优先检查：标题与 DOI、主文/SI 配对、章节顺序、所有拟引用图表、关键公式、定量结论所在页、计算参数和限制条件。人工修订后，在 `extraction-report.md` 追加审核记录，不重写历史记录。
