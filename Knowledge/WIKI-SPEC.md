# 科研 Wiki 规范（Research Wiki Specification）

Wiki 是面向人类阅读、由 Agent 维护、以证据约束的“活知识”。它不是 PDF 摘要堆栈；每个页面应跨来源累积，并随着查询、审核和新论文持续改善。

## 十一种互斥页面类型

每页只能选择一个中心对象作为 `type`：

| type | 目录 | 中心对象 |
|---|---|---|
| source | `Sources/` | 一份来源的忠实科研摘要 |
| project | `Projects/` | 一个科研问题 |
| system | `Systems/` | 材料、结构、缺陷、表面或异质结 |
| phenomenon | `Phenomena/` | 物理行为或效应 |
| physical-model | `Physical-Models/` | 理论模型或有效哈密顿量 |
| computational-method | `Computational-Methods/` | 非 AI 数值或理论方法 |
| ai-method | `AI-Methods/` | ML 模型、训练策略、RAG 或 Agent 方法 |
| dataset | `Datasets/` | 有谱系的数据集合 |
| observable | `Observables/` | 可计算或可测量的量 |
| workflow | `Workflows/` | 有输入、步骤、产物和验收标准的过程 |
| software | `Software/` | 软件、服务、调度器或基础设施 |

材料页面提到 EPC，并不使它同时成为 phenomenon；用 `phenomena` 关系连接。禁止通过复制段落解决分类问题。

## 页面生命周期

```text
draft → reviewed → stable
  └────→ disputed
reviewed/stable → superseded
```

- `draft`：Agent 已编译，尚未完成人类审核。
- `reviewed`：人类检查了表达、证据和分类。
- `stable`：可作为后续综合的可靠基线。
- `disputed`：来源之间存在未决冲突。
- `superseded`：保留历史，但由更准确页面或结论替代。

Agent 不能自行把页面晋升为 `reviewed` 或 `stable`。

## 人类阅读体验

知识页按以下阅读顺序组织：

1. **一句话认识**：快速说明页面讲什么。
2. **为什么重要**：与凝聚态物理或当前项目的关系。
3. **物理图像**：先给直觉，再给数学形式。
4. **核心主张**：每条主张标记类型和证据。
5. **定量结果**：数值、单位、参数和比较基准。
6. **方法与适用边界**：近似、计算设置、失败条件。
7. **图表与公式**：嵌入必要截图，不把页面变成图片墙。
8. **知识连接**：链接到体系、现象、模型、方法和项目。
9. **开放问题**：未知、冲突、下一步来源。
10. **人类审核**：审核状态、修改意见和受保护内容。

页面开头适合快速阅读，证据和审核信息置于后部；不要求读者先穿过 frontmatter 或日志才能理解物理内容。

## 主张规范

每条关键主张使用紧凑格式：

```markdown
- **[reported]** 量子几何对 EPC 提供非平凡贡献。
  证据：[[来源页]]；PDF p.4，Fig. 2；状态：`verified`
```

允许的 `claim_type`：

- `reported`：来源直接陈述或展示。
- `derived`：从已列出的证据推导，必须展示推导输入。
- `hypothesis`：尚未验证的解释。
- `conflict`：来源不一致，必须并列双方证据。

## 人类修改权

人工内容优先于 Agent。需要防止覆盖的段落使用：

```markdown
<!-- human:lock reason="保留课题组判断" -->
这里是人类确认或修订的内容。
<!-- /human:lock -->
```

Agent 可以指出锁定内容与新证据冲突，但不得修改锁定块。所有人工审核记录追加到“人类审核”章节，并由 Git 保存完整历史。

## Frontmatter

```yaml
---
id: stable-id
type: phenomenon
title: 电子-声子耦合
aliases: [electron-phonon coupling, EPC]
status: draft
evidence_status: needs-review
sources: []
projects: []
systems: []
phenomena: []
physical_models: []
computational_methods: []
ai_methods: []
datasets: []
observables: []
workflows: []
software: []
created: YYYY-MM-DD
updated: YYYY-MM-DD
reviewed_by: []
reviewed_at:
---
```

`title` 使用中文规范名；英文名称和缩写进入 `aliases`。写入前必须查询 `indexes/terminology.md`。

## 导航与概览

- `index.md`：DataviewJS 动态目录，只读展示。
- `overview.md`：人类可读的活综合，只有新来源改变全局认识时才更新。
- `indexes/review-queue.md`：自动汇总 draft、needs-review 和 disputed 页面。
- `log.md`：追加记录 ingest、compile、query、review、health、lint。

概览页不能成为所有页面都链接的“重力井”；具体知识应直接互链。
