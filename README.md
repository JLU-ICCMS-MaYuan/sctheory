# WorkPlaner

这是一个可直接作为 Obsidian vault 根目录使用的项目工作模板。克隆仓库后，在
Obsidian 中选择 WorkPlaner 仓库根目录作为 vault；不要把仓库嵌套复制到另一个 vault。

## 已包含

- `Journal/Daily`：每日原始日记与工作进展。
- `Journal/Weekly`、`Journal/Monthly`、`Journal/Quarterly`、`Journal/Yearly`：周期导航。
- `Templates`：周期笔记模板。
- `.obsidian`：已安装并启用的 Calendar、Periodic Notes、Templater、BRAT、Claudian、AnuPpuccin 主题，以及每日星期配色样式。
- `.claudian`：不含凭据和会话历史的 Codex Agent 初始配置。

## 使用

在 Calendar 中选择日期，或通过 Periodic Notes 命令创建周期笔记。每日文件使用
`YYYY-MM-DD` 命名，并会自动渲染日期、前后日链接和周期导航。

Daily 模板依赖 Templater 自动触发：当新文件创建在 `Journal/Daily` 中时，Templater 会
执行 `Templates/Daily-TEMPLATE.md`，生成日期、前后日链接与周期导航。核心插件 Daily
Notes 保持关闭，Daily、Weekly、Monthly、Quarterly 和 Yearly 笔记统一由 Periodic Notes
管理。

WorkPlaner 默认使用 AnuPpuccin 主题，并启用 `Daily Note Themes` 与 `General Tweaks`
CSS 片段，使每日深色页面中的内嵌表格与背景保持一致。

Daily 中的 `## 工作灵感` 用于记录尚未完成、尚未落地并有待后续评估的碎片想法；
“工作进展”用于已进入执行阶段的目标、进展与产出、问题与决策、下一步；其余个人记录
写在“日记”部分。

### 新 vault 验收

1. 完全退出旧 vault，使用仓库根目录打开新 vault；不要只复制内容子目录。
2. 确认 Calendar、Periodic Notes 和 Templater 已启用，核心 Daily Notes 已关闭；在
   Templater 设置中确认“Trigger Templater on new file creation”设为 Folder templates，且
   `Journal/Daily` 映射到 `Templates/Daily-TEMPLATE.md`。
3. 运行 `node scripts/validate-journal-config.mjs` 检查模板、插件配置和目录结构。
4. 在 Calendar 中选择一个尚无日记的日期。
5. 确认文件生成在 `Journal/Daily/YYYY-MM-DD.md`，包含 `## 工作灵感`、居中的前后日导航，
   且没有残留 `<% ... %>` 模板源码。

如果新建日记仍显示 `<% ... %>`，说明 Templater 自动触发未执行。请检查文件是否位于
`Journal/Daily`，并在 Templater 设置中确认文件夹映射；对于已创建的文件，可执行
`Templater: Replace templates in the active file` 重新渲染。

## Agent

Claudian 默认使用 Codex、中文界面、`gpt-5.6-terra` 和中等推理强度。打开右侧栏的
Claudian 视图即可开始新会话；BRAT 已订阅 `YishenTu/claudian`，用于后续更新插件。

WorkPlaner 不携带 API Key、Base URL、设备 ID、CLI 绝对路径或历史会话。首次在一台新机器
使用前，需要安装 Codex CLI 并通过 Codex 自身完成认证；不要把密钥写入可提交文件。
`.gitignore` 会排除 Claudian 本地设置、会话记录和 Obsidian 工作区状态。

## 依赖版本

- Obsidian 1.13.0 或更高版本
- Templater 2.24.3
- Periodic Notes 0.0.17
- Calendar 1.5.10
- BRAT 2.2.0
- Claudian 2.1.2（桌面端）
- Codex CLI（由用户环境提供）

Obsidian 只读取 vault 根目录的 `.obsidian` 配置。若必须将 WorkPlaner 并入现有项目，
应先退出 Obsidian 并备份目标 `.obsidian`，再把仓库内容合并到目标 vault 根目录；重新打开
后执行上述验收，避免运行中的旧插件设置覆盖新配置。
