# WorkPlaner Snapshot Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Do not dispatch subagents unless the user explicitly authorizes delegation.

**Goal:** Import WorkPlaner as a history-free snapshot at the `sctheory` repository root, reorganize the existing research workspace under `Research/`, retain tracked source PDFs, and preserve only the existing `sctheory` remotes and Git history.

**Architecture:** WorkPlaner becomes the Obsidian vault shell at repository root. Research code, calculations, notebooks, and project documents move to four independent directories under `Research/`; evidence PDFs move to `Knowledge/Sources/PDF/`; external Git submodules remain at `external/`. WorkPlaner's `Knowledge/Wiki/Projects/` remains reserved for schema-compliant knowledge pages and is not used for research code or data.

**Tech Stack:** Git, Obsidian/WorkPlaner, Markdown, Node.js validation, existing Python/Mathematica/research assets.

---

## Locked Decisions

- Import WorkPlaner files only; do not merge or retain its Git history.
- Use a temporary shallow clone solely to obtain the snapshot, then remove its `.git` metadata and temporary directory.
- Do not add a persistent `workplaner` remote.
- Preserve `origin` and `gitee` exactly as currently configured for `sctheory`.
- Replace the root `AGENTS.md` and `README.md` with WorkPlaner versions.
- Delete `CLAUDE.md`, `GEMINI.md`, and `request.md`.
- Keep `external/` and all entries in `.gitmodules` at their existing paths.
- Track PDFs under `Knowledge/Sources/PDF/`; do not keep WorkPlaner's PDF ignore rule.
- Do not create compatibility symlinks or retain old top-level research paths.
- Preserve the user's pre-existing `.vscode/settings.json` worktree modification and exclude it from every migration commit.

## Target Structure

```text
sctheory/
├── .agents/                          # WorkPlaner skills
├── .obsidian/                        # WorkPlaner vault configuration
├── AGENTS.md                         # WorkPlaner version
├── Bases/
├── Journal/
├── Knowledge/
│   ├── Sources/PDF/                  # tracked original PDFs
│   ├── Extracted/
│   ├── Wiki/Projects/                # schema-compliant Wiki pages only
│   ├── Templates/
│   ├── RAW-SPEC.md
│   ├── WIKI-SPEC.md
│   ├── COMPILE-SPEC.md
│   └── log.md
├── Research/
│   ├── molecular-orbital/
│   │   ├── H3S/
│   │   ├── MH6/
│   │   ├── MH9/
│   │   ├── MH10/
│   │   ├── MY2H24/
│   │   ├── scripts/
│   │   └── 分子轨道mathematica代码/
│   ├── quantum-geometry-epc/
│   │   ├── quantum_geometry/
│   │   └── learning_ppt/
│   ├── chemical-bonding/
│   │   ├── mathematica_GroupTheory1.4/
│   │   ├── reference-notes/
│   │   ├── 04_LaH10_CeH10_research_plan.md
│   │   └── QUICKSTART.md
│   └── migdal-eliashberg-theory/
│       ├── EliashbergEquation/
│       ├── 01_fortran_to_python_guide.md
│       ├── 02_questions_and_answers.md
│       └── 03_theory_comprehensive.md
├── Templates/
├── external/                         # unchanged submodule paths
│   └── docs/
│       ├── deeptb_cookbook.md
│       └── deeptb_training_troubleshooting.md
├── scripts/                          # WorkPlaner validation scripts only
├── .gitignore
├── .gitmodules
└── README.md                         # WorkPlaner version
```

### Task 1: Capture a Reproducible Preflight Baseline

**Files:**
- Read: `.git/config`
- Read: `.gitmodules`
- Read: `.vscode/settings.json`
- Create temporarily: `/tmp/workplaner-import.XXXXXX/`

- [ ] **Step 1: Verify branch, worktree, remotes, and current commit**

Run:

```bash
git status --short --branch
git remote -v
git rev-parse HEAD
git submodule status
```

Expected: branch `master`; only the known `.vscode/settings.json` user modification may be present; `origin` and `gitee` both still target `sctheory`; submodule paths remain under `external/`.

- [ ] **Step 2: Record the WorkPlaner snapshot commit without configuring a remote**

Run:

```bash
git ls-remote --symref "git@gitee.com:mayuan_JLUPHY/WorkPlaner.git" HEAD
```

Expected: `HEAD` resolves to `refs/heads/master` and prints a commit SHA. Record that SHA in the migration commit body.

- [ ] **Step 3: Obtain an isolated shallow snapshot**

Run:

```bash
workplaner_tmp="$(mktemp -d -p /tmp workplaner-import.XXXXXX)"
git clone --depth 1 "git@gitee.com:mayuan_JLUPHY/WorkPlaner.git" "$workplaner_tmp/WorkPlaner"
git -C "$workplaner_tmp/WorkPlaner" rev-parse HEAD
```

Expected: the cloned SHA matches Step 2. Do not add this repository as a `sctheory` remote.

### Task 2: Move Existing Research Assets into Their Final Ownership Boundaries

**Files:**
- Move: `H3S/`, `MH6/`, `MH9/`, `MH10/`, `MY2H24/`, `scripts/`, `分子轨道mathematica代码/`
- Move: `quantum_geometry/`, `learning_ppt/`
- Move: `mathematica_GroupTheory1.4/`, `reference/`
- Move: `EliashbergEquation/`, all six files currently under `docs/`
- Move: `deeptb_cookbook.md`
- Delete: `CLAUDE.md`, `GEMINI.md`, `request.md`

- [ ] **Step 1: Create only the agreed destination directories**

Run:

```bash
mkdir -p "Research/molecular-orbital"
mkdir -p "Research/quantum-geometry-epc"
mkdir -p "Research/chemical-bonding"
mkdir -p "Research/migdal-eliashberg-theory"
mkdir -p "Knowledge/Sources/PDF"
mkdir -p "external/docs"
```

- [ ] **Step 2: Move the molecular-orbital research workspace**

Run each path explicitly:

```bash
git mv "H3S" "Research/molecular-orbital/H3S"
git mv "MH6" "Research/molecular-orbital/MH6"
git mv "MH9" "Research/molecular-orbital/MH9"
git mv "MH10" "Research/molecular-orbital/MH10"
git mv "MY2H24" "Research/molecular-orbital/MY2H24"
git mv "scripts" "Research/molecular-orbital/scripts"
git mv "分子轨道mathematica代码" "Research/molecular-orbital/分子轨道mathematica代码"
```

Expected: no old top-level directory from this list remains.

- [ ] **Step 3: Move the quantum-geometry and chemical-bonding workspaces**

Run:

```bash
git mv "quantum_geometry" "Research/quantum-geometry-epc/quantum_geometry"
git mv "learning_ppt" "Research/quantum-geometry-epc/learning_ppt"
git mv "mathematica_GroupTheory1.4" "Research/chemical-bonding/mathematica_GroupTheory1.4"
git mv "reference" "Research/chemical-bonding/reference-notes"
```

Expected: each directory has exactly one destination; the duplicated `learning_ppt/` entry from the original request does not create duplicate content.

- [ ] **Step 4: Move PDFs from notes into the WorkPlaner source layer**

Run:

```bash
find "Research/chemical-bonding/reference-notes" -maxdepth 1 -type f -iname '*.pdf' -print0 |
while IFS= read -r -d '' pdf_path; do
  git mv "$pdf_path" "Knowledge/Sources/PDF/"
done
```

Expected: `Knowledge/Sources/PDF/` contains the formerly tracked `reference/*.pdf` files; Markdown, Python, and text files remain in `Research/chemical-bonding/reference-notes/`.

- [ ] **Step 5: Move theory and project documentation**

Run:

```bash
git mv "EliashbergEquation" "Research/migdal-eliashberg-theory/EliashbergEquation"
git mv "docs/01_fortran_to_python_guide.md" "Research/migdal-eliashberg-theory/01_fortran_to_python_guide.md"
git mv "docs/02_questions_and_answers.md" "Research/migdal-eliashberg-theory/02_questions_and_answers.md"
git mv "docs/03_theory_comprehensive.md" "Research/migdal-eliashberg-theory/03_theory_comprehensive.md"
git mv "docs/04_LaH10_CeH10_research_plan.md" "Research/chemical-bonding/04_LaH10_CeH10_research_plan.md"
git mv "docs/QUICKSTART.md" "Research/chemical-bonding/QUICKSTART.md"
git mv "docs/deeptb_training_troubleshooting.md" "external/docs/deeptb_training_troubleshooting.md"
git mv "deeptb_cookbook.md" "external/docs/deeptb_cookbook.md"
```

- [ ] **Step 6: Delete the explicitly retired root documents**

Run:

```bash
git rm "CLAUDE.md" "GEMINI.md" "request.md"
```

Expected: the files are absent from the worktree but remain recoverable from earlier `sctheory` commits.

### Task 3: Import the WorkPlaner Snapshot without Its Git Identity

**Files:**
- Replace: `AGENTS.md`, `README.md`, `.gitignore`
- Create: `.agents/`, `.obsidian/`, `Bases/`, `Journal/`, `Knowledge/` template/spec files, `Templates/`, `scripts/`
- Preserve: `.git/`, `.gitmodules`, `external/`, `.vscode/settings.json`

- [ ] **Step 1: Remove Git identity from the temporary snapshot**

Run after validating the explicit temporary path from Task 1:

```bash
rm -rf "$workplaner_tmp/WorkPlaner/.git"
```

Expected: only the temporary clone's `.git` directory is removed. Never target the current repository's `.git` directory.

- [ ] **Step 2: Copy the WorkPlaner snapshot into the repository root**

Run:

```bash
rsync -a "$workplaner_tmp/WorkPlaner/" "./"
```

Expected: WorkPlaner `AGENTS.md`, `README.md`, `.gitignore`, vault configuration, templates, journal placeholders, knowledge specifications, and validation script exist at root. Current `.gitmodules`, `external/`, and `.vscode/` remain present because the snapshot has no conflicting entries.

- [ ] **Step 3: Confirm no WorkPlaner remote or history was introduced**

Run:

```bash
git remote -v
git log --all --oneline --decorate -5
```

Expected: only `origin` and `gitee` are configured; no `workplaner` remote exists; current ancestry is still exclusively the `sctheory` history.

### Task 4: Merge Ignore Rules and Repair Active References

**Files:**
- Modify: `.gitignore`
- Modify: `Research/chemical-bonding/04_LaH10_CeH10_research_plan.md`
- Modify: `Research/chemical-bonding/QUICKSTART.md`
- Modify: `Research/molecular-orbital/scripts/LaH10_CeH10_analysis/README.md`
- Modify: `Research/molecular-orbital/MH6/deeptb_training/README.md`
- Inspect and modify where operational: `Research/molecular-orbital/MH6/deeptb_training/data_generation/*.py`, `Research/molecular-orbital/MH6/verify_with_pythtb.py`
- Inspect only unless safely editable: `.nb`, generated `.out`, `.fls`, and archived calculation outputs

- [ ] **Step 1: Make WorkPlaner's `.gitignore` retain project-specific ignores and track PDFs**

Use `apply_patch` to remove:

```gitignore
Knowledge/Sources/PDF/**/*.pdf
```

Append these rewritten rules:

```gitignore
# Research workspace generated data.
Research/molecular-orbital/MH6/deeptb_training/h6_pure/data
Research/molecular-orbital/MH6/deeptb_training/h6_pure/output
Research/molecular-orbital/MH6/deeptb_training/h6_pure/train1
Research/molecular-orbital/MH6/deeptb_training/h6_pure/train2
Research/molecular-orbital/MH6/deeptb_training/h6_pure/train3
Research/quantum-geometry-epc/quantum_geometry/graphene/runs
```

Expected: `git check-ignore "Knowledge/Sources/PDF/<existing-pdf-name>.pdf"` exits nonzero, while WorkPlaner local state and generated extraction assets remain ignored.

- [ ] **Step 2: Update active Markdown links and commands**

Use `rg` to find old repository-root references:

```bash
rg -n --hidden --glob '!.git/**' --glob '!external/*/.git/**' \
  '(^|[^[:alnum:]_])(EliashbergEquation|H3S|MH6|MH9|MH10|MY2H24|quantum_geometry|learning_ppt|mathematica_GroupTheory1\.4|分子轨道mathematica代码|docs|reference|scripts)/' \
  "Research" "external/docs" "Knowledge"
```

Apply these canonical mappings only where the path is operational:

```text
EliashbergEquation/               -> Research/migdal-eliashberg-theory/EliashbergEquation/
H3S/                              -> Research/molecular-orbital/H3S/
MH6/                              -> Research/molecular-orbital/MH6/
MH9/                              -> Research/molecular-orbital/MH9/
MH10/                             -> Research/molecular-orbital/MH10/
MY2H24/                           -> Research/molecular-orbital/MY2H24/
scripts/LaH10_CeH10_analysis/     -> Research/molecular-orbital/scripts/LaH10_CeH10_analysis/
quantum_geometry/                 -> Research/quantum-geometry-epc/quantum_geometry/
learning_ppt/                     -> Research/quantum-geometry-epc/learning_ppt/
mathematica_GroupTheory1.4/       -> Research/chemical-bonding/mathematica_GroupTheory1.4/
reference/<paper>.pdf             -> Knowledge/Sources/PDF/<paper>.pdf
docs/01_fortran_to_python_guide.md -> Research/migdal-eliashberg-theory/01_fortran_to_python_guide.md
docs/02_questions_and_answers.md   -> Research/migdal-eliashberg-theory/02_questions_and_answers.md
docs/03_theory_comprehensive.md    -> Research/migdal-eliashberg-theory/03_theory_comprehensive.md
docs/04_LaH10_CeH10_research_plan.md -> Research/chemical-bonding/04_LaH10_CeH10_research_plan.md
docs/QUICKSTART.md                -> Research/chemical-bonding/QUICKSTART.md
```

Do not alter historical paths embedded in raw calculation output merely to silence the search. Report them as historical artifacts. Do not create old-path symlinks.

- [ ] **Step 3: Verify source and configuration files contain no actionable old root paths**

Run:

```bash
rg -n --glob '*.py' --glob '*.sh' --glob '*.json' --glob '*.yaml' --glob '*.yml' --glob '*.md' \
  '(^|[^[:alnum:]_])(EliashbergEquation|H3S|MH6|MH9|MH10|MY2H24|quantum_geometry|learning_ppt|mathematica_GroupTheory1\.4|分子轨道mathematica代码|docs|reference)/' \
  "Research" "Knowledge" "external/docs"
```

Expected: remaining matches are descriptive historical text or paths internal to a moved directory; no executable command depends on an old repository-root path.

### Task 5: Validate the Vault, Git Boundaries, Moves, and Research Entry Points

**Files:**
- Test: `scripts/validate-journal-config.mjs`
- Test: `.gitmodules`
- Test: moved Python entry points
- Test: all staged paths

- [ ] **Step 1: Validate WorkPlaner's Obsidian configuration**

Run:

```bash
node "scripts/validate-journal-config.mjs"
```

Expected: exit code 0 with templates, plugins, and journal directory structure accepted.

- [ ] **Step 2: Verify required paths and forbidden legacy paths**

Run:

```bash
test -f "AGENTS.md"
test -f "README.md"
test -d "Research/molecular-orbital"
test -d "Research/quantum-geometry-epc"
test -d "Research/chemical-bonding"
test -d "Research/migdal-eliashberg-theory"
test -d "Knowledge/Sources/PDF"
test -d "external/docs"
test ! -e "H3S"
test ! -e "MH6"
test ! -e "MH9"
test ! -e "MH10"
test ! -e "MY2H24"
test ! -e "quantum_geometry"
test ! -e "learning_ppt"
test ! -e "EliashbergEquation"
test ! -e "docs"
test ! -e "reference"
test ! -e "CLAUDE.md"
test ! -e "GEMINI.md"
test ! -e "request.md"
```

Expected: all checks exit 0. The implementation plan itself must be moved from `docs/superpowers/plans/` to `Research/plans/` before asserting that `docs/` is absent.

- [ ] **Step 3: Move this plan to its durable post-migration location**

Run:

```bash
mkdir -p "Research/plans"
git mv "docs/superpowers/plans/2026-08-09-workplaner-migration.md" "Research/plans/2026-08-09-workplaner-migration.md"
```

Expected: the plan remains tracked, and the obsolete root `docs/` hierarchy can disappear.

- [ ] **Step 4: Verify submodule paths and status are unchanged**

Run:

```bash
git config --file ".gitmodules" --get-regexp 'submodule\..*\.path'
git submodule status
```

Expected paths: `external/pythtb`, `external/DeePTB`, `external/dftio`, and `external/HamEPC`; no submodule is relocated or absorbed into the parent repository.

- [ ] **Step 5: Run the smallest safe research smoke checks**

Run:

```bash
python "Research/migdal-eliashberg-theory/EliashbergEquation/eliashberg_solver.py" --help
python -m compileall -q \
  "Research/migdal-eliashberg-theory/EliashbergEquation" \
  "Research/molecular-orbital/scripts" \
  "Research/molecular-orbital/MH6/deeptb_training/data_generation"
```

Expected: solver help exits 0; Python sources compile without syntax errors. Do not run expensive VASP, QE, Mathematica, or DeePTB training jobs during migration validation.

- [ ] **Step 6: Audit PDF tracking and ignored local state**

Run:

```bash
git ls-files "Knowledge/Sources/PDF/*.pdf"
git check-ignore -v ".obsidian/workspace.json" ".claudian/sessions/example.json" "Knowledge/Extracted/example/assets/example.png"
```

Expected: all migrated PDFs are listed by `git ls-files`; local WorkPlaner state and generated extraction assets are ignored.

### Task 6: Stage Only the Migration and Create One Atomic Commit

**Files:**
- Stage: all WorkPlaner snapshot files, agreed moves, deletions, path repairs, and this relocated plan
- Exclude: `.vscode/settings.json` and all unrelated pre-existing worktree changes

- [ ] **Step 1: Review the complete change set and detect accidental data loss**

Run:

```bash
git status --short
git diff --stat
git diff --summary
git diff --check
```

Expected: Git reports renames where possible; deleted root documents are limited to `CLAUDE.md`, `GEMINI.md`, and `request.md`; the root `README.md`, `AGENTS.md`, and `.gitignore` are replacements; `.vscode/settings.json` remains unstaged.

- [ ] **Step 2: Stage explicit migration path groups**

Run:

```bash
git add ".agents" ".obsidian" "Bases" "Journal" "Knowledge" "Research" "Templates" "scripts"
git add "external/docs" "AGENTS.md" "README.md" ".gitignore"
git add -u -- "H3S" "MH6" "MH9" "MH10" "MY2H24" "quantum_geometry" "learning_ppt" "EliashbergEquation" "mathematica_GroupTheory1.4" "分子轨道mathematica代码" "reference" "docs" "CLAUDE.md" "GEMINI.md" "request.md" "deeptb_cookbook.md"
```

Do not use `git add .` or `git add -A`.

- [ ] **Step 3: Confirm the user-owned editor change is excluded**

Run:

```bash
git diff --cached --name-only
git diff --cached -- ".vscode/settings.json"
```

Expected: `.vscode/settings.json` is absent from the staged name list and has no staged diff.

- [ ] **Step 4: Commit using the repository-mandated structured message**

Construct the exact timestamp and staged statistics immediately before committing, then commit with this structure:

```text
[时间] YYYY-MM-DD HH:MM:SS
[修改目的] 部署 WorkPlaner 并重组科研工作区
[修改内容]
- WorkPlaner: 导入 Obsidian vault 配置、模板、Journal 与 Knowledge 规范
- Research: 按四个研究方向迁移现有代码、数据、Notebook 与项目文档
- Knowledge/Sources/PDF: 迁移并继续跟踪原始参考文献
- external/docs: 归档 DeePTB 使用经验，保持子模块路径不变
- 根目录: 使用 WorkPlaner 的 AGENTS.md、README.md 和合并后的 .gitignore
[影响范围] Obsidian vault、科研目录路径、文档链接、脚本入口和参考文献位置
[验证方式] WorkPlaner 配置校验、路径断言、子模块检查、Python 编译与 Git 跟踪审计
已更改 X 个文件, XX 行插入(+), XX 行删除(-)
```

- [ ] **Step 5: Verify the committed result and retained remotes**

Run:

```bash
git status --short --branch
git remote -v
git show --stat --oneline --decorate HEAD
```

Expected: the migration is committed; `.vscode/settings.json` remains as the user's pre-existing unstaged change; only `origin` and `gitee` exist; no push is performed.

- [ ] **Step 6: Remove only the validated temporary snapshot directory**

Resolve and print the exact `workplaner_tmp` path, confirm it is under `/tmp/workplaner-import.*`, then remove that explicit directory. Report that this temporary copy is not recoverable but all imported content is now tracked in the `sctheory` commit.

## Acceptance Criteria

- WorkPlaner opens from the `sctheory` root as an Obsidian vault and passes `node scripts/validate-journal-config.mjs`.
- `sctheory` retains its original commit ancestry and its `origin`/`gitee` remotes.
- No WorkPlaner remote, `.git` metadata, or Git ancestry remains in the repository.
- All agreed research assets live under `Research/`; no old top-level compatibility paths exist.
- `external/` submodules and `.gitmodules` paths are unchanged.
- Root `AGENTS.md` and `README.md` are WorkPlaner versions.
- `CLAUDE.md`, `GEMINI.md`, `request.md`, and the old research `README.md` content are not retained in the new tree.
- Existing and future PDFs under `Knowledge/Sources/PDF/` are trackable by default.
- Active scripts and documents use new paths; historical paths in immutable/generated outputs are reported rather than rewritten.
- The user's pre-existing `.vscode/settings.json` modification is untouched and excluded from the migration commit.
