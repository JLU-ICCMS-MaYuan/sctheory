---
name: maintain-research-wiki
description: Compile, query, index, and audit an Obsidian research wiki for condensed-matter computation and AI-for-science. Use after PDF extraction is approved, when saving evidence-backed knowledge, asking what the wiki contains, rebuilding Dataview indexes, checking terminology, detecting conflicts or orphan pages, or linting Knowledge/Wiki. Never reads Journal and never writes from a plain query or lint without approval.
---

# Maintain Research Wiki

This is a vendored and substantially adapted derivative of Astro-Han/karpathy-llm-wiki under MIT. It has no runtime dependency on that repository or another skill.

## Non-negotiable boundaries

- Read sources only from `Knowledge/Sources/**` and approved evidence from `Knowledge/Extracted/**`.
- Never read or search `Journal/**`.
- Human wiki edits are authoritative. Never silently overwrite them.
- Query and lint are read-only by default.
- Show a create/update/conflict/deprecate plan and obtain approval before publication or repair.
- Append every approved write to `Knowledge/log.md`.

## Knowledge types

Assign exactly one primary `type`; represent other relationships as frontmatter links:

1. `source`: one source record.
2. `project`: one research question.
3. `system`: material or physical system.
4. `phenomenon`: physical behavior.
5. `physical-model`: theoretical representation.
6. `computational-method`: non-AI computational method.
7. `ai-method`: data-driven model or agent method.
8. `dataset`: reusable data asset and provenance.
9. `observable`: measurable or computable quantity.
10. `workflow`: ordered reproducible procedure.
11. `software`: executable tool or infrastructure.

Read `references/taxonomy.md` when classification is ambiguous. Never create a broad catch-all AI directory outside `AI-Methods/`.

## Compile

1. Read `Knowledge/RAW-SPEC.md`, `Knowledge/WIKI-SPEC.md`, and `Knowledge/COMPILE-SPEC.md`.
2. Require an approved `Knowledge/Extracted/<source-id>/extraction-report.md`; an explicitly requested end-to-end test may publish draft-only pages.
3. Read `Knowledge/Wiki/indexes/terminology.md` first.
4. Search IDs, titles, English names, abbreviations, and aliases across the full wiki.
5. Classify the source as New, Update, Disputed, or No material.
6. Draft a change plan. Compile sources sequentially using `Knowledge/Templates/`.
7. Give every claim a `claim_type` (`reported`, `derived`, `hypothesis`, `conflict`) and evidence locator.
8. Preserve older claims and every `human:lock` block; never silently rewrite history.
9. Update indexes and append the log only after all pages validate.

## Query

Read `Knowledge/Wiki/index.md`, then full-text search using Chinese terms, English names, abbreviations, and aliases. Answer from the wiki with links. Do not write unless the user explicitly asks to save and approves the plan.

## Lint

Run `python scripts/lint_wiki.py <project-root>`. Report deterministic and semantic findings without fixes. Check links, IDs, fields, evidence locators, terminology, local assets, orphans, conflicts, unsupported claims, and the `Journal/**` boundary. Semantic fixes require individual review.

## DataviewJS

Dataview and DataviewJS may read `Knowledge/**` only. They must not access the network, execute commands, edit files, or scan `Journal/**`. Markdown/frontmatter remains the source of truth.

## Language

Write Chinese-first prose. Preserve English titles and captions in parentheses. Preserve formulas. Define symbols and units bilingually. Use one canonical Chinese title and English/abbreviation aliases.

## Upstream notice

Read `references/UPSTREAM-LICENSE.md` when redistributing or updating this skill.
