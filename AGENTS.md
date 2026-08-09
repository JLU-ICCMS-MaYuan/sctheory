# WorkPlaner Agent Instructions

## Project skills

- `.agents/skills/ingest-research-pdf/SKILL.md`: extract and stage research PDFs.
- `.agents/skills/maintain-research-wiki/SKILL.md`: compile, query, and lint the research wiki.

## Hard boundaries

- Never read, search, summarize, ingest, or cite `Journal/**` for research-knowledge tasks.
- Only files explicitly placed in `Knowledge/Sources/**` are eligible source material.
- Treat `Knowledge/Sources/PDF/**` as immutable.
- PDF extraction writes only to `Knowledge/Extracted/**`; publishing to `Knowledge/Wiki/**` is a separate confirmed operation.
- Wiki queries and lint are read-only unless the user explicitly approves a write plan.
- Human edits in `Knowledge/Wiki/**` are authoritative and must not be silently overwritten.

## Language and evidence

- Write Chinese-first prose with English originals in parentheses.
- Preserve equations exactly. Explain symbols as Chinese meaning (English definition).
- Every load-bearing scientific claim needs a source ID and page/figure/table/equation locator.
- Distinguish `reported`, `derived`, `hypothesis`, and `conflict` claims.
- Consult `Knowledge/Wiki/indexes/terminology.md` before creating or renaming scientific terms.

## Knowledge schema

Before publishing extracted material, read these files in order:

1. `Knowledge/RAW-SPEC.md`
2. `Knowledge/WIKI-SPEC.md`
3. `Knowledge/COMPILE-SPEC.md`

Use the templates under `Knowledge/Templates/`. Do not invent an incompatible page shape.
