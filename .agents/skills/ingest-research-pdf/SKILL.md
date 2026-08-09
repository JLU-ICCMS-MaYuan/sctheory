---
name: ingest-research-pdf
description: Self-contained local research-PDF ingestion for condensed-matter computation and AI-for-science. Use when asked to inspect, split, extract, reprocess, or quality-check user-provided main papers and supplementary PDFs under Knowledge/Sources/PDF. Uses pypdf, pdfplumber, and Poppler without OCR models or network retrieval; PyMuPDF4LLM is an explicit optional backend. Produces auditable Markdown and evidence under Knowledge/Extracted; never publishes directly to the wiki.
---

# Ingest Research PDF

Perform PDF-only, evidence-preserving local extraction. Do not invoke or require another skill.

## Boundary

- Read only explicitly selected PDFs under `Knowledge/Sources/PDF/`.
- Use only PDF files supplied or explicitly selected by the user. Do not search for arXiv, DOI, publisher HTML, source data, or replacement copies.
- Never read `Journal/**`.
- Never modify source PDFs.
- Write only to `Knowledge/Extracted/<source-id>/`.
- Do not upload PDF content or call cloud parsing/OCR services.
- Stop after extraction. `maintain-research-wiki` owns publication.
- Process one source at a time unless the user explicitly approves a batch.

## Workflow

1. Show the selected PDF path, SHA-256, size, page count, and planned output directory. State that processing is local and model-free.
2. Verify that the selected file is a readable PDF. Reject URLs and non-PDF substitutes.
3. Prefer the project `.venv` interpreter. Run `python scripts/doctor.py` before extraction.
4. Run `python scripts/extract_pdf.py <pdf> --output <directory> --source-id <id>` for the base path: pypdf metadata, pdfplumber text/coordinates and candidates, and Poppler rendering/crops.
5. Use PyMuPDF4LLM only when the user explicitly selects it and it is already installed. Keep it optional because of its PyMuPDF licensing implications; never auto-install or silently switch to it.
6. Inspect representative rendered crops visually. For a long SI, inspect the first, middle, and last pages plus pages flagged by extraction heuristics.
7. Read `extraction-report.md`; never repair uncertainty with guesses. Mark scanned or unreadable pages `needs-review` because this workflow has no OCR.
8. Produce Chinese-first notes using English originals in parentheses. Preserve equations, units, identifiers, and source wording.
9. Present the extraction result and wait for explicit approval before wiki compilation.

## Evidence contract

Use claim types `reported`, `derived`, `hypothesis`, and `conflict`. Every scientific claim identifies the stable source ID and at least one page, figure, table, equation, or section anchor. Unknown locators remain `needs-review`.

Figures and tables use Obsidian embeds pointing to local assets. Keep original captions in English and add Chinese translations. Never reconstruct an unreadable table cell.

## Language

- Chinese prose first; retain English titles and terminology in parentheses.
- Preserve formula bodies exactly.
- Define symbols as Chinese meaning (English definition).
- Format units as symbol: Chinese name (English name).
- Use the repository terminology registry during later compilation.

## Environment

- Install enhanced dependencies only in the project-root `.venv`.
- Never modify the public/system Python when `.venv` exists.
- Do not commit `.venv` or downloaded model caches.
- Do not install or invoke Docling, Marker, MinerU, Tesseract, cloud OCR clients, or local OCR models.
- Missing OCR is an intentional limitation, not a dependency error.

## Progressive references

- Read `references/output-schema.md` before changing output formats.
- Read `references/dependencies.md` when doctor reports missing capabilities.
- Read `references/backend-routing.md` before changing local backend selection or quality thresholds.
