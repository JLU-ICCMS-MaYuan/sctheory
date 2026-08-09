# PDF extraction output schema

Each `Knowledge/Extracted/<source-id>/` contains:

- `source.md`: metadata, section map, extracted page text, and proposed knowledge items.
- `pages/page-NNNN.md`: page-scoped Markdown with stable PDF-page anchors.
- `figures.md`: figure crops, bilingual captions, page locators, and supported claims.
- `tables.md`: table crops and optional faithfully reconstructed Markdown tables.
- `equations.md`: equation crops, original equation text when recoverable, symbols, and assumptions.
- `references.md`: extracted bibliography and identifiers.
- `extraction-report.md`: tool versions, coverage, warnings, and unresolved items.
- `manifest.json`: hashes, page dimensions, artifacts, and crop boxes.
- `assets/`: local PNG crops excluded from ordinary Git.

Frontmatter requires `id`, `type`, `source_file`, `source_sha256`, `title`, `aliases`, `status`, `extracted_at`, `page_count`, and `evidence_status`.

The manifest records the local backend and tool versions when available. Optional PyMuPDF4LLM Markdown never overwrites the pypdf/pdfplumber/Poppler evidence baseline.

Use PDF page numbers as `pdf_page`. Preserve printed page labels separately when detectable.
