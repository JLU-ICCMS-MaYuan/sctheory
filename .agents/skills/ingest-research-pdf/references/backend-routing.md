# Local backend routing

## Decision order

1. Accept only the explicitly selected local PDF.
2. Use pypdf for metadata and page accounting.
3. Use pdfplumber for text, coordinates, captions, table candidates, and evidence anchors.
4. Use Poppler for page rendering and figure/table/equation crops.
5. Use PyMuPDF4LLM only after explicit user selection and only when already installed.

Never retrieve a replacement source or send the PDF over the network. Preserve extraction under `Knowledge/Extracted` and require human approval before wiki compilation.

## Base backend

The base backend is authoritative for page numbering, coordinates, hashes, and rendered evidence. It is model-free and must remain the default.

## Optional PyMuPDF4LLM

Use it only to improve Markdown reading order for born-digital, multi-column PDFs. Do not treat its Markdown as stronger evidence than the base page coordinates and crops. Record the package version and retain base extraction for comparison.

Do not auto-install it. PyMuPDF/PyMuPDF4LLM licensing must be reviewed before redistribution or use in a distributed product.

## Unsupported pages

If a page has no usable text, mark it `needs-review` and retain a Poppler rendering. Do not add OCR, infer missing formulae, or reconstruct unreadable table cells.

## Quality gate

Block compilation when any of these remain:

- missing pages;
- blank or severely garbled native text;
- unresolved formula or table placeholders;
- Markdown references a remote asset;
- source hash or local backend provenance is absent.
