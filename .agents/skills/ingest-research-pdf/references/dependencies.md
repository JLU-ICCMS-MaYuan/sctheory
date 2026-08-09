# Capability tiers

## Base

Require Python 3.10+, `pypdf`, `pdfplumber`, Pillow, Poppler `pdfinfo`, and `pdftoppm`. Extract text, metadata, page geometry, captions, and PNG crops.

## Optional multi-column Markdown

PyMuPDF4LLM is optional and must be installed only after explicit user selection:

```powershell
& ".venv/Scripts/python.exe" -m pip install "pymupdf4llm"
```

- Use it for born-digital, multi-column reading order only.
- Keep pypdf/pdfplumber/Poppler output as the evidence baseline.
- Review PyMuPDF licensing before distributing the project or service.

It does not make OCR part of the default workflow.

## Deliberately unsupported

Do not use Docling, Marker, MinerU, OCRmyPDF, Tesseract, Mistral, Mathpix, or other OCR/model services in this workflow. Report scanned or unreadable pages for human review.

Install allowed capabilities in the project-root `.venv`. Run `scripts/doctor.py` to report capabilities. Missing PyMuPDF4LLM is informational; missing base tools blocks extraction.
