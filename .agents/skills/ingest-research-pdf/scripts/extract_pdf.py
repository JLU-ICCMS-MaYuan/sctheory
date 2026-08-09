from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pdfplumber
from PIL import Image
from pypdf import PdfReader


CAPTION_RE = re.compile(r"^(fig(?:ure)?\.?|extended data fig(?:ure)?\.?|table)\s*([s\d][\w.-]*)", re.I)
EQUATION_RE = re.compile(r"\(([A-Z]?\d{1,3})\)\s*$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def slug(value: str) -> str:
    value = value.lower().replace("/", "-")
    value = re.sub(r"[^a-z0-9._-]+", "-", value)
    return re.sub(r"-+", "-", value).strip("-") or "source"


def lines_with_boxes(page) -> list[dict]:
    words = page.extract_words(use_text_flow=True, keep_blank_chars=False)
    groups: dict[int, list[dict]] = defaultdict(list)
    for word in words:
        groups[round(float(word["top"]) / 3) * 3].append(word)
    lines = []
    for _, group in sorted(groups.items()):
        ordered = sorted(group, key=lambda item: float(item["x0"]))
        lines.append({
            "text": " ".join(item["text"] for item in ordered),
            "x0": min(float(item["x0"]) for item in ordered),
            "x1": max(float(item["x1"]) for item in ordered),
            "top": min(float(item["top"]) for item in ordered),
            "bottom": max(float(item["bottom"]) for item in ordered),
        })
    return lines


def find_poppler(name: str) -> str:
    root = Path(os.environ.get("CODEX_PRIMARY_RUNTIME", Path.home() / ".cache/codex-runtimes/codex-primary-runtime"))
    candidate = root / "dependencies/native/poppler/Library/bin" / f"{name}.exe"
    if candidate.is_file():
        return str(candidate)
    found = shutil.which(name)
    if found:
        return found
    raise FileNotFoundError(f"{name} not found; run doctor.py")


def render_page(pdf: Path, page_number: int, output: Path, dpi: int) -> Path:
    prefix = output / f"page-{page_number:04d}"
    subprocess.run([
        find_poppler("pdftoppm"), "-f", str(page_number), "-l", str(page_number),
        "-r", str(dpi), "-png", "-singlefile", str(pdf), str(prefix)
    ], check=True, capture_output=True)
    return prefix.with_suffix(".png")


def crop_box(kind: str, line: dict, width: float, height: float) -> tuple[float, float, float, float]:
    if kind == "table":
        return (0, max(0, line["top"] - 24), width, min(height, line["top"] + height * 0.55))
    if kind == "equation":
        return (0, max(0, line["top"] - 48), width, min(height, line["bottom"] + 48))
    return (0, max(0, line["top"] - height * 0.55), width, min(height, line["bottom"] + 44))


def md_escape(value: str) -> str:
    return value.replace("\n", " ").replace('"', "'").strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("pdf", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args()

    pdf = args.pdf.resolve()
    output = args.output.resolve()
    assets = output / "assets"
    output.mkdir(parents=True, exist_ok=True)
    assets.mkdir(parents=True, exist_ok=True)

    reader = PdfReader(str(pdf))
    metadata = {str(k).lstrip("/"): str(v) for k, v in (reader.metadata or {}).items()}
    source_hash = sha256(pdf)
    title = metadata.get("Title") or pdf.stem
    now = datetime.now(timezone.utc).isoformat()
    pages_data = []
    artifacts = []
    references_lines = []
    rendered: dict[int, Path] = {}

    with pdfplumber.open(str(pdf)) as document:
        for page_index, page in enumerate(document.pages, start=1):
            text = page.extract_text(x_tolerance=2, y_tolerance=3) or ""
            lines = lines_with_boxes(page)
            if re.search(r"\bReferences\b|参考文献", text, re.I):
                references_lines.append(f"\n## PDF page {page_index}\n\n{text}")
            candidates = []
            for line in lines:
                stripped = line["text"].strip()
                caption = CAPTION_RE.match(stripped)
                if caption:
                    kind = "table" if caption.group(1).lower().startswith("table") else "figure"
                    candidates.append((kind, caption.group(2), line, stripped))
                else:
                    equation = EQUATION_RE.search(stripped)
                    if equation and len(stripped) > 5:
                        candidates.append(("equation", equation.group(1), line, stripped))

            for candidate_index, (kind, label, line, evidence_text) in enumerate(candidates, start=1):
                if page_index not in rendered:
                    rendered[page_index] = render_page(pdf, page_index, assets, args.dpi)
                box = crop_box(kind, line, float(page.width), float(page.height))
                scale = args.dpi / 72.0
                pixel_box = tuple(round(value * scale) for value in box)
                filename = f"{kind}-p{page_index:04d}-{slug(label)}-{candidate_index:02d}.png"
                target = assets / filename
                with Image.open(rendered[page_index]) as image:
                    image.crop(pixel_box).save(target)
                artifacts.append({
                    "kind": kind,
                    "label": label,
                    "pdf_page": page_index,
                    "caption": evidence_text,
                    "crop_box_pdf_points": list(box),
                    "asset": f"assets/{filename}",
                    "sha256": sha256(target),
                    "confidence": "needs-review" if kind == "equation" else "base-detected",
                })
            pages_data.append({"pdf_page": page_index, "width": page.width, "height": page.height, "text": text})

    for temp_page in rendered.values():
        temp_page.unlink(missing_ok=True)

    manifest = {
        "schema_version": 1,
        "source_id": args.source_id,
        "source_file": pdf.name,
        "source_path_at_ingest": str(pdf),
        "source_sha256": source_hash,
        "title": title,
        "metadata": metadata,
        "page_count": len(pages_data),
        "extracted_at": now,
        "dpi": args.dpi,
        "artifacts": artifacts,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    source_lines = [
        "---", f"id: {args.source_id}", "type: source-record", f"source_file: \"{pdf.name}\"",
        f"source_sha256: {source_hash}", f"title: \"{md_escape(title)}\"", "aliases: []",
        "status: extracted", f"extracted_at: {now}", f"page_count: {len(pages_data)}",
        "evidence_status: needs-review", "---", "", f"# {title}", "", "## 元数据（Metadata）", "",
        f"- 主题（Subject）: {metadata.get('Subject', 'Unknown')}", f"- 作者（Author）: {metadata.get('Author', 'Unknown')}",
        f"- SHA-256: `{source_hash}`", "", "## 全文逐页提取（Full-text extraction）", ""
    ]
    for page in pages_data:
        source_lines.extend([f"### PDF page {page['pdf_page']}", "", page["text"], ""])
    (output / "source.md").write_text("\n".join(source_lines), encoding="utf-8")

    for kind in ("figure", "table", "equation"):
        rows = [item for item in artifacts if item["kind"] == kind]
        md = [f"# {kind.title()} evidence", ""]
        for item in rows:
            md.extend([
                f"## {kind.title()} {item['label']} - PDF page {item['pdf_page']}", "",
                f"![[{item['asset']}]]", "", f"- 原文（Original）: {item['caption']}",
                f"- 裁剪区域（Crop box）: `{item['crop_box_pdf_points']}`", f"- 状态（Status）: `{item['confidence']}`", ""
            ])
        (output / f"{kind}s.md").write_text("\n".join(md), encoding="utf-8")

    (output / "references.md").write_text("# References\n" + "\n".join(references_lines), encoding="utf-8")
    report = [
        "# Extraction report", "", f"- Source ID: `{args.source_id}`", f"- Pages: {len(pages_data)}",
        f"- Figures detected: {sum(a['kind'] == 'figure' for a in artifacts)}",
        f"- Tables detected: {sum(a['kind'] == 'table' for a in artifacts)}",
        f"- Equation candidates: {sum(a['kind'] == 'equation' for a in artifacts)}",
        "- Mode: base (`pypdf` + `pdfplumber` + Poppler)",
        "- Network access: not used",
        "- OCR/model inference: not used",
        "- Optional PyMuPDF4LLM: not used by the base extractor",
        "- Publication gate: blocked until human review", "", "## Required review", "",
        "- Verify bilingual title and DOI metadata.", "- Inspect every figure/table crop before citing it.",
        "- Equation candidates are heuristic and must be checked against the rendered PDF.",
        "- Confirm references section boundaries.",
    ]
    (output / "extraction-report.md").write_text("\n".join(report), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
