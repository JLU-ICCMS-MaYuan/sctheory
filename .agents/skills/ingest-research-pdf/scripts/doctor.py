from __future__ import annotations

import importlib.util
import json
import os
import shutil
import sys
from pathlib import Path


def present(module: str) -> bool:
    return importlib.util.find_spec(module) is not None


def command_present(name: str) -> bool:
    if shutil.which(name):
        return True
    root = Path(os.environ.get("CODEX_PRIMARY_RUNTIME", Path.home() / ".cache/codex-runtimes/codex-primary-runtime"))
    return (root / "dependencies/native/poppler/Library/bin" / f"{name}.exe").is_file()


report = {
    "python": sys.version.split()[0],
    "base": {
        "pypdf": present("pypdf"),
        "pdfplumber": present("pdfplumber"),
        "PIL": present("PIL"),
        "pdfinfo": command_present("pdfinfo"),
        "pdftoppm": command_present("pdftoppm"),
    },
    "optional": {
        "pymupdf4llm": present("pymupdf4llm"),
    },
}
report["recommended_backend"] = "native"
report["network_or_ocr_required"] = False
report["base_ready"] = all(report["base"].values())
print(json.dumps(report, ensure_ascii=False, indent=2))
raise SystemExit(0 if report["base_ready"] else 2)
