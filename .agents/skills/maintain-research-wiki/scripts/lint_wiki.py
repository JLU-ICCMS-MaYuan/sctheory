from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path


REQUIRED = {"id", "type", "title", "aliases", "status", "updated", "evidence_status"}
ALLOWED = {"source", "project", "system", "phenomenon", "physical-model", "computational-method", "ai-method", "dataset", "observable", "workflow", "software"}


def frontmatter(text: str) -> dict[str, str]:
    if not text.startswith("---\n"):
        return {}
    end = text.find("\n---\n", 4)
    if end < 0:
        return {}
    result = {}
    for line in text[4:end].splitlines():
        if ":" in line and not line.startswith((" ", "-")):
            key, value = line.split(":", 1)
            result[key.strip()] = value.strip()
    return result


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    wiki = root / "Knowledge" / "Wiki"
    issues = []
    ids = Counter()
    pages = list(wiki.rglob("*.md"))
    for page in pages:
        rel = page.relative_to(root).as_posix()
        text = page.read_text(encoding="utf-8")
        if page.name == "index.md" or "indexes" in page.parts:
            continue
        meta = frontmatter(text)
        missing = REQUIRED - set(meta)
        if missing:
            issues.append(f"FRONTMATTER {rel}: missing {sorted(missing)}")
        if meta.get("type") and meta["type"] not in ALLOWED:
            issues.append(f"TYPE {rel}: {meta['type']}")
        if meta.get("id"):
            ids[meta["id"]] += 1
        for target in re.findall(r"\[\[([^]|#]+)", text):
            if target.startswith("Journal/"):
                issues.append(f"BOUNDARY-LINK {rel}: {target}")
        for match in re.findall(r"!\[\[([^]|#]+)", text):
            candidate = root / match
            if not candidate.exists():
                issues.append(f"MISSING-ASSET {rel}: {match}")
    for item, count in ids.items():
        if count > 1:
            issues.append(f"DUPLICATE-ID {item}: {count}")
    print(f"pages={len(pages)} issues={len(issues)}")
    for issue in issues:
        print(issue)
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
