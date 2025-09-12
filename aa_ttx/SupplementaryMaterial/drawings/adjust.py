#!/usr/bin/env python3
import re
import sys
from pathlib import Path

FNAME_RE = re.compile(r'^GL(\d+)\.dot$')
# Accept optional spaces, and optional sign in the parenthetical integer.
LABEL_RE = re.compile(
    r'^(?P<pre>\s*label\s*=\s*"FSG\s+GL\d+\s+#)'
    r'(?P<num>\d+)'
    r'(?P<post>\s+x\s*\(\s*[-+]?\d+\s*\)"\s*;)'
    r'(?P<eol>\r?\n)?$'
)

def extract_j(p: Path):
    m = FNAME_RE.match(p.name)
    return int(m.group(1)) if m else None

def main():
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path.cwd()
    files = [p for p in root.iterdir() if p.is_file() and FNAME_RE.match(p.name)]
    if not files:
        print("No GL<J>.dot files found", file=sys.stderr)
        sys.exit(1)

    files.sort(key=lambda p: extract_j(p))

    print("Ordered files (by J):")
    for i, p in enumerate(files, start=1):
        print(f"{i:4d}  {p.name}")

    for i, p in enumerate(files, start=1):
        lines = p.read_text(encoding="utf-8").splitlines(keepends=True)
        if len(lines) < 2 or not lines[1].lstrip().startswith("label"):
            raise ValueError(f"{p.name}: expected the label as the second line")
        #new_label = LABEL_RE.sub(lambda m: f"{m.group(1)}{i}{m.group(3)}", lines[1], count=1)
        new_label = LABEL_RE.sub(
           lambda m: f"{m.group('pre')}{i}{m.group('post')}{m.group('eol') or ''}",
            lines[1],
            count=1,
        )
        if new_label == lines[1]:
            raise ValueError(f"{p.name}: label line did not match expected format:\n{lines[1]}")
        lines[1] = new_label
        p.write_text("".join(lines), encoding="utf-8")
        print(f"Updated {p.name}: set # to {i}")

if __name__ == "__main__":
    main()

