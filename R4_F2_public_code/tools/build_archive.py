"""Build a code-only local review archive. Never upload or publish anything."""
from pathlib import Path
import hashlib
import json
import re
import zipfile

root = Path(__file__).resolve().parents[1]
destination = root.parent / "R4_F2_results"
destination.mkdir(exist_ok=True)
directories = ("R", "metadata", "docs", "environment", "tests", "tools")
allowed = {root / "README.md", root / "run.R", root / "config.R", root / ".gitignore", root / "inputs" / "README.md"}
for directory in directories:
    for path in (root / directory).rglob("*"):
        if path.is_file() and path.suffix in {".R", ".md", ".csv", ".json", ".txt", ".py"}:
            allowed.add(path)
paths = sorted(allowed, key=lambda path: path.relative_to(root).as_posix())
if not all(path.is_file() for path in paths):
    raise RuntimeError("Missing release file")
inventory = []
for path in paths:
    relative = path.relative_to(root).as_posix()
    if "config.local" in relative or "__pycache__" in relative:
        raise RuntimeError(f"Private/generated file in archive: {relative}")
    data = path.read_bytes()
    content = data.decode("utf-8-sig")
    if re.search(r"[A-Za-z]:[\\/](?:Users|paperfanxiu)[\\/]", content, re.IGNORECASE):
        raise RuntimeError(f"Personal absolute path in public file: {relative}")
    inventory.append({"file": relative, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
archive = destination / "R4_F2_public_code.zip"
with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as output:
    for path in paths:
        output.write(path, root.name + "/" + path.relative_to(root).as_posix())
with zipfile.ZipFile(archive) as check:
    if check.testzip() is not None:
        raise RuntimeError("Archive integrity check failed")
    for entry, item in zip(check.infolist(), inventory):
        if hashlib.sha256(check.read(entry)).hexdigest() != item["sha256"]:
            raise RuntimeError(f"Archive file hash mismatch: {entry.filename}")
manifest = {"scope": "Local replacement prepared 2026-10-20; original Figure 6 definitions; remote repository is not updated by this tool", "files": inventory,
            "archive": archive.name, "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest()}
(destination / "archive_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
print(f"Verified {len(paths)} code/documentation files; archive size {archive.stat().st_size} bytes")
print(archive)
