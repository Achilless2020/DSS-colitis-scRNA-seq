"""Verify that the distributable contains a runnable entry point and no private inputs."""
from pathlib import Path
import subprocess
import sys
import zipfile

root = Path(__file__).resolve().parents[1]
subprocess.run([sys.executable, str(root / "tools/build_archive.py")], check=True)
with zipfile.ZipFile(root.parent / "R4_F2_results/R4_F2_public_code.zip") as archive:
    files = archive.namelist()
    assert f"{root.name}/run.R" in files, "Public ZIP is missing the documented CLI entry point"
    assert f"{root.name}/config.R" in files
    assert not any("config.local" in name or "/outputs/" in name for name in files)
    assert not any(name.endswith((".rds", ".RData", ".mtx.gz")) for name in files)
print("PASS: public archive has CLI/configuration and excludes private inputs/outputs")
