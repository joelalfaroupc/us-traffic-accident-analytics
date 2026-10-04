"""Restore local academic data without redistributing dataset files."""
import argparse
import shutil
import zipfile
from pathlib import Path


def restore(archive, stage, root):
    with zipfile.ZipFile(archive) as z:
        for info in z.infolist():
            parts = Path(info.filename).parts
            if "__MACOSX" in parts or not parts:
                continue
            name = parts[-1]
            if name.startswith(".") or Path(name).suffix.lower() not in {".csv", ".rds"}:
                continue
            # Fixed basename only: archive paths cannot escape local data directory.
            destination = root / "data" / stage / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            content = z.read(info)
            if destination.exists() and destination.read_bytes() != content:
                # Some stages contain duplicate basenames. Preserve first occurrence.
                print("Skipped different duplicate:", info.filename)
                continue
            destination.write_bytes(content)
            datasets = root / "workspace" / stage / "Datasets"
            datasets.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(destination, datasets / name)
            if name in {"US_Accidents_Final_FE.csv", "US_Accidents_Final_FE.rds"}:
                shutil.copyfile(destination, root / "workspace" / stage / name)
                shutil.copyfile(destination, root / "data" / name)
    print("Restored", stage, "inputs; missing intermediate files must still be regenerated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--d3", type=Path, required=True)
    parser.add_argument("--d4", type=Path, required=True)
    options = parser.parse_args()
    repository = Path(__file__).resolve().parents[1]
    restore(options.d3, "d3", repository)
    restore(options.d4, "d4", repository)
