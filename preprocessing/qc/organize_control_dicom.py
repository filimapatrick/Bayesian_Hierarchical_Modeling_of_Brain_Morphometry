#!/usr/bin/env python3
"""
Organize control DICOM data into a canonical subject-based layout.

Expected source layout:
  /Volumes/MyHDD/Controls_from_Rivers_State_University_Teaching_Hospital/control_rsuth/
      Ahmad Mohammed/
      Anubo Korie/
      ...
Each person folder may contain:
  - a nested folder such as AHM/, or
  - DICOM files and/or nested DICOM folders directly.

This script copies each top-level person folder into:
  data/raw/control/sub-XXX/dicom/
and preserves the internal folder structure under dicom/.
"""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


REPO_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry")
DEFAULT_SOURCE_ROOT = Path("/Volumes/MyHDD/Controls_from_Rivers_State_University_Teaching_Hospital/control_rsuth")
DEFAULT_RAW_ROOT = REPO_ROOT / "data" / "raw" / "control"
DEFAULT_MAPPING_TSV = Path("/Volumes/MyHDD/ABDN_DATA/control_subject_mapping.tsv")


def is_data_file(path: Path) -> bool:
    return path.is_file() and path.name != ".DS_Store" and not path.name.startswith("._")


def collect_data_files(root: Path) -> list[Path]:
    return sorted(
        [p for p in root.rglob("*") if is_data_file(p)],
        key=lambda p: str(p),
    )


def discover_subject_folders(source_root: Path) -> list[Path]:
    if not source_root.exists():
        return []
    return [
        p for p in sorted(source_root.iterdir())
        if p.is_dir() and not p.name.startswith(".")
    ]


def copy_subject(source_dir: Path, subject_id: str, raw_root: Path, dry_run: bool) -> int:
    files = collect_data_files(source_dir)
    target_root = raw_root / subject_id / "dicom"

    if dry_run:
        print(f"[DRY-RUN] {subject_id} ({source_dir.name}): would copy {len(files)} files -> {target_root}")
        return len(files)

    target_root.mkdir(parents=True, exist_ok=True)

    copied = 0
    for src in files:
        rel = src.relative_to(source_dir)
        dst = target_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists():
            continue
        shutil.copy2(src, dst)
        copied += 1

    print(f"[OK] {subject_id} ({source_dir.name}): copied {copied} files -> {target_root}")
    return copied


def write_mapping(rows: list[dict[str, str]], mapping_tsv: Path) -> None:
    mapping_tsv.parent.mkdir(parents=True, exist_ok=True)
    with mapping_tsv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "canonical_subject_id",
                "source_name",
                "source_path",
                "file_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Organize control DICOM folders into canonical subject directories.")
    parser.add_argument(
        "--source-root",
        type=Path,
        default=DEFAULT_SOURCE_ROOT,
        help="Root containing control person folders",
    )
    parser.add_argument(
        "--raw-root",
        type=Path,
        default=DEFAULT_RAW_ROOT,
        help="Destination raw root, e.g. .../data/raw/control",
    )
    parser.add_argument(
        "--mapping-tsv",
        type=Path,
        default=DEFAULT_MAPPING_TSV,
        help="Where to save the mapping TSV",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Actually copy files instead of dry-run",
    )
    args = parser.parse_args()

    subject_dirs = discover_subject_folders(args.source_root)
    if not subject_dirs:
        print(f"[WARN] No subject folders found under {args.source_root}")
        return 0

    print(f"Discovered {len(subject_dirs)} control subject folders under {args.source_root}")

    rows: list[dict[str, str]] = []
    total_files = 0
    total_copied = 0

    for idx, source_dir in enumerate(subject_dirs, start=1):
        subject_id = f"sub-{idx:03d}"
        files = collect_data_files(source_dir)
        total_files += len(files)

        rows.append(
            {
                "canonical_subject_id": subject_id,
                "source_name": source_dir.name,
                "source_path": str(source_dir),
                "file_count": str(len(files)),
            }
        )

        total_copied += copy_subject(
            source_dir=source_dir,
            subject_id=subject_id,
            raw_root=args.raw_root,
            dry_run=not args.apply,
        )

    write_mapping(rows, args.mapping_tsv)
    print(f"Saved mapping to: {args.mapping_tsv}")
    print(f"\nFinished. Matched files: {total_files}. Copied: {total_copied}.")
    if not args.apply:
        print("Dry-run only. Re-run with --apply to copy files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())