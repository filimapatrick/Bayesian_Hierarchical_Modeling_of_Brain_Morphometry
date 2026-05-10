#!/usr/bin/env python3
"""
Organize epilepsy DICOM exports into a canonical subject-based layout.

This script is tailored to the EPILEPSY1 folder structure where each top-level
folder appears to be a subject group (e.g. CC, CG, FO, GB, ID, KO, LC, NA, NJ, OE, OM, PF).

What it does:
- scans the epilepsy source root
- assigns canonical IDs sub-001, sub-002, ...
- copies each subject group's DICOM files into:
    data/raw/sub-XXX/dicom/
- writes a mapping TSV for downstream steps

It preserves the relative path under each top-level subject group so that series
folders remain intact.

Expected usage:
    python organize_epilepsy_dicom.py
    python organize_epilepsy_dicom.py --apply
    python organize_epilepsy_dicom.py --source-root /path/to/EPILEPSY1 --apply
"""

from __future__ import annotations

import argparse
import csv
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List


REPO_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry")
DEFAULT_SOURCE_ROOT = Path("/Volumes/MyHDD/ABDN_DATA/EPILEPSY1")
DEFAULT_RAW_ROOT = REPO_ROOT / "data" / "raw"
DEFAULT_MAPPING_CSV = Path("/Volumes/MyHDD/ABDN_DATA/epilepsy_subject_mapping.tsv")


@dataclass(frozen=True)
class SubjectGroup:
    canonical_subject_id: str
    source: str
    group_name: str
    source_dir: Path
    file_count: int


def iter_data_files(root: Path) -> Iterable[Path]:
    """Yield all non-sidecar files under root."""
    if not root.exists():
        return []

    for p in root.rglob("*"):
        if not p.is_file():
            continue
        if p.name == ".DS_Store" or p.name.startswith("._"):
            continue
        yield p


def discover_subject_groups(source_root: Path) -> List[Path]:
    """Return top-level folders that look like subject groups."""
    groups: List[Path] = []
    for p in sorted(source_root.iterdir()):
        if p.is_dir() and not p.name.startswith("."):
            groups.append(p)
    return groups


def assign_subject_ids(group_dirs: List[Path]) -> List[SubjectGroup]:
    """Assign sub-001, sub-002, ... in sorted top-level folder order."""
    rows: List[SubjectGroup] = []
    for idx, group_dir in enumerate(group_dirs, start=1):
        files = list(iter_data_files(group_dir))
        rows.append(
            SubjectGroup(
                canonical_subject_id=f"sub-{idx:03d}",
                source="EPILEPSY1",
                group_name=group_dir.name,
                source_dir=group_dir,
                file_count=len(files),
            )
        )
    return rows


def copy_group_files(
    row: SubjectGroup,
    raw_root: Path,
    dry_run: bool = True,
) -> int:
    """
    Copy files from one epilepsy subject group into:
        data/raw/sub-XXX/dicom/

    Relative paths under the group directory are preserved.
    """
    target_root = raw_root / row.canonical_subject_id / "dicom"

    files = list(iter_data_files(row.source_dir))
    if dry_run:
        print(
            f"[DRY-RUN] {row.canonical_subject_id} ({row.group_name}): "
            f"would copy {len(files)} files -> {target_root}"
        )
        return len(files)

    copied = 0
    for src in files:
        rel = src.relative_to(row.source_dir)
        dst = target_root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)

        if dst.exists():
            continue

        shutil.copy2(src, dst)
        copied += 1

    print(
        f"[OK] {row.canonical_subject_id} ({row.group_name}): "
        f"copied {copied} files -> {target_root}"
    )
    return copied


def write_mapping(rows: List[SubjectGroup], out_tsv: Path) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["canonical_subject_id", "source", "group_name", "file_count", "source_dir"])
        for r in rows:
            w.writerow(
                [
                    r.canonical_subject_id,
                    r.source,
                    r.group_name,
                    r.file_count,
                    str(r.source_dir),
                ]
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Organize epilepsy DICOM data into canonical subject folders.")
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT, help="Root folder containing EPILEPSY1")
    parser.add_argument("--raw-root", type=Path, default=DEFAULT_RAW_ROOT, help="Target raw data directory")
    parser.add_argument("--mapping-out", type=Path, default=DEFAULT_MAPPING_CSV, help="Output TSV mapping file")
    parser.add_argument("--apply", action="store_true", help="Actually copy files instead of dry-run")
    args = parser.parse_args()

    if not args.source_root.exists():
        print(f"[ERROR] Source root not found: {args.source_root}")
        return 2

    group_dirs = discover_subject_groups(args.source_root)
    if not group_dirs:
        print(f"[ERROR] No subject group folders found under: {args.source_root}")
        return 2

    rows = assign_subject_ids(group_dirs)
    print(f"Discovered {len(rows)} subject groups under {args.source_root}")

    write_mapping(rows, args.mapping_out)
    print(f"Saved mapping to: {args.mapping_out}")

    total_found = 0
    total_copied = 0

    for row in rows:
        total_found += row.file_count
        copied = copy_group_files(row, args.raw_root, dry_run=not args.apply)
        total_copied += copied

    print(f"\nFinished. Matched files: {total_found}. Copied: {total_copied}.")
    if not args.apply:
        print("Dry-run only. Re-run with --apply to copy files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())