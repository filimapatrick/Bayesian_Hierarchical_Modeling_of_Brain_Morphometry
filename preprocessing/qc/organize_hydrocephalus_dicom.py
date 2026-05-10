#!/usr/bin/env python3
"""
Organize hydrocephalus DICOM exports into a canonical subject-based layout.

This script handles the messy Hydrocephalus source tree:

/Volumes/MyHDD/ABDN_DATA/Hydrocephalus/
  hydrocephalus_data_bmh/
    AR/
    FG/
    ...
  hydrocephalus_data_upth/
    AzS/
    CGlo/
    ...

It discovers leaf folders that contain DICOM-like files, assigns canonical
subject IDs (sub-001, sub-002, ...), and copies them to:

data/raw/hydrocephalus/sub-XXX/dicom/
"""

from __future__ import annotations

import argparse
import csv
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry")
DEFAULT_SOURCE_ROOT = Path("/Volumes/MyHDD/ABDN_DATA/Hydrocephalus")
DEFAULT_RAW_ROOT = REPO_ROOT / "data" / "raw" / "hydrocephalus"
DEFAULT_MAPPING_TSV = Path("/Volumes/MyHDD/ABDN_DATA/hydrocephalus_subject_mapping.tsv")


@dataclass(frozen=True)
class SourceGroup:
    site: str
    group_dir: Path
    file_count: int


def is_data_file(path: Path) -> bool:
    return path.is_file() and path.name != ".DS_Store" and not path.name.startswith("._")


def dir_has_data_files(root: Path) -> bool:
    if not root.exists():
        return False
    return any(is_data_file(p) for p in root.rglob("*"))


def iter_directories(root: Path) -> list[Path]:
    if not root.exists():
        return []
    return sorted((p for p in root.rglob("*") if p.is_dir()), key=lambda p: (len(p.parts), str(p)))


def discover_site_roots(source_root: Path) -> list[Path]:
    if not source_root.exists():
        return []
    return [p for p in sorted(source_root.iterdir()) if p.is_dir() and not p.name.startswith(".")]


def site_label_from_name(name: str) -> str:
    lower = name.lower()
    if "bmh" in lower:
        return "BMH"
    if "upth" in lower:
        return "UPTH"
    return name.upper()


def discover_group_dirs(site_root: Path) -> list[Path]:
    """
    Find leaf-ish directories that contain data files.

    A directory qualifies if:
    - it contains data files somewhere under it, and
    - it does not have a child directory that also contains data files.

    This works for:
    - BMH: site/AR/AR_02_A/
    - UPTH: site/AzS/
    """
    candidates: list[Path] = []
    dirs = [site_root] + iter_directories(site_root)

    for d in dirs:
        if not dir_has_data_files(d):
            continue
        if any(child.is_dir() and dir_has_data_files(child) for child in d.iterdir()):
            continue
        candidates.append(d)

    # Remove duplicates while preserving order
    seen = set()
    unique: list[Path] = []
    for d in candidates:
        if d in seen:
            continue
        seen.add(d)
        unique.append(d)
    return unique


def collect_data_files(root: Path) -> list[Path]:
    return sorted(
        [
            p
            for p in root.rglob("*")
            if is_data_file(p)
        ],
        key=lambda p: str(p),
    )


def load_existing_subjects(mapping_tsv: Path) -> set[str]:
    if not mapping_tsv.exists():
        return set()
    subjects: set[str] = set()
    with mapping_tsv.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sid = row.get("canonical_subject_id", "").strip()
            if sid:
                subjects.add(sid)
    return subjects


def write_mapping(rows: list[dict[str, str]], mapping_tsv: Path) -> None:
    mapping_tsv.parent.mkdir(parents=True, exist_ok=True)
    with mapping_tsv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "canonical_subject_id",
                "site",
                "source_group",
                "source_group_path",
                "file_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def copy_group_files(
    subject_id: str,
    group_dir: Path,
    raw_root: Path,
    dry_run: bool = True,
) -> int:
    files = collect_data_files(group_dir)
    target_dir = raw_root / subject_id / "dicom"

    if dry_run:
        print(f"[DRY-RUN] {subject_id}: would copy {len(files)} files -> {target_dir}")
        return len(files)

    target_dir.mkdir(parents=True, exist_ok=True)
    copied = 0
    for src in files:
        rel = src.relative_to(group_dir)
        dst = target_dir / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists():
            continue
        shutil.copy2(src, dst)
        copied += 1

    print(f"[OK] {subject_id}: copied {copied} files -> {target_dir}")
    return copied


def main() -> int:
    parser = argparse.ArgumentParser(description="Organize hydrocephalus DICOM data into canonical subject folders.")
    parser.add_argument(
        "--source-root",
        type=Path,
        default=DEFAULT_SOURCE_ROOT,
        help="Root containing hydrocephalus source trees, e.g. .../ABDN_DATA/Hydrocephalus",
    )
    parser.add_argument(
        "--raw-root",
        type=Path,
        default=DEFAULT_RAW_ROOT,
        help="Destination raw root, e.g. .../bayesian-brain-morphometry/data/raw/hydrocephalus",
    )
    parser.add_argument(
        "--mapping-tsv",
        type=Path,
        default=DEFAULT_MAPPING_TSV,
        help="Where to save the discovered mapping TSV",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Actually copy files instead of dry-run",
    )
    args = parser.parse_args()

    site_roots = discover_site_roots(args.source_root)
    if not site_roots:
        print(f"[WARN] No site roots found under {args.source_root}")
        return 0

    discovered: list[SourceGroup] = []
    for site_root in site_roots:
        site = site_label_from_name(site_root.name)
        group_dirs = discover_group_dirs(site_root)
        for group_dir in group_dirs:
            file_count = len(collect_data_files(group_dir))
            discovered.append(SourceGroup(site=site, group_dir=group_dir, file_count=file_count))

    discovered.sort(key=lambda x: (x.site, str(x.group_dir)))

    print(f"Discovered {len(discovered)} subject groups under {args.source_root}")

    existing_subjects = load_existing_subjects(args.mapping_tsv)
    rows: list[dict[str, str]] = []
    total_copied = 0
    subject_index = 1

    for group in discovered:
        subject_id = f"sub-{subject_index:03d}"
        while subject_id in existing_subjects:
            subject_index += 1
            subject_id = f"sub-{subject_index:03d}"

        rows.append(
            {
                "canonical_subject_id": subject_id,
                "site": group.site,
                "source_group": group.group_dir.name,
                "source_group_path": str(group.group_dir),
                "file_count": str(group.file_count),
            }
        )

        total_copied += copy_group_files(
            subject_id=subject_id,
            group_dir=group.group_dir,
            raw_root=args.raw_root,
            dry_run=not args.apply,
        )

        subject_index += 1

    write_mapping(rows, args.mapping_tsv)
    print(f"Saved mapping to: {args.mapping_tsv}")
    print(f"\nFinished. Matched files: {sum(int(r['file_count']) for r in rows)}. Copied: {total_copied}.")
    if not args.apply:
        print("Dry-run only. Re-run with --apply to copy files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())