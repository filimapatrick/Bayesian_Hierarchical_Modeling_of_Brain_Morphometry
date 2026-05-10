#!/usr/bin/env python3
"""Organize dementia DICOM exports into a canonical subject-based layout."""

from __future__ import annotations

import argparse
import csv
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_ROOT = REPO_ROOT / "data"
DEFAULT_RAW_ROOT = DEFAULT_DATA_ROOT / "raw"
DEFAULT_MAPPING_CSV = Path("/Volumes/MyHDD/ABDN_DATA/dementia_subject_mapping.csv")
DEFAULT_SOURCES_ROOT = Path("/Volumes/MyHDD/ABDN_DATA")


@dataclass(frozen=True)
class SubjectRow:
    canonical_subject_id: str
    source: str
    batch: str
    export_tag: str
    diagnosis: str
    subject_group: str
    file_count: int


@dataclass(frozen=True)
class SubjectMatch:
    matched_root: Path
    files: list[Path]


def load_mapping(mapping_csv: Path) -> list[SubjectRow]:
    rows: list[SubjectRow] = []

    with mapping_csv.open(newline="") as f:
        reader = csv.DictReader(f)

        required = {
            "canonical_subject_id",
            "source",
            "batch",
            "export_tag",
            "diagnosis",
            "subject_group",
            "file_count",
        }

        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Mapping CSV is missing columns: {sorted(missing)}")

        for row in reader:
            rows.append(
                SubjectRow(
                    canonical_subject_id=row["canonical_subject_id"].strip(),
                    source=row["source"].strip(),
                    batch=row["batch"].strip(),
                    export_tag=row["export_tag"].strip(),
                    diagnosis=row["diagnosis"].strip(),
                    subject_group=row["subject_group"].strip(),
                    file_count=int(row["file_count"]),
                )
            )

    return rows


def source_candidates(row: SubjectRow, sources_root: Path) -> list[Path]:
    candidates: list[Path] = []

    if row.source.upper() == "UPTH":
        if row.batch:
            candidates.append(
                sources_root
                / "Dementia 2"
                / "UPTH"
                / row.batch
                / "UPTH"
                / "SL"
                / row.diagnosis
                / row.subject_group
            )
        candidates.append(
            sources_root / "Dementia 2" / "UPTH" / "SL" / row.diagnosis / row.subject_group
        )

    elif row.source.upper() == "BMH":
        candidates.append(
            sources_root
            / "Dementia 2"
            / "bmh"
            / row.batch
            / row.source
            / row.export_tag
            / row.diagnosis
            / row.subject_group
        )
        candidates.append(
            sources_root
            / "Dementia 2"
            / "bmh"
            / row.batch
            / row.source
            / row.export_tag
            / row.diagnosis
        )

    return candidates


def iter_dicom_files(root: Path) -> Iterable[Path]:
    """Yield real files only (skip macOS sidecars)."""
    if not root.exists():
        return []

    for p in root.rglob("*"):
        if not p.is_file():
            continue
        if p.name == ".DS_Store" or p.name.startswith("._"):
            continue
        yield p


def collect_subject_files(row: SubjectRow, sources_root: Path) -> SubjectMatch | None:
    for candidate in source_candidates(row, sources_root):
        files = list(iter_dicom_files(candidate))
        if files:
            return SubjectMatch(matched_root=candidate, files=files)
    return None


def copy_subject_files(
    subject_id: str,
    matched_root: Path,
    files: list[Path],
    raw_root: Path,
    dry_run: bool = True,
) -> int:
    target_dir = raw_root / subject_id / "dicom"

    if dry_run:
        print(f"[DRY-RUN] {subject_id}: would copy {len(files)} files -> {target_dir}")
        return len(files)

    target_dir.mkdir(parents=True, exist_ok=True)

    copied = 0
    for src in files:
        if src.name == ".DS_Store" or src.name.startswith("._"):
            continue

        rel = src.relative_to(matched_root)
        dst = target_dir / rel
        dst.parent.mkdir(parents=True, exist_ok=True)

        if dst.exists():
            continue

        shutil.copyfile(src, dst)
        copied += 1

    return copied


def remove_sidecar_files(root: Path) -> int:
    """Delete macOS sidecar files from a copied tree."""
    removed = 0
    if not root.exists():
        return removed

    for p in root.rglob("*"):
        if p.is_file() and (p.name.startswith("._") or p.name == ".DS_Store"):
            p.unlink()
            removed += 1

    return removed


def main() -> int:
    parser = argparse.ArgumentParser(description="Organize dementia DICOM data")
    parser.add_argument("--mapping-csv", type=Path, default=DEFAULT_MAPPING_CSV)
    parser.add_argument("--sources-root", type=Path, default=DEFAULT_SOURCES_ROOT)
    parser.add_argument("--raw-root", type=Path, default=DEFAULT_RAW_ROOT)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    if not args.mapping_csv.exists():
        print(f"Mapping CSV not found: {args.mapping_csv}")
        return 2

    rows = load_mapping(args.mapping_csv)
    print(f"Loaded {len(rows)} subject rows")

    total_found = 0
    total_copied = 0
    total_cleaned = 0

    for row in rows:
        match = collect_subject_files(row, args.sources_root)

        if match is None:
            print(f"[WARN] {row.canonical_subject_id}: no files found")
            continue

        total_found += len(match.files)

        copied = copy_subject_files(
            row.canonical_subject_id,
            match.matched_root,
            match.files,
            args.raw_root,
            dry_run=not args.apply,
        )
        total_copied += copied

        if args.apply:
            cleaned = remove_sidecar_files(args.raw_root / row.canonical_subject_id / "dicom")
            total_cleaned += cleaned
            if cleaned:
                print(f"[CLEAN] {row.canonical_subject_id}: removed {cleaned} sidecar files")

    print(f"Finished. Matched: {total_found}, Copied: {total_copied}, Cleaned: {total_cleaned}")

    if not args.apply:
        print("Dry-run only. Use --apply to execute.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())