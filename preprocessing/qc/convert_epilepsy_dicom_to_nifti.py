#!/usr/bin/env python3
"""
Convert organized DICOM folders to NIfTI using dcm2niix.

This version works with both layouts:
- flat folders:
    data/raw/epilepsy/sub-001/dicom/*.dcm
- nested series folders:
    data/raw/dementia/sub-013/dicom/DF_02/*.dcm

It scans each subject folder under --raw-root, finds conversion roots, and
runs dcm2niix for each one.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry")
DEFAULT_RAW_ROOT = REPO_ROOT / "data" / "raw"
DEFAULT_NIFTI_ROOT = REPO_ROOT / "data" / "derivatives" / "nifti"


def is_data_file(path: Path) -> bool:
    return path.is_file() and not path.name.startswith("._") and path.name != ".DS_Store"


def dir_has_data_files(root: Path) -> bool:
    if not root.exists():
        return False
    return any(is_data_file(p) for p in root.rglob("*"))


def discover_subject_dirs(raw_root: Path) -> list[Path]:
    if not raw_root.exists():
        return []
    return [p for p in sorted(raw_root.iterdir()) if p.is_dir() and not p.name.startswith(".")]


def find_conversion_roots(subject_dicom_root: Path) -> list[Path]:
    """
    Return the folders that should be passed to dcm2niix.

    Priority:
    1. Immediate child directories that contain data files
    2. The dicom folder itself if it contains data files directly
    """
    if not subject_dicom_root.exists():
        return []

    child_dirs = [
        p for p in sorted(subject_dicom_root.iterdir())
        if p.is_dir() and dir_has_data_files(p)
    ]
    if child_dirs:
        return child_dirs

    if dir_has_data_files(subject_dicom_root):
        return [subject_dicom_root]

    return []


def run_dcm2niix(input_root: Path, output_dir: Path, filename_prefix: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "dcm2niix",
        "-z", "y",
        "-b", "y",
        "-ba", "y",
        "-o", str(output_dir),
        "-f", filename_prefix,
        str(input_root),
    ]

    print(f"[RUN] {input_root} -> {output_dir}/{filename_prefix}")
    subprocess.run(cmd, check=True)


def process_subject(subject_dir: Path, nifti_root: Path) -> tuple[int, int]:
    subject_id = subject_dir.name
    dicom_root = subject_dir / "dicom"

    conversion_roots = find_conversion_roots(dicom_root)
    if not conversion_roots:
        print(f"[WARN] No series folders or DICOM files found under {dicom_root}")
        return 0, 0

    subject_out = nifti_root / subject_id
    ok = 0
    failed = 0

    for idx, root in enumerate(conversion_roots, start=1):
        prefix = f"{subject_id}_run{idx:02d}"
        try:
            run_dcm2niix(root, subject_out, prefix)
            ok += 1
        except subprocess.CalledProcessError as e:
            failed += 1
            print(f"[WARN] dcm2niix failed for {root} with exit code {e.returncode}")

    return ok, failed


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert DICOM folders to NIfTI using dcm2niix.")
    parser.add_argument(
        "--raw-root",
        type=Path,
        default=DEFAULT_RAW_ROOT,
        help="Path to the pathology root, e.g. .../data/raw/dementia or .../data/raw/epilepsy",
    )
    parser.add_argument(
        "--nifti-root",
        type=Path,
        default=DEFAULT_NIFTI_ROOT,
        help="Output root for converted NIfTI files, e.g. .../data/derivatives/nifti/dementia",
    )
    parser.add_argument(
        "--subject",
        nargs="*",
        default=None,
        help="Optional subject IDs to process, e.g. sub-001 sub-002",
    )
    args = parser.parse_args()

    subject_dirs = discover_subject_dirs(args.raw_root)
    if args.subject:
        wanted = set(args.subject)
        subject_dirs = [d for d in subject_dirs if d.name in wanted]

    if not subject_dirs:
        print(f"[WARN] No subject folders found under {args.raw_root}")
        return 0

    total_ok = 0
    total_failed = 0

    for subject_dir in subject_dirs:
        ok, failed = process_subject(subject_dir, args.nifti_root)
        total_ok += ok
        total_failed += failed

    print(f"\nDone. Planned/ran {total_ok + total_failed} conversion jobs. OK={total_ok}, failed={total_failed}.")
    return 0 if total_failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())