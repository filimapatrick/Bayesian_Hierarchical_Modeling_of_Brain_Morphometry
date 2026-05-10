#!/usr/bin/env python3
"""Convert organized dementia DICOM folders to NIfTI and pick the best T1w series.

This script:
- walks each subject's raw DICOM tree,
- runs dcm2niix on each leaf DICOM series folder,
- avoids filename collisions by building a stable prefix from the folder path,
- writes a per-subject conversion log,
- selects the best structural T1-weighted series per subject using JSON metadata,
- writes a summary TSV of the selected T1 series.

Expected input layout:
    data/raw/sub-001/dicom/...

Expected output layout:
    data/derivatives/nifti/sub-001/*.nii.gz
    data/derivatives/nifti/sub-001/*.json
    data/derivatives/nifti/selected_t1w.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RAW_ROOT = REPO_ROOT / "data" / "raw"
DEFAULT_NIFTI_ROOT = REPO_ROOT / "data" / "derivatives" / "nifti"


@dataclass(frozen=True)
class SeriesResult:
    subject_id: str
    series_root: Path
    nifti_stem: str
    json_path: Path | None
    nii_path: Path | None
    score: int = 0
    selected: bool = False


def natural_key(text: str) -> list[object]:
    """Sort strings so run2 comes before run10."""
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", text)]


def find_subjects(raw_root: Path, subject: str | None) -> list[Path]:
    if subject:
        subject_dir = raw_root / subject / "dicom"
        if subject_dir.exists():
            return [subject_dir]
        return []

    subjects = []
    if not raw_root.exists():
        return subjects

    for subj_dir in sorted(raw_root.glob("sub-*")):
        dicom_dir = subj_dir / "dicom"
        if dicom_dir.exists():
            subjects.append(dicom_dir)
    return subjects


def iter_series_roots(dicom_root: Path) -> Iterator[Path]:
    """
    Yield candidate series roots.

    A series root is treated as a folder that contains DICOM files directly.
    We ignore Apple sidecars and .DS_Store.
    """
    if not dicom_root.exists():
        return
    for p in sorted(dicom_root.rglob("*")):
        if not p.is_dir():
            continue
        try:
            children = list(p.iterdir())
        except PermissionError:
            continue
        has_dicom = any(
            c.is_file() and not c.name.startswith("._") and c.name != ".DS_Store"
            for c in children
        )
        if has_dicom:
            yield p


def is_dicom_series_folder(folder: Path) -> bool:
    try:
        return any(
            p.is_file() and not p.name.startswith("._") and p.name != ".DS_Store"
            for p in folder.iterdir()
        )
    except PermissionError:
        return False


def make_prefix(subject_id: str, series_root: Path, dicom_root: Path) -> str:
    """Create a stable, filesystem-safe prefix from the relative path."""
    rel = series_root.relative_to(dicom_root)
    parts = [subject_id] + [part for part in rel.parts if part not in {"", "."}]
    safe = "_".join(parts)
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", safe)
    safe = re.sub(r"_+", "_", safe).strip("_")
    return safe[:180] if len(safe) > 180 else safe


def run_dcm2niix(series_root: Path, output_dir: Path, prefix: str) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "dcm2niix",
        "-z",
        "y",
        "-b",
        "y",
        "-ba",
        "y",
        "-o",
        str(output_dir),
        "-f",
        prefix,
        str(series_root),
    ]

    print(f"[RUN] {series_root} -> {output_dir / prefix}_%s")
    try:
        proc = subprocess.run(cmd, check=False)
    except FileNotFoundError as e:
        raise RuntimeError("dcm2niix was not found on PATH.") from e

    return proc.returncode


def load_json_metadata(json_path: Path) -> dict:
    try:
        with json_path.open() as f:
            return json.load(f)
    except Exception:
        return {}


def score_series(meta: dict) -> int:
    """
    Score a series for structural T1 selection.

    Higher is better.
    """
    desc = f"{meta.get('SeriesDescription', '')} {meta.get('ProtocolName', '')}".lower()

    score = 0

    # Strong positives
    if "mprage" in desc:
        score += 100
    if "t1" in desc:
        score += 30
    if "sag" in desc:
        score += 20

    # Prefer non-contrast if possible
    if "post" in desc or "+c" in desc or "contrast" in desc:
        score -= 10

    # Strong negatives
    if "dwi" in desc or "diff" in desc or "adc" in desc:
        score -= 200
    if "flair" in desc or "t2" in desc:
        score -= 100
    if "tof" in desc:
        score -= 150
    if "cal " in desc or desc.startswith("cal") or "screen save" in desc:
        score -= 200
    if "mip" in desc:
        score -= 150
    if "hemo" in desc or "hemo" in desc:
        score -= 100

    # Small bonus for likely structural sagittal T1 even if not MPRAGE
    if "sag" in desc and "t1" in desc and "mprage" not in desc:
        score += 25

    return score


def locate_converted_pairs(output_dir: Path) -> list[tuple[Path, Path]]:
    """
    Return pairs of (nii_path, json_path) from an output directory.
    """
    pairs: list[tuple[Path, Path]] = []
    if not output_dir.exists():
        return pairs

    for nii in sorted(output_dir.glob("*.nii.gz")):
        stem = nii.name[:-7]  # strip .nii.gz
        js = output_dir / f"{stem}.json"
        if js.exists():
            pairs.append((nii, js))
    return pairs


def pick_best_t1_for_subject(subject_id: str, subject_out: Path) -> SeriesResult | None:
    """
    Select the best structural T1 series for a subject from all converted series.
    """
    candidates: list[SeriesResult] = []

    for nii, js in locate_converted_pairs(subject_out):
        meta = load_json_metadata(js)
        score = score_series(meta)
        candidates.append(
            SeriesResult(
                subject_id=subject_id,
                series_root=subject_out,
                nifti_stem=nii.name[:-7],
                json_path=js,
                nii_path=nii,
                score=score,
                selected=False,
            )
        )

    if not candidates:
        return None

    candidates.sort(key=lambda x: (x.score, x.nifti_stem), reverse=True)
    best = candidates[0]
    return SeriesResult(
        subject_id=best.subject_id,
        series_root=best.series_root,
        nifti_stem=best.nifti_stem,
        json_path=best.json_path,
        nii_path=best.nii_path,
        score=best.score,
        selected=True,
    )


def write_selected_t1_manifest(rows: list[SeriesResult], nifti_root: Path) -> Path:
    out_path = nifti_root / "selected_t1w.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "subject_id",
                "nifti_stem",
                "nii_path",
                "json_path",
                "score",
                "series_description",
                "protocol_name",
            ]
        )
        for r in rows:
            meta = load_json_metadata(r.json_path) if r.json_path else {}
            w.writerow(
                [
                    r.subject_id,
                    r.nifti_stem,
                    str(r.nii_path) if r.nii_path else "",
                    str(r.json_path) if r.json_path else "",
                    r.score,
                    meta.get("SeriesDescription", ""),
                    meta.get("ProtocolName", ""),
                ]
            )
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert dementia DICOM series to NIfTI and select the best T1w scan."
    )
    parser.add_argument(
        "--raw-root",
        type=Path,
        default=DEFAULT_RAW_ROOT,
        help="Root data/raw directory.",
    )
    parser.add_argument(
        "--nifti-root",
        type=Path,
        default=DEFAULT_NIFTI_ROOT,
        help="Root output directory for converted NIfTI files.",
    )
    parser.add_argument(
        "--subject",
        type=str,
        default=None,
        help="Convert only one subject, e.g. sub-001.",
    )
    parser.add_argument(
        "--pick-only",
        action="store_true",
        help="Do not run dcm2niix; only score existing NIfTI JSON pairs and write selected_t1w.tsv.",
    )
    args = parser.parse_args()

    subjects = find_subjects(args.raw_root, args.subject)
    if not subjects:
        if args.subject:
            print(f"[WARN] No DICOM roots found under: {args.raw_root / args.subject / 'dicom'}")
        else:
            print(f"[WARN] No DICOM roots found under: {args.raw_root}")
        return 0

    all_selected: list[SeriesResult] = []
    planned_or_run = 0
    ok = 0
    failed = 0

    for dicom_root in subjects:
        subject_id = dicom_root.parent.name
        subject_out = args.nifti_root / subject_id
        subject_out.mkdir(parents=True, exist_ok=True)

        if not args.pick_only:
            series_roots = [p for p in iter_series_roots(dicom_root) if is_dicom_series_folder(p)]
            series_roots = sorted(series_roots, key=lambda p: natural_key(str(p.relative_to(dicom_root))))

            if not series_roots:
                print(f"[WARN] No series folders found under {dicom_root}")
            for idx, series_root in enumerate(series_roots, start=1):
                prefix = make_prefix(subject_id, series_root, dicom_root)
                rc = run_dcm2niix(series_root, subject_out, prefix)
                planned_or_run += 1
                if rc == 0:
                    ok += 1
                else:
                    failed += 1
                    print(f"[WARN] dcm2niix failed for {series_root} with exit code {rc}")

        best = pick_best_t1_for_subject(subject_id, subject_out)
        if best is not None:
            all_selected.append(best)
            meta = load_json_metadata(best.json_path) if best.json_path else {}
            print(
                f"[SELECT] {subject_id}: {best.nifti_stem} "
                f"(score={best.score}, SeriesDescription={meta.get('SeriesDescription', '')!r})"
            )
        else:
            print(f"[WARN] {subject_id}: no converted NIfTI/JSON pairs found for T1 selection")

    manifest_path = write_selected_t1_manifest(all_selected, args.nifti_root)
    print(f"[WRITE] {manifest_path}")

    if args.pick_only:
        print(f"Done. Selected {len(all_selected)} subject T1 series.")
    else:
        print(f"Done. Planned/ran {planned_or_run} conversion jobs. OK={ok}, failed={failed}.")
        print(f"Selected T1 series for {len(all_selected)} subject(s).")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())