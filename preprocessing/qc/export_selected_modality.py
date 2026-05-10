#!/usr/bin/env python3
"""Export one selected modality per subject into a clean analysis-ready layout.

This script is meant to be the next step after DICOM organization and DICOM -> NIfTI
conversion. It copies the chosen image for each subject into a canonical folder such as:

    data/derivatives/analysis_ready/sub-001/T1w.nii.gz
    data/derivatives/analysis_ready/sub-001/T1w.json

It is intentionally flexible so it can work with slightly different manifest formats.

Supported manifest columns (case-insensitive):
- subject_id / canonical_subject_id / subject
- nifti_path / nii_path / image_path / source_nifti
- json_path / sidecar_path / metadata_path / source_json
- stem / series_id / nifti_stem (optional fallback)

If nifti_path is not provided, the script will try to resolve a matching NIfTI inside
    <nifti_root>/<subject_id>/
by using the stem or by searching for a single .nii.gz file.

Usage examples:

    python export_selected_modality.py \
        --selection-manifest /path/to/selected_t1w.tsv \
        --nifti-root /path/to/data/derivatives/nifti \
        --output-root /path/to/data/derivatives/analysis_ready \
        --modality-name T1w --apply

    python export_selected_modality.py --selection-manifest selected_t1w.tsv

Recommended next step:
- Generate a selection manifest with one row per subject containing the chosen T1w.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional


DEFAULT_MANIFEST = Path("/Volumes/MyHDD/ABDN_DATA/selected_t1w.tsv")
DEFAULT_NIFTI_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry/data/derivatives/nifti")
DEFAULT_OUTPUT_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry/data/derivatives/analysis_ready")
DEFAULT_MODALITY_NAME = "T1w"


@dataclass(frozen=True)
class SelectionRow:
    subject_id: str
    nifti_path: Optional[Path]
    json_path: Optional[Path]
    stem: Optional[str]
    raw: dict[str, str]


@dataclass(frozen=True)
class ExportResult:
    subject_id: str
    status: str
    nifti_source: str
    json_source: str
    nifti_target: str
    json_target: str
    message: str


def _normalize_key(key: str) -> str:
    return key.strip().lower().replace(" ", "_")


def _first_existing_path(candidates: Iterable[Path]) -> Optional[Path]:
    for candidate in candidates:
        if candidate and candidate.exists():
            return candidate
    return None


def read_manifest(manifest_path: Path) -> list[SelectionRow]:
    if not manifest_path.exists():
        raise FileNotFoundError(f"Selection manifest not found: {manifest_path}")

    suffix = manifest_path.suffix.lower()
    delimiter = "\t" if suffix in {".tsv", ".tab"} else ","

    rows: list[SelectionRow] = []
    with manifest_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Manifest has no header: {manifest_path}")

        normalized_fields = {_normalize_key(name): name for name in reader.fieldnames}

        def pick(row: dict[str, str], *names: str) -> str:
            for name in names:
                key = normalized_fields.get(_normalize_key(name))
                if key and row.get(key) is not None:
                    value = row[key].strip()
                    if value:
                        return value
            return ""

        for row in reader:
            subject_id = pick(row, "subject_id", "canonical_subject_id", "subject")
            if not subject_id:
                continue

            nifti_text = pick(row, "nifti_path", "nii_path", "image_path", "source_nifti")
            json_text = pick(row, "json_path", "sidecar_path", "metadata_path", "source_json")
            stem = pick(row, "stem", "series_id", "nifti_stem")

            rows.append(
                SelectionRow(
                    subject_id=subject_id,
                    nifti_path=Path(nifti_text) if nifti_text else None,
                    json_path=Path(json_text) if json_text else None,
                    stem=stem or None,
                    raw={k: (v or "").strip() for k, v in row.items()},
                )
            )
    return rows


def subject_root(nifti_root: Path, subject_id: str) -> Path:
    return nifti_root / subject_id


def infer_nifti_and_json(row: SelectionRow, nifti_root: Path) -> tuple[Optional[Path], Optional[Path], str]:
    """Resolve the source NIfTI and JSON for one subject.

    Returns:
        (nifti_path, json_path, message)
    """
    # 1) Explicit paths in manifest.
    if row.nifti_path and row.nifti_path.exists():
        nifti_path = row.nifti_path
        json_path = row.json_path if row.json_path and row.json_path.exists() else nifti_path.with_suffix(".json")
        if not json_path.exists():
            return nifti_path, None, "JSON sidecar missing"
        return nifti_path, json_path, "Resolved from manifest paths"

    # 2) Search inside the subject folder.
    sroot = subject_root(nifti_root, row.subject_id)
    if not sroot.exists():
        return None, None, f"Subject folder not found: {sroot}"

    # A stem may be given, e.g. sub-001_series_10 or sub-015_DH_01_DG_01_1_1.
    if row.stem:
        stem = row.stem
        candidates = [
            sroot / f"{stem}.nii.gz",
            sroot / f"{stem}.nii",
            sroot / f"{stem}.json",
        ]
        nifti_path = _first_existing_path(candidates[:2])
        if nifti_path:
            json_path = sroot / f"{nifti_path.stem}.json"
            if json_path.exists():
                return nifti_path, json_path, f"Resolved from stem={stem}"
            return nifti_path, None, f"NIfTI found for stem={stem}, but JSON missing"

    # 3) Single obvious NIfTI.
    nifti_candidates = sorted(sroot.rglob("*.nii.gz")) + sorted(sroot.rglob("*.nii"))
    nifti_candidates = [p for p in nifti_candidates if p.is_file()]

    if len(nifti_candidates) == 1:
        nifti_path = nifti_candidates[0]
        json_path = nifti_path.with_suffix(".json")
        if json_path.exists():
            return nifti_path, json_path, "Resolved by single-file fallback"
        return nifti_path, None, "Single NIfTI found, but JSON missing"

    # 4) If multiple candidates, prefer files that look like T1w.
    t1_candidates = [
        p for p in nifti_candidates
        if "t1" in p.name.lower() or "mprage" in p.name.lower() or "t1w" in p.name.lower()
    ]
    if len(t1_candidates) == 1:
        nifti_path = t1_candidates[0]
        json_path = nifti_path.with_suffix(".json")
        if json_path.exists():
            return nifti_path, json_path, "Resolved by T1-like filename fallback"
        return nifti_path, None, "T1-like NIfTI found, but JSON missing"

    return None, None, f"Ambiguous NIfTI selection under {sroot}"


def copy_pair(
    subject_id: str,
    nifti_path: Path,
    json_path: Path,
    output_root: Path,
    modality_name: str,
    overwrite: bool,
    dry_run: bool,
) -> ExportResult:
    out_dir = output_root / subject_id
    out_dir.mkdir(parents=True, exist_ok=True)

    target_nii = out_dir / f"{modality_name}.nii.gz"
    target_json = out_dir / f"{modality_name}.json"

    if dry_run:
        return ExportResult(
            subject_id=subject_id,
            status="dry_run",
            nifti_source=str(nifti_path),
            json_source=str(json_path),
            nifti_target=str(target_nii),
            json_target=str(target_json),
            message="Would copy",
        )

    if target_nii.exists() and not overwrite:
        return ExportResult(
            subject_id=subject_id,
            status="skipped_exists",
            nifti_source=str(nifti_path),
            json_source=str(json_path),
            nifti_target=str(target_nii),
            json_target=str(target_json),
            message="Target already exists; use --overwrite to replace",
        )

    shutil.copy2(nifti_path, target_nii)
    shutil.copy2(json_path, target_json)

    return ExportResult(
        subject_id=subject_id,
        status="copied",
        nifti_source=str(nifti_path),
        json_source=str(json_path),
        nifti_target=str(target_nii),
        json_target=str(target_json),
        message="Copied",
    )


def write_manifest(results: list[ExportResult], output_root: Path) -> Path:
    manifest_path = output_root / "selected_modality_export_manifest.csv"
    output_root.mkdir(parents=True, exist_ok=True)

    with manifest_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "subject_id",
            "status",
            "nifti_source",
            "json_source",
            "nifti_target",
            "json_target",
            "message",
        ])
        for r in results:
            writer.writerow([
                r.subject_id,
                r.status,
                r.nifti_source,
                r.json_source,
                r.nifti_target,
                r.json_target,
                r.message,
            ])
    return manifest_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Export one selected modality per subject into analysis-ready folders.")
    parser.add_argument("--selection-manifest", type=Path, default=DEFAULT_MANIFEST, help="TSV/CSV manifest with selected modality rows")
    parser.add_argument("--nifti-root", type=Path, default=DEFAULT_NIFTI_ROOT, help="Root directory containing subject NIfTI outputs")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT, help="Destination root for analysis-ready outputs")
    parser.add_argument("--modality-name", default=DEFAULT_MODALITY_NAME, help="Name used for exported files, e.g. T1w")
    parser.add_argument("--apply", action="store_true", help="Actually copy files; otherwise dry-run")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing exported files")
    args = parser.parse_args()

    rows = read_manifest(args.selection_manifest)
    print(f"Loaded {len(rows)} selection row(s) from {args.selection_manifest}")

    results: list[ExportResult] = []
    copied = 0
    skipped = 0
    failed = 0

    for row in rows:
        nifti_path, json_path, message = infer_nifti_and_json(row, args.nifti_root)
        if nifti_path is None or json_path is None:
            print(f"[WARN] {row.subject_id}: {message}")
            results.append(
                ExportResult(
                    subject_id=row.subject_id,
                    status="failed",
                    nifti_source=str(nifti_path or ""),
                    json_source=str(json_path or ""),
                    nifti_target=str(args.output_root / row.subject_id / f"{args.modality_name}.nii.gz"),
                    json_target=str(args.output_root / row.subject_id / f"{args.modality_name}.json"),
                    message=message,
                )
            )
            failed += 1
            continue

        result = copy_pair(
            subject_id=row.subject_id,
            nifti_path=nifti_path,
            json_path=json_path,
            output_root=args.output_root,
            modality_name=args.modality_name,
            overwrite=args.overwrite,
            dry_run=not args.apply,
        )
        results.append(result)

        if result.status in {"copied", "dry_run"}:
            copied += 1
            print(f"[{result.status.upper()}] {row.subject_id}: {message} -> {result.nifti_target}")
        elif result.status == "skipped_exists":
            skipped += 1
            print(f"[SKIP] {row.subject_id}: {result.message}")
        else:
            failed += 1
            print(f"[WARN] {row.subject_id}: {result.message}")

    manifest_out = write_manifest(results, args.output_root)
    print(
        f"Finished. Subjects processed: {len(rows)}. Exported: {copied}. Skipped: {skipped}. Failed: {failed}."
    )
    print(f"Wrote export manifest: {manifest_out}")

    if not args.apply:
        print("Dry-run only. Re-run with --apply to copy files.")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
