#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional


REPO_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry")
DEFAULT_ANALYSIS_ROOT = REPO_ROOT / "data" / "derivatives" / "nifti" / "control"
DEFAULT_OUT = Path("/Volumes/MyHDD/ABDN_DATA/selected_t1w_control.tsv")


BAD_PATTERNS = [
    "flair",
    "t2",
    "t2*",
    "t2star",
    "dwi",
    "diff",
    "adc",
    "tof",
    "mip",
    "cal",
    "screen save",
    "localizer",
    "scout",
    "survey",
    "sag t2",
    "cor t2",
    "ax t2",
]

GOOD_PATTERNS = [
    "t1",
    "mprage",
    "mpr",
    "spgr",
    "cube",
    "fse",
    "se",
    "post contrast",
    "postcontrast",
    "+c",
    "+ c",
    "contrast",
    "pc",
    "optimized",
]


@dataclass(frozen=True)
class Candidate:
    subject_id: str
    nifti_path: Path
    json_path: Path
    series_description: str
    protocol_name: str
    score: int


def normalize_text(value: str) -> str:
    return re.sub(r"\s+", " ", (value or "").strip().lower())


def score_series(series_description: str, protocol_name: str, nifti_path: Path) -> Optional[int]:
    text = normalize_text(series_description + " " + protocol_name + " " + nifti_path.name)

    if not text:
        return None

    # Hard exclusions for obvious non-T1s or derived images.
    for bad in BAD_PATTERNS:
        if bad in text:
            return None

    score = 0

    # Strong T1 signals.
    if "t1" in text:
        score += 60
    if "mprage" in text:
        score += 120
    if "mpr" in text:
        score += 40
    if "spgr" in text:
        score += 50
    if "cube" in text:
        score += 40
    if "fse" in text:
        score += 20
    if "se" in text:
        score += 10

    # Orientation cues that often correspond to structural T1s.
    if "sag" in text:
        score += 15
    if "ax" in text:
        score += 10
    if "cor" in text:
        score += 10

    # Contrast / post-contrast structural T1s are still valid anatomical T1s.
    if "post contrast" in text or "postcontrast" in text:
        score += 35
    if "+c" in text or "+ c" in text:
        score += 30
    if "contrast" in text:
        score += 25
    if "pc" in text:
        score += 20

    if "optimized" in text:
        score += 10

    # Some helpful penalties.
    if "localizer" in text or "scout" in text or "survey" in text:
        score -= 100

    # Need at least some T1 signal.
    if "t1" not in text and "mprage" not in text and "spgr" not in text:
        return None

    return score


def iter_subject_dirs(root: Path) -> list[Path]:
    if not root.exists():
        return []
    return [p for p in sorted(root.iterdir()) if p.is_dir() and p.name.startswith("sub-")]


def find_json_files(subject_dir: Path) -> list[Path]:
    return sorted(subject_dir.rglob("*.json"))


def read_candidate(json_path: Path) -> Optional[Candidate]:
    try:
        with json_path.open() as f:
            meta = json.load(f)
    except Exception:
        return None

    series_description = str(meta.get("SeriesDescription", "") or "")
    protocol_name = str(meta.get("ProtocolName", "") or "")
    nifti_path = json_path.with_suffix(".nii.gz")

    if not nifti_path.exists():
        return None

    score = score_series(series_description, protocol_name, nifti_path)
    if score is None:
        return None

    subject_id = json_path.parents[1].name if json_path.parent.name != json_path.parents[0].name else json_path.parents[0].name
    # More robust fallback: use the nearest sub-* folder in the path.
    for part in json_path.parts:
        if part.startswith("sub-"):
            subject_id = part
            break

    return Candidate(
        subject_id=subject_id,
        nifti_path=nifti_path,
        json_path=json_path,
        series_description=series_description,
        protocol_name=protocol_name,
        score=score,
    )


def choose_best(candidates: list[Candidate]) -> Optional[Candidate]:
    if not candidates:
        return None
    # Highest score wins, then stable path order.
    candidates = sorted(
        candidates,
        key=lambda c: (-c.score, c.json_path.name.lower(), str(c.json_path)),
    )
    return candidates[0]


def main() -> int:
    parser = argparse.ArgumentParser(description="Select the best T1 series for each control subject.")
    parser.add_argument(
        "--analysis-root",
        type=Path,
        default=DEFAULT_ANALYSIS_ROOT,
        help="Root of converted control NIfTI files",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT,
        help="Output TSV for selected T1s",
    )
    args = parser.parse_args()

    rows: list[dict[str, str]] = []
    subject_dirs = iter_subject_dirs(args.analysis_root)

    if not subject_dirs:
        print(f"[WARN] No subjects found in {args.analysis_root}")
        return 0

    for subject_dir in subject_dirs:
        subject_id = subject_dir.name
        candidates: list[Candidate] = []

        for json_path in find_json_files(subject_dir):
            cand = read_candidate(json_path)
            if cand is not None:
                candidates.append(cand)

        best = choose_best(candidates)
        if best is None:
            print(f"[WARN] {subject_id}: no plausible T1 candidate found")
            continue

        rows.append(
            {
                "subject_id": best.subject_id,
                "nifti_path": str(best.nifti_path),
                "json_path": str(best.json_path),
                "series_description": best.series_description,
                "protocol_name": best.protocol_name,
                "score": str(best.score),
            }
        )
        print(
            f"[OK] {subject_id}: selected {best.json_path.name} "
            f"(score={best.score}, desc='{best.series_description}')"
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "subject_id",
                "nifti_path",
                "json_path",
                "series_description",
                "protocol_name",
                "score",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nSaved {len(rows)} T1 selections to: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())