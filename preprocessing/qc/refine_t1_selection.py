#!/usr/bin/env python3
"""
Refine T1 selection from dcm2niix outputs.

This script scans each subject folder in:
    /Volumes/MyHDD/bayesian-brain-morphometry/data/derivatives/nifti/sub-XXX

and chooses the best structural T1 candidate based on metadata keywords.

It writes:
    /Volumes/MyHDD/ABDN_DATA/selected_t1w.tsv

Output columns:
    subject_id    nifti_path    json_path    series_description    protocol_name    score

The goal is to prefer:
- MPRAGE / MP-RAGE
- SPGR / BRAVO
- T1 SE / T1 TSE structural scans

and avoid:
- localizers / scouts
- MIP / projections
- TOF / angiography
- DWI / ADC / trace
- screen saves / cal scans / derived images
- obvious non-structural series
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

ANALYSIS_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry/data/derivatives/nifti")
DEFAULT_OUT = Path("/Volumes/MyHDD/ABDN_DATA/selected_t1w.tsv")


NEGATIVE_PATTERNS = [
    "localizer",
    "scout",
    "survey",
    "mip",
    "projection",
    "screen save",
    "cal ",
    "calibration",
    "adc",
    "dwi",
    "diff",
    "tracew",
    "trace",
    "tof",
    "angio",
    "angiography",
    "venous",
    "venogram",
    "mra",
    "screenshot",
    "derived",
    "secondary",
    "reformat",
    "mpr ",
    "mpr_",
    "mip_sag",
    "mip_cor",
    "mip_tra",
]

POSITIVE_3D_T1 = [
    "mprage",
    "mp-rage",
    "bravo",
    "spgr",
    "t1_mprage",
    "t1_mpr",
    "3d t1",
    "3d-t1",
    "3d_t1",
    "t1 cube",
    "cube t1",
    "tfe",
    "gre t1",
]

POSITIVE_T1 = [
    "t1",
    "t1w",
    "t1 weighted",
    "t1-weighted",
    "t1 se",
    "t1 tse",
    "t1_fse",
    "t1 se ",
    "t1 se_",
    "t1 se-",
]

CONTRAST_HINTS = [
    "post contrast",
    "post-contrast",
    "postcontrast",
    "pc",
    "+c",
    " + c",
    "c+",
    "contrast",
]

DERIVED_HINTS = [
    "mip",
    "screen save",
    "cal ",
    "calibration",
    "localizer",
    "scout",
]


@dataclass(frozen=True)
class Candidate:
    subject_id: str
    nifti_path: Path
    json_path: Path
    series_description: str
    protocol_name: str
    score: int


def norm_text(value: object) -> str:
    return str(value or "").strip().lower()


def candidate_jsons(subject_dir: Path) -> Iterable[Path]:
    for js in sorted(subject_dir.glob("*.json")):
        yield js


def load_metadata(js_path: Path) -> dict:
    with js_path.open() as f:
        return json.load(f)


def score_series(meta: dict, nifti_path: Path, json_path: Path) -> int:
    """
    Heuristic scoring:
    - strong positive for structural T1 keywords
    - penalty for non-structural / derived series
    - mild penalty for post-contrast when a cleaner pre-contrast T1 exists
    """
    desc = norm_text(meta.get("SeriesDescription"))
    prot = norm_text(meta.get("ProtocolName"))
    image_type = " ".join(
        norm_text(x) for x in (
            meta.get("ImageType"),
            meta.get("ScanningSequence"),
            meta.get("SequenceName"),
        )
    )
    name = " ".join([desc, prot, image_type])

    score = 0

    # Strong evidence of a structural T1.
    if any(k in name for k in POSITIVE_3D_T1):
        score += 120
    if any(k in name for k in POSITIVE_T1):
        score += 60

    # Bonus for common structural T1 markers.
    if "t1" in name:
        score += 20
    if "mprage" in name or "spgr" in name or "bravo" in name:
        score += 40

    # Mild bonus if it looks like a full structural acquisition rather than an image-derived series.
    if "brain" in name:
        score += 5

    # Penalties for obviously wrong series.
    if any(k in name for k in NEGATIVE_PATTERNS):
        score -= 120

    # Do not pick obvious MIPs or projection images.
    if "mip" in desc or "mip" in prot:
        score -= 80

    # Prefer pre-contrast over post-contrast, but still keep contrast T1 if nothing better exists.
    if any(k in name for k in CONTRAST_HINTS):
        score -= 15

    # If it says trace/diffusion/adc, it is not your T1.
    if "adc" in name or "trace" in name or "dwi" in name or "diffusion" in name:
        score -= 150

    # If it says TOF or angio, not T1.
    if "tof" in name or "angi" in name:
        score -= 150

    # Prefer series with "t1" over series with only a generic structural name.
    if "t1" in name:
        score += 10

    # Slight bonus for likely 3D structural language.
    if "3d" in name:
        score += 10

    return score


def best_candidate_for_subject(subject_dir: Path) -> Optional[Candidate]:
    subject_id = subject_dir.name
    best: Optional[Candidate] = None

    for js in candidate_jsons(subject_dir):
        nii = js.with_suffix(".nii.gz")
        if not nii.exists():
            continue

        try:
            meta = load_metadata(js)
        except Exception as exc:
            print(f"[WARN] {subject_id}: could not read {js.name}: {exc}")
            continue

        score = score_series(meta, nii, js)
        desc = str(meta.get("SeriesDescription") or "")
        prot = str(meta.get("ProtocolName") or "")

        # Require at least some T1 evidence.
        name = f"{desc} {prot}".lower()
        if "t1" not in name and not any(k in name for k in POSITIVE_3D_T1):
            continue

        cand = Candidate(
            subject_id=subject_id,
            nifti_path=nii,
            json_path=js,
            series_description=desc,
            protocol_name=prot,
            score=score,
        )

        if best is None or cand.score > best.score:
            best = cand

    return best


def write_manifest(rows: list[Candidate], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "subject_id",
            "nifti_path",
            "json_path",
            "series_description",
            "protocol_name",
            "score",
        ])
        for r in rows:
            w.writerow([
                r.subject_id,
                str(r.nifti_path),
                str(r.json_path),
                r.series_description,
                r.protocol_name,
                r.score,
            ])


def main() -> int:
    parser = argparse.ArgumentParser(description="Refine T1 selection from dcm2niix outputs.")
    parser.add_argument(
        "--analysis-root",
        type=Path,
        default=ANALYSIS_ROOT,
        help="Folder containing sub-*/ with NIfTI/JSON outputs.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT,
        help="Output TSV manifest.",
    )
    parser.add_argument(
        "--subject",
        nargs="*",
        default=None,
        help="Optional subject IDs to process, e.g. sub-001 sub-002.",
    )
    parser.add_argument(
        "--preview",
        action="store_true",
        help="Print ranked candidates per subject.",
    )
    args = parser.parse_args()

    subjects = sorted(args.analysis_root.glob("sub-*"))
    if args.subject:
        wanted = set(args.subject)
        subjects = [p for p in subjects if p.name in wanted]

    selected: list[Candidate] = []

    for subj_dir in subjects:
        best = best_candidate_for_subject(subj_dir)
        if best is None:
            print(f"[WARN] {subj_dir.name}: no plausible T1 candidate found")
            continue

        selected.append(best)
        print(
            f"[OK] {subj_dir.name}: selected {best.json_path.name} "
            f"(score={best.score}, desc={best.series_description!r})"
        )

        if args.preview:
            print(f"      nifti: {best.nifti_path}")
            print(f"      json : {best.json_path}")

    write_manifest(selected, args.out)
    print(f"\nSaved {len(selected)} T1 selections to: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())