#!/usr/bin/env python3
"""
Extract simple T1-derived features from preprocessed subjects.

Input folder structure:
    data/derivatives/preprocessed/sub-XXX/

Expected files per subject:
    - T1w.nii.gz
    - T1w_brain.nii.gz
    - T1w_brain_mask.nii.gz
    - T1w_brain_fast_pve_0.nii.gz
    - T1w_brain_fast_pve_1.nii.gz
    - T1w_brain_fast_pve_2.nii.gz

Outputs:
    data/derivatives/features/t1_features.csv

Features computed:
    - voxel_volume_mm3
    - brain_volume_ml
    - csf_volume_ml
    - gm_volume_ml
    - wm_volume_ml
    - tissue_sum_volume_ml
    - tissue_balance_ml
    - brain_mean_intensity
    - brain_std_intensity
    - brain_median_intensity
    - brain_min_intensity
    - brain_max_intensity

This is intentionally simple and robust, so it can be reused later for other pathologies.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import nibabel as nib
import numpy as np


REPO_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry")
DEFAULT_PREPROCESSED_ROOT = REPO_ROOT / "data" / "derivatives" / "preprocessed"
DEFAULT_OUT_CSV = REPO_ROOT / "data" / "derivatives" / "features" / "t1_features.csv"


@dataclass(frozen=True)
class SubjectFeatures:
    subject_id: str
    voxel_volume_mm3: float
    brain_volume_ml: float
    csf_volume_ml: float
    gm_volume_ml: float
    wm_volume_ml: float
    tissue_sum_volume_ml: float
    tissue_balance_ml: float
    brain_mean_intensity: float
    brain_std_intensity: float
    brain_median_intensity: float
    brain_min_intensity: float
    brain_max_intensity: float


def load_nii(path: Path) -> np.ndarray:
    img = nib.load(str(path))
    return img.get_fdata(dtype=np.float32), img.affine, img.header


def voxel_volume_mm3(header) -> float:
    zooms = header.get_zooms()[:3]
    return float(np.prod(zooms))


def safe_mean(x: np.ndarray) -> float:
    return float(np.mean(x)) if x.size else float("nan")


def safe_std(x: np.ndarray) -> float:
    return float(np.std(x)) if x.size else float("nan")


def safe_median(x: np.ndarray) -> float:
    return float(np.median(x)) if x.size else float("nan")


def safe_min(x: np.ndarray) -> float:
    return float(np.min(x)) if x.size else float("nan")


def safe_max(x: np.ndarray) -> float:
    return float(np.max(x)) if x.size else float("nan")


def compute_subject_features(subject_dir: Path) -> Optional[SubjectFeatures]:
    subject_id = subject_dir.name

    t1_path = subject_dir / "T1w.nii.gz"
    brain_path = subject_dir / "T1w_brain.nii.gz"
    mask_path = subject_dir / "T1w_brain_mask.nii.gz"
    pve0_path = subject_dir / "T1w_brain_fast_pve_0.nii.gz"
    pve1_path = subject_dir / "T1w_brain_fast_pve_1.nii.gz"
    pve2_path = subject_dir / "T1w_brain_fast_pve_2.nii.gz"

    required = [t1_path, brain_path, mask_path, pve0_path, pve1_path, pve2_path]
    missing = [p.name for p in required if not p.exists()]
    if missing:
        print(f"[WARN] {subject_id}: missing files: {', '.join(missing)}")
        return None

    t1, _, t1_header = load_nii(t1_path)
    brain, _, brain_header = load_nii(brain_path)
    mask, _, mask_header = load_nii(mask_path)
    pve0, _, pve0_header = load_nii(pve0_path)
    pve1, _, pve1_header = load_nii(pve1_path)
    pve2, _, pve2_header = load_nii(pve2_path)

    vv = voxel_volume_mm3(t1_header)
    mm3_to_ml = 1.0 / 1000.0

    # Use the brain mask as the default volume estimate.
    brain_voxels = mask > 0.5
    brain_volume_ml = float(np.sum(brain_voxels) * vv * mm3_to_ml)

    # FAST partial-volume estimates.
    csf_volume_ml = float(np.sum(pve0) * vv * mm3_to_ml)
    gm_volume_ml = float(np.sum(pve1) * vv * mm3_to_ml)
    wm_volume_ml = float(np.sum(pve2) * vv * mm3_to_ml)

    tissue_sum_volume_ml = csf_volume_ml + gm_volume_ml + wm_volume_ml
    tissue_balance_ml = brain_volume_ml - tissue_sum_volume_ml

    brain_values = brain[brain_voxels]
    if brain_values.size == 0:
        print(f"[WARN] {subject_id}: empty brain mask")
        return None

    return SubjectFeatures(
        subject_id=subject_id,
        voxel_volume_mm3=vv,
        brain_volume_ml=brain_volume_ml,
        csf_volume_ml=csf_volume_ml,
        gm_volume_ml=gm_volume_ml,
        wm_volume_ml=wm_volume_ml,
        tissue_sum_volume_ml=tissue_sum_volume_ml,
        tissue_balance_ml=tissue_balance_ml,
        brain_mean_intensity=safe_mean(brain_values),
        brain_std_intensity=safe_std(brain_values),
        brain_median_intensity=safe_median(brain_values),
        brain_min_intensity=safe_min(brain_values),
        brain_max_intensity=safe_max(brain_values),
    )


def write_csv(rows: List[SubjectFeatures], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "subject_id",
                "voxel_volume_mm3",
                "brain_volume_ml",
                "csf_volume_ml",
                "gm_volume_ml",
                "wm_volume_ml",
                "tissue_sum_volume_ml",
                "tissue_balance_ml",
                "brain_mean_intensity",
                "brain_std_intensity",
                "brain_median_intensity",
                "brain_min_intensity",
                "brain_max_intensity",
            ]
        )
        for r in rows:
            writer.writerow(
                [
                    r.subject_id,
                    f"{r.voxel_volume_mm3:.6f}",
                    f"{r.brain_volume_ml:.6f}",
                    f"{r.csf_volume_ml:.6f}",
                    f"{r.gm_volume_ml:.6f}",
                    f"{r.wm_volume_ml:.6f}",
                    f"{r.tissue_sum_volume_ml:.6f}",
                    f"{r.tissue_balance_ml:.6f}",
                    f"{r.brain_mean_intensity:.6f}",
                    f"{r.brain_std_intensity:.6f}",
                    f"{r.brain_median_intensity:.6f}",
                    f"{r.brain_min_intensity:.6f}",
                    f"{r.brain_max_intensity:.6f}",
                ]
            )


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract basic T1 features from preprocessed subjects.")
    parser.add_argument(
        "--preprocessed-root",
        type=Path,
        default=DEFAULT_PREPROCESSED_ROOT,
        help="Root folder containing sub-*/ preprocessed outputs.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=DEFAULT_OUT_CSV,
        help="Output CSV path.",
    )
    parser.add_argument(
        "--subject",
        nargs="*",
        default=None,
        help="Optional subject IDs to process, e.g. sub-009 sub-010.",
    )
    args = parser.parse_args()

    subjects = sorted(args.preprocessed_root.glob("sub-*"))
    if args.subject:
        wanted = set(args.subject)
        subjects = [s for s in subjects if s.name in wanted]

    if not subjects:
        print(f"[WARN] No subjects found in {args.preprocessed_root}")
        return 0

    rows: List[SubjectFeatures] = []
    for subj_dir in subjects:
        features = compute_subject_features(subj_dir)
        if features is None:
            continue
        rows.append(features)
        print(f"[OK] {features.subject_id}: brain={features.brain_volume_ml:.2f} mL, GM={features.gm_volume_ml:.2f} mL, WM={features.wm_volume_ml:.2f} mL")

    if not rows:
        print("[WARN] No features extracted.")
        return 0

    write_csv(rows, args.out_csv)
    print(f"\nSaved features for {len(rows)} subjects to: {args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())