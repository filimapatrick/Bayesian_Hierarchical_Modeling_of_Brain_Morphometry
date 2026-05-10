#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd


REPO_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry")
DEFAULT_PREPROCESSED_ROOT = REPO_ROOT / "data" / "derivatives" / "preprocessed" / "control"
DEFAULT_OUT_CSV = REPO_ROOT / "data" / "derivatives" / "features" / "control" / "t1_features.csv"


def volume_ml_from_mask(mask_path: Path, voxel_volume_mm3: float) -> float:
    img = nib.load(str(mask_path))
    data = img.get_fdata()
    voxels = np.count_nonzero(data > 0)
    return float(voxels * voxel_volume_mm3 / 1000.0)


def load_subject_dir(root: Path) -> list[Path]:
    if not root.exists():
        return []
    return sorted([p for p in root.iterdir() if p.is_dir() and p.name.startswith("sub-")])


def compute_features(subject_dir: Path) -> dict[str, object] | None:
    brain_img_path = subject_dir / "T1w_brain.nii.gz"
    mask_path = subject_dir / "T1w_brain_mask.nii.gz"
    pve0_path = subject_dir / "T1w_brain_fast_pve_0.nii.gz"
    pve1_path = subject_dir / "T1w_brain_fast_pve_1.nii.gz"
    pve2_path = subject_dir / "T1w_brain_fast_pve_2.nii.gz"

    required = [brain_img_path, mask_path, pve0_path, pve1_path, pve2_path]
    missing = [p.name for p in required if not p.exists()]
    if missing:
        print(f"[WARN] {subject_dir.name}: missing files: {', '.join(missing)}")
        return None

    brain_img = nib.load(str(brain_img_path))
    voxel_volume_mm3 = float(np.prod(brain_img.header.get_zooms()[:3]))

    brain_volume_ml = volume_ml_from_mask(mask_path, voxel_volume_mm3)
    gm_volume_ml = volume_ml_from_mask(pve1_path, voxel_volume_mm3)
    wm_volume_ml = volume_ml_from_mask(pve2_path, voxel_volume_mm3)
    csf_volume_ml = volume_ml_from_mask(pve0_path, voxel_volume_mm3)

    tissue_balance_ml = brain_volume_ml - (gm_volume_ml + wm_volume_ml + csf_volume_ml)

    print(
        f"[OK] {subject_dir.name}: "
        f"brain={brain_volume_ml:.2f} mL, "
        f"GM={gm_volume_ml:.2f} mL, "
        f"WM={wm_volume_ml:.2f} mL"
    )

    return {
        "subject_id": subject_dir.name,
        "brain_volume_ml": round(brain_volume_ml, 6),
        "gm_volume_ml": round(gm_volume_ml, 6),
        "wm_volume_ml": round(wm_volume_ml, 6),
        "csf_volume_ml": round(csf_volume_ml, 6),
        "tissue_balance_ml": round(tissue_balance_ml, 6),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract simple volumetric features from control T1 preprocessing outputs.")
    parser.add_argument(
        "--preprocessed-root",
        type=Path,
        default=DEFAULT_PREPROCESSED_ROOT,
        help="Root containing subject preprocessing folders",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=DEFAULT_OUT_CSV,
        help="Output CSV for extracted features",
    )
    args = parser.parse_args()

    subject_dirs = load_subject_dir(args.preprocessed_root)
    if not subject_dirs:
        print(f"[WARN] No subjects found in {args.preprocessed_root}")
        return 0

    rows = []
    for subject_dir in subject_dirs:
        feat = compute_features(subject_dir)
        if feat is not None:
            rows.append(feat)

    if not rows:
        print(f"[WARN] No features could be extracted from {args.preprocessed_root}")
        return 0

    df = pd.DataFrame(rows)
    df = df.sort_values("subject_id").reset_index(drop=True)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)

    print(f"\nSaved features for {len(df)} subjects to: {args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())