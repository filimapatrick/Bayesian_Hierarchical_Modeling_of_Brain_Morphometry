#!/usr/bin/env python3
"""
Macro-Morphometric Feature Extraction Pipeline for Heterogeneous Clinical MRI

Extracts robust macro-structural biomarkers designed for real-world clinical scans:
- Brain Parenchymal Fraction (BPF = Brain_Volume / TIV)
- Ventricle-to-Brain Ratio (VBR = Ventricle_Volume / Brain_Volume)
- Evans' Index (Frontal horn width / inner skull diameter)
- Hemispheric Asymmetry Index (Left vs Right volume asymmetry)
- Total Intracranial Volume (TIV in cm3)
- Brain Parenchymal Volume (V_brain in cm3)
- Ventricular CSF Volume (V_ventricle in cm3)
- Signal-to-Noise Ratio (SNR proxy)

Merges all extracted morphometry with clinical and scanner covariates from participants.tsv.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage


def extract_macro_morphometry(nii_path: Path) -> Dict[str, float]:
    """
    Extracts macro-morphometric biomarkers from a 3D or 2D multi-slice anatomical T1w scan.
    """
    img = nib.load(str(nii_path))
    data = img.get_fdata(dtype=np.float32)
    zooms = img.header.get_zooms()[:3]
    vox_vol_cm3 = float(np.prod(zooms) / 1000.0)
    shape = data.shape

    # Handle single-slice / defective scans
    min_dim = min(shape[:3])
    if min_dim < 5:
        return {
            "dim_x": shape[0], "dim_y": shape[1], "dim_z": shape[2],
            "voxel_volume_mm3": float(np.prod(zooms)),
            "brain_volume_cm3": np.nan, "tiv_cm3": np.nan, "ventricle_volume_cm3": np.nan,
            "bpf": np.nan, "vbr": np.nan, "evans_index": np.nan, "asymmetry_index": np.nan,
            "snr_proxy": np.nan, "is_volumetric_valid": 0.0
        }

    # 1. Background filtering & Brain extraction (Otsu-inspired percentile thresholding)
    pos_voxels = data[data > 0]
    if len(pos_voxels) == 0:
        return {
            "dim_x": shape[0], "dim_y": shape[1], "dim_z": shape[2],
            "voxel_volume_mm3": float(np.prod(zooms)),
            "brain_volume_cm3": np.nan, "tiv_cm3": np.nan, "ventricle_volume_cm3": np.nan,
            "bpf": np.nan, "vbr": np.nan, "evans_index": np.nan, "asymmetry_index": np.nan,
            "snr_proxy": np.nan, "is_volumetric_valid": 0.0
        }

    intensity_thresh = np.percentile(pos_voxels, 28)
    cand_mask = data > intensity_thresh

    # Connected component labeling to isolate cerebral mass
    labeled, num_feat = ndimage.label(cand_mask)
    if num_feat == 0:
        return {
            "dim_x": shape[0], "dim_y": shape[1], "dim_z": shape[2],
            "voxel_volume_mm3": float(np.prod(zooms)),
            "brain_volume_cm3": np.nan, "tiv_cm3": np.nan, "ventricle_volume_cm3": np.nan,
            "bpf": np.nan, "vbr": np.nan, "evans_index": np.nan, "asymmetry_index": np.nan,
            "snr_proxy": np.nan, "is_volumetric_valid": 0.0
        }

    sizes = ndimage.sum(cand_mask, labeled, range(num_feat + 1))
    largest_label = int(np.argmax(sizes[1:]) + 1)
    brain_mask = (labeled == largest_label)
    brain_mask = ndimage.binary_fill_holes(brain_mask)

    brain_vox_count = float(np.sum(brain_mask))
    brain_vol_cm3 = brain_vox_count * vox_vol_cm3

    # 2. Total Intracranial Volume (TIV) estimation (Brain + CSF envelope)
    struct = ndimage.generate_binary_structure(3, 1)
    icv_mask = ndimage.binary_dilation(brain_mask, structure=struct, iterations=2)
    icv_mask = ndimage.binary_fill_holes(icv_mask)
    tiv_cm3 = float(np.sum(icv_mask) * vox_vol_cm3)

    # Brain Parenchymal Fraction (BPF)
    bpf = float(brain_vol_cm3 / tiv_cm3) if tiv_cm3 > 0 else np.nan

    # 3. Ventricular CSF Extraction & Ventricle-to-Brain Ratio (VBR)
    brain_intensities = data[brain_mask]
    csf_thresh = np.percentile(brain_intensities, 25)

    # Restrict to periventricular anatomical bounding box
    x_c, y_c, z_c = [s // 2 for s in shape]
    dx, dy, dz = int(shape[0] * 0.28), int(shape[1] * 0.28), int(shape[2] * 0.35)
    central_box = np.zeros_like(brain_mask, dtype=bool)
    central_box[
        max(0, x_c - dx) : min(shape[0], x_c + dx),
        max(0, y_c - dy) : min(shape[1], y_c + dy),
        max(0, z_c - dz) : min(shape[2], z_c + dz)
    ] = True

    ventricle_cand = (data < csf_thresh) & brain_mask & central_box
    v_labeled, v_num = ndimage.label(ventricle_cand)
    if v_num > 0:
        v_sizes = ndimage.sum(ventricle_cand, v_labeled, range(v_num + 1))
        # Exclude tiny isolated noise voxels (< 8 voxels)
        v_mask = np.isin(v_labeled, np.where(v_sizes > 8)[0])
        ventricle_vol_cm3 = float(np.sum(v_mask) * vox_vol_cm3)
    else:
        ventricle_vol_cm3 = 0.0

    vbr = float(ventricle_vol_cm3 / brain_vol_cm3) if brain_vol_cm3 > 0 else np.nan

    # 4. Evans' Index (Frontal horn width / inner skull diameter)
    axial_sums = np.sum(ventricle_cand, axis=(0, 1))
    if np.max(axial_sums) > 0:
        max_slice_idx = int(np.argmax(axial_sums))
        v_slice = ventricle_cand[:, :, max_slice_idx]
        b_slice = icv_mask[:, :, max_slice_idx]

        v_proj_x = np.where(np.sum(v_slice, axis=1) > 0)[0]
        b_proj_x = np.where(np.sum(b_slice, axis=1) > 0)[0]

        if len(v_proj_x) > 1 and len(b_proj_x) > 1:
            w_v = (v_proj_x[-1] - v_proj_x[0]) * zooms[0]
            w_b = (b_proj_x[-1] - b_proj_x[0]) * zooms[0]
            evans_index = float(w_v / w_b) if w_b > 0 else np.nan
        else:
            evans_index = np.nan
    else:
        evans_index = np.nan

    # 5. Hemispheric Asymmetry Index (AI)
    mid_x = shape[0] // 2
    left_brain = float(np.sum(brain_mask[:mid_x, :, :]))
    right_brain = float(np.sum(brain_mask[mid_x:, :, :]))
    asymmetry_index = float((left_brain - right_brain) / (0.5 * (left_brain + right_brain))) if (left_brain + right_brain) > 0 else 0.0

    # 6. Technical SNR Proxy
    bg_mask = (data == 0) | (~icv_mask)
    bg_std = np.std(data[bg_mask]) if np.sum(bg_mask) > 10 else 1.0
    mean_brain_signal = np.mean(data[brain_mask]) if brain_vox_count > 0 else 0.0
    snr_proxy = float(mean_brain_signal / bg_std) if bg_std > 0 else np.nan

    return {
        "dim_x": shape[0],
        "dim_y": shape[1],
        "dim_z": shape[2],
        "voxel_volume_mm3": float(np.prod(zooms)),
        "brain_volume_cm3": brain_vol_cm3,
        "tiv_cm3": tiv_cm3,
        "ventricle_volume_cm3": ventricle_vol_cm3,
        "bpf": bpf,
        "vbr": vbr,
        "evans_index": evans_index,
        "asymmetry_index": asymmetry_index,
        "snr_proxy": snr_proxy,
        "is_volumetric_valid": 1.0
    }


def main():
    parser = argparse.ArgumentParser(description="Extract macro-morphometric biomarkers from multi-center BIDS archive.")
    parser.add_argument("--bids_dir", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/data/bids", help="Path to BIDS root directory.")
    parser.add_argument("--output_csv", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results/tables/macro_features.csv", help="Path for output features CSV.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output CSV.")
    args = parser.parse_args()

    bids_path = Path(args.bids_dir)
    out_csv = Path(args.output_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    tsv_file = bids_path / "participants.tsv"
    if not tsv_file.exists():
        sys.exit(f"Error: participants.tsv not found in {bids_path}")

    participants_df = pd.read_csv(tsv_file, sep="\t")
    print("=" * 70)
    print(f"Macro-Morphometric Feature Extraction Pipeline")
    print(f"Dataset: {len(participants_df)} subjects from {bids_path}")
    print(f"Target Output: {out_csv}")
    print("=" * 70)

    feature_rows = []
    t_start = time.time()

    for idx, (_, row) in enumerate(participants_df.iterrows(), 1):
        pid = row["participant_id"]
        nii_path = bids_path / pid / "anat" / f"{pid}_T1w.nii.gz"

        if not nii_path.exists():
            print(f"[{idx:03d}/{len(participants_df):03d}] ⚠️ Missing T1w file for {pid}")
            continue

        try:
            morph = extract_macro_morphometry(nii_path)
            # Combine image morphometry with participant clinical metadata
            record = dict(row)
            record.update(morph)
            feature_rows.append(record)

            valid_str = "VALID" if morph["is_volumetric_valid"] == 1.0 else "SCOUT"
            print(f"[{idx:03d}/{len(participants_df):03d}] {pid} ({row['diagnosis']:<13}) [{valid_str}]: BPF={morph['bpf']:.3f} | VBR={morph['vbr']:.4f} | Evans={morph['evans_index']:.3f} | BrainVol={morph['brain_volume_cm3']:.1f}cm³")
        except Exception as e:
            print(f"[{idx:03d}/{len(participants_df):03d}] ❌ Error processing {pid}: {e}")

    df_out = pd.DataFrame(feature_rows)
    df_out.to_csv(out_csv, index=False)
    elapsed = time.time() - t_start

    print("\n" + "=" * 70)
    print(f"Extraction Complete in {elapsed:.1f} seconds!")
    print(f"Saved: {out_csv} ({len(df_out)} subjects)")
    print("=" * 70)

    # Print Clinical Group Summary Table
    print("\n📈 Macro-Morphometric Group Summary (Mean ± SD):")
    valid_df = df_out[df_out["is_volumetric_valid"] == 1.0]
    summary = valid_df.groupby("diagnosis")[["bpf", "vbr", "evans_index", "brain_volume_cm3"]].agg(["count", "mean", "std"])
    print(summary.to_string())


if __name__ == "__main__":
    main()
