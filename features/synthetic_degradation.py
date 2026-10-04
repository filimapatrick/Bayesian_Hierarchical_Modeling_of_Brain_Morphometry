#!/usr/bin/env python3
"""
Experiment 1: Controlled Synthetic Degradation Lab

Tests the hypothesis: "What survives the blur?"
Empirically quantifies how morphometric biomarkers degrade as acquisition resolution
is progressively downsampled and blurred from research-grade (1.0 mm³) to routine
clinical slice thicknesses (3.0 mm, 4.0 mm, 5.0 mm, 6.0 mm).

Uses high-resolution 3D acquisitions (N=35) as an internal ground-truth laboratory,
comparing the resilience of coarse macro-morphometry (Evans' Index, BPF, VBR)
against micro-structural/subcortical segmentation proxies.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage


def extract_all_biomarkers(data: np.ndarray, zooms: Tuple[float, float, float]) -> Dict[str, float]:
    """
    Extracts coarse macro-morphometry and a micro-structural proxy from a 3D MRI volume.
    """
    vox_vol_cm3 = float(np.prod(zooms) / 1000.0)
    shape = data.shape
    
    pos_voxels = data[data > 0]
    if len(pos_voxels) == 0:
        return {
            "bpf": np.nan, "vbr": np.nan, "evans_index": np.nan,
            "asymmetry_index": np.nan, "brain_vol_cm3": np.nan,
            "ventricle_vol_cm3": np.nan, "subcortical_proxy_vol_cm3": np.nan
        }

    # 1. Brain extraction mask via adaptive intensity thresholding
    intensity_thresh = np.percentile(pos_voxels, 28)
    cand_mask = data > intensity_thresh
    labeled, num_feat = ndimage.label(cand_mask)
    if num_feat == 0:
        return {
            "bpf": np.nan, "vbr": np.nan, "evans_index": np.nan,
            "asymmetry_index": np.nan, "brain_vol_cm3": np.nan,
            "ventricle_vol_cm3": np.nan, "subcortical_proxy_vol_cm3": np.nan
        }
        
    sizes = ndimage.sum(cand_mask, labeled, range(num_feat + 1))
    largest_label = int(np.argmax(sizes[1:]) + 1)
    brain_mask = (labeled == largest_label)
    brain_mask = ndimage.binary_fill_holes(brain_mask)
    brain_vol_cm3 = float(np.sum(brain_mask) * vox_vol_cm3)

    # 2. Total Intracranial Volume (TIV) estimation
    struct = ndimage.generate_binary_structure(3, 1)
    icv_mask = ndimage.binary_dilation(brain_mask, structure=struct, iterations=2)
    icv_mask = ndimage.binary_fill_holes(icv_mask)
    tiv_cm3 = float(np.sum(icv_mask) * vox_vol_cm3)
    bpf = float(brain_vol_cm3 / tiv_cm3) if tiv_cm3 > 0 else np.nan

    # 3. Ventricular CSF extraction & VBR
    brain_intensities = data[brain_mask]
    csf_thresh = np.percentile(brain_intensities, 25)

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
        v_mask = np.isin(v_labeled, np.where(v_sizes > 8)[0])
        ventricle_vol_cm3 = float(np.sum(v_mask) * vox_vol_cm3)
    else:
        ventricle_vol_cm3 = 0.0
    vbr = float(ventricle_vol_cm3 / brain_vol_cm3) if brain_vol_cm3 > 0 else np.nan

    # 4. Evans' Index (Frontal horn width / inner skull diameter on max ventricular slice)
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

    # 6. Micro-structural proxy: Medial temporal subcortical cluster
    # Represents a fine anatomical structure (approx. 2-4 cm3) vulnerable to partial volume averaging
    mt_box = np.zeros_like(brain_mask, dtype=bool)
    mt_dx, mt_dy, mt_dz = int(shape[0] * 0.12), int(shape[1] * 0.12), int(shape[2] * 0.08)
    mt_box[
        max(0, x_c - mt_dx) : min(shape[0], x_c + mt_dx),
        max(0, y_c - mt_dy) : min(shape[1], y_c + mt_dy),
        max(0, z_c - 2 * mt_dz) : min(shape[2], z_c)
    ] = True
    subcortical_cand = brain_mask & mt_box & (data > np.percentile(brain_intensities, 40)) & (data < np.percentile(brain_intensities, 60))
    subcortical_vol_cm3 = float(np.sum(subcortical_cand) * vox_vol_cm3)

    return {
        "bpf": bpf,
        "vbr": vbr,
        "evans_index": evans_index,
        "asymmetry_index": asymmetry_index,
        "brain_vol_cm3": brain_vol_cm3,
        "ventricle_vol_cm3": ventricle_vol_cm3,
        "subcortical_proxy_vol_cm3": subcortical_vol_cm3
    }


def degrade_scan(data: np.ndarray, native_zooms: Tuple[float, float, float], target_thickness_mm: float) -> Tuple[np.ndarray, Tuple[float, float, float]]:
    """
    Simulates a thick 2D clinical multi-slice acquisition via slice-profile Gaussian filtering
    and downsampling along the slice-select (z) axis.
    """
    native_sz = float(native_zooms[2])
    if target_thickness_mm <= native_sz:
        return data, native_zooms

    # 1. Apply Gaussian slice profile blur (FWHM = target_thickness_mm)
    sigma = (target_thickness_mm / native_sz) / 2.355
    blurred = ndimage.gaussian_filter1d(data, sigma=sigma, axis=2)

    # 2. Downsample along z-axis to simulate 2D multi-slice stepping
    step = max(1, int(round(target_thickness_mm / native_sz)))
    downsampled = blurred[:, :, ::step]
    degraded_zooms = (float(native_zooms[0]), float(native_zooms[1]), float(native_sz * step))

    return downsampled, degraded_zooms


def run_experiment(bids_dir: Path, output_csv: Path, summary_csv: Path, plot_path: Path):
    tsv_path = bids_dir / "participants.tsv"
    if not tsv_path.exists():
        sys.exit(f"Error: {tsv_path} not found.")

    participants_df = pd.read_csv(tsv_path, sep="\t")
    # Identify high-resolution ground truth candidates (slice thickness <= 1.4 mm, 3D volumetric)
    high_res_df = participants_df[
        (participants_df["slice_thickness_mm"] <= 1.4) &
        (participants_df["mr_acquisition_type"] == "3D")
    ].copy()

    print("=" * 75)
    print("EXPERIMENT 1: CONTROLLED SYNTHETIC DEGRADATION LAB")
    print(f"Cohort: {len(high_res_df)} pristine 1.0 mm³ 3D scans from {bids_dir}")
    print("Degradation levels: 1.0 mm (reference) -> 3.0 mm -> 4.0 mm -> 5.0 mm -> 6.0 mm")
    print("=" * 75)

    target_thicknesses = [3.0, 4.0, 5.0, 6.0]
    records = []
    t_start = time.time()

    for idx, (_, row) in enumerate(high_res_df.iterrows(), 1):
        pid = row["participant_id"]
        diag = row["diagnosis"]
        nii_path = bids_dir / pid / "anat" / f"{pid}_T1w.nii.gz"

        if not nii_path.exists():
            print(f"[{idx:02d}/{len(high_res_df):02d}] ⚠️ Missing {nii_path}")
            continue

        try:
            img = nib.load(str(nii_path))
            data = img.get_fdata(dtype=np.float32)
            zooms = tuple(img.header.get_zooms()[:3])

            # Reference measurements at 1.0 mm
            ref_metrics = extract_all_biomarkers(data, zooms)

            # Record reference row
            ref_record = {
                "participant_id": pid,
                "diagnosis": diag,
                "target_thickness_mm": 1.0,
                "is_reference": True,
                **ref_metrics
            }
            # Bias and relative error are zero for reference
            for m in ["bpf", "vbr", "evans_index", "asymmetry_index", "subcortical_proxy_vol_cm3"]:
                ref_record[f"{m}_bias"] = 0.0
                ref_record[f"{m}_rel_error_pct"] = 0.0
            records.append(ref_record)

            # Evaluate each degraded level
            deg_summary = []
            for target_h in target_thicknesses:
                deg_data, deg_zooms = degrade_scan(data, zooms, target_h)
                deg_metrics = extract_all_biomarkers(deg_data, deg_zooms)

                deg_record = {
                    "participant_id": pid,
                    "diagnosis": diag,
                    "target_thickness_mm": target_h,
                    "is_reference": False,
                    **deg_metrics
                }

                for m in ["bpf", "vbr", "evans_index", "asymmetry_index", "subcortical_proxy_vol_cm3"]:
                    ref_val = ref_metrics[m]
                    deg_val = deg_metrics[m]
                    if np.isnan(ref_val) or np.isnan(deg_val) or ref_val == 0:
                        bias = np.nan
                        rel_err = np.nan
                    else:
                        bias = deg_val - ref_val
                        rel_err = abs(deg_val - ref_val) / abs(ref_val) * 100.0

                    deg_record[f"{m}_bias"] = bias
                    deg_record[f"{m}_rel_error_pct"] = rel_err

                records.append(deg_record)
                ei_err = deg_record.get("evans_index_rel_error_pct", np.nan)
                bpf_err = deg_record.get("bpf_rel_error_pct", np.nan)
                deg_summary.append(f"{target_h}mm(EI:{ei_err:.1f}%,BPF:{bpf_err:.1f}%)")

            print(f"[{idx:02d}/{len(high_res_df):02d}] ✓ {pid} ({diag:<13}): {' | '.join(deg_summary)}")

        except Exception as e:
            print(f"[{idx:02d}/{len(high_res_df):02d}] ❌ Error processing {pid}: {e}")

    df_results = pd.DataFrame(records)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df_results.to_csv(output_csv, index=False)
    print(f"\n✓ Saved trial results: {output_csv} ({len(df_results)} rows)")

    # Compute summary statistics by target thickness
    degraded_subset = df_results[df_results["target_thickness_mm"] > 1.0]
    metrics = ["evans_index", "bpf", "vbr", "subcortical_proxy_vol_cm3"]
    summary_rows = []

    for h, group in degraded_subset.groupby("target_thickness_mm"):
        row = {"thickness_mm": h}
        for m in metrics:
            err_col = f"{m}_rel_error_pct"
            vals = group[err_col].dropna()
            row[f"{m}_mean_err"] = float(np.mean(vals))
            row[f"{m}_std_err"] = float(np.std(vals))
            row[f"{m}_median_err"] = float(np.median(vals))
            row[f"{m}_q25_err"] = float(np.percentile(vals, 25))
            row[f"{m}_q75_err"] = float(np.percentile(vals, 75))
        summary_rows.append(row)

    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(summary_csv, index=False)
    print(f"✓ Saved summary table: {summary_csv}")

    print("\n" + "=" * 75)
    print("EXPERIMENT 1 DEGRADATION SUMMARY: RELATIVE PERCENTAGE ERROR (% ± SD)")
    print("=" * 75)
    print(f"{'Thickness':<10} | {'Evans Index (%)':<18} | {'BPF (%)':<16} | {'VBR (%)':<16} | {'Subcortical (%)':<18}")
    print("-" * 85)
    for _, r in df_summary.iterrows():
        print(f"{r['thickness_mm']:<10.1f} | {r['evans_index_mean_err']:>6.2f} ± {r['evans_index_std_err']:<6.2f}%   | {r['bpf_mean_err']:>5.2f} ± {r['bpf_std_err']:<5.2f}% | {r['vbr_mean_err']:>5.2f} ± {r['vbr_std_err']:<5.2f}% | {r['subcortical_proxy_vol_cm3_mean_err']:>6.2f} ± {r['subcortical_proxy_vol_cm3_std_err']:<6.2f}%")
    print("=" * 75)

    # Generate Publication Figure
    generate_degradation_plot(df_results, df_summary, plot_path)
    print(f"✓ Saved publication figure: {plot_path}")
    print(f"Experiment 1 completed in {time.time() - t_start:.1f}s")


def generate_degradation_plot(df_results: pd.DataFrame, df_summary: pd.DataFrame, plot_path: Path):
    """
    Renders high-resolution publication-quality figure:
    Figure 1: Degradation Error Curves and Macro vs Micro Biomarker Resilience
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        print(f"⚠️ Matplotlib or Seaborn not installed. Skipping figure generation to {plot_path}.")
        print("   To generate the plot, run: pip install matplotlib seaborn")
        return

    plot_path.parent.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", font_scale=1.1)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={"width_ratios": [1.2, 1.0]})

    # Panel A: Degradation Error Curves (Mean ± SE vs Slice Thickness)
    ax1 = axes[0]
    thicknesses = [1.0] + list(df_summary["thickness_mm"])

    # Line colors & markers
    palette = {
        "Evans' Index": ("#2ca02c", "o", "-"),
        "Brain Parenchymal Fraction (BPF)": ("#1f77b4", "s", "-"),
        "Ventricle-to-Brain Ratio (VBR)": ("#ff7f0e", "^", "--"),
        "Subcortical Micro-Proxy": ("#d62728", "D", ":")
    }

    # Evans Index
    ei_errs = [0.0] + list(df_summary["evans_index_mean_err"])
    ei_se = [0.0] + list(df_summary["evans_index_std_err"] / np.sqrt(35))
    ax1.errorbar(thicknesses, ei_errs, yerr=ei_se, label="Evans' Index (In-Plane Transverse)",
                 color=palette["Evans' Index"][0], marker=palette["Evans' Index"][1],
                 linestyle=palette["Evans' Index"][2], linewidth=2.5, markersize=8, capsize=4)

    # BPF
    bpf_errs = [0.0] + list(df_summary["bpf_mean_err"])
    bpf_se = [0.0] + list(df_summary["bpf_std_err"] / np.sqrt(35))
    ax1.errorbar(thicknesses, bpf_errs, yerr=bpf_se, label="Brain Parenchymal Fraction (BPF)",
                 color=palette["Brain Parenchymal Fraction (BPF)"][0], marker=palette["Brain Parenchymal Fraction (BPF)"][1],
                 linestyle=palette["Brain Parenchymal Fraction (BPF)"][2], linewidth=2.5, markersize=8, capsize=4)

    # VBR
    vbr_errs = [0.0] + list(df_summary["vbr_mean_err"])
    vbr_se = [0.0] + list(df_summary["vbr_std_err"] / np.sqrt(35))
    ax1.errorbar(thicknesses, vbr_errs, yerr=vbr_se, label="Ventricle-to-Brain Ratio (VBR)",
                 color=palette["Ventricle-to-Brain Ratio (VBR)"][0], marker=palette["Ventricle-to-Brain Ratio (VBR)"][1],
                 linestyle=palette["Ventricle-to-Brain Ratio (VBR)"][2], linewidth=2.5, markersize=8, capsize=4)

    # Subcortical Micro-Proxy
    sub_errs = [0.0] + list(df_summary["subcortical_proxy_vol_cm3_mean_err"])
    sub_se = [0.0] + list(df_summary["subcortical_proxy_vol_cm3_std_err"] / np.sqrt(35))
    ax1.errorbar(thicknesses, sub_errs, yerr=sub_se, label="Subcortical Micro-Structure Proxy",
                 color=palette["Subcortical Micro-Proxy"][0], marker=palette["Subcortical Micro-Proxy"][1],
                 linestyle=palette["Subcortical Micro-Proxy"][2], linewidth=2.5, markersize=8, capsize=4)

    # Annotations & styling
    ax1.axvspan(4.0, 6.0, color="#f0f0f0", alpha=0.6, label="Routine Nigerian Clinical Regime (4-6 mm)")
    ax1.set_xlabel("Acquisition Slice Thickness (mm)", fontweight="bold")
    ax1.set_ylabel("Mean Relative Measurement Error (%)", fontweight="bold")
    ax1.set_title("A. Error Degradation Curves across Scan Resolutions", fontweight="bold", pad=12)
    ax1.set_xlim(0.8, 6.2)
    ax1.set_ylim(-2, 105)
    ax1.legend(loc="upper left", frameon=True, fontsize=9)

    # Panel B: Distribution of Measurement Error at 5.0 mm Routine Clinical Thickness
    ax2 = axes[1]
    df_5mm = df_results[df_results["target_thickness_mm"] == 5.0].copy()

    plot_data = pd.DataFrame({
        "Biomarker": (
            ["Evans' Index"] * len(df_5mm) +
            ["BPF"] * len(df_5mm) +
            ["VBR"] * len(df_5mm) +
            ["Subcortical Proxy"] * len(df_5mm)
        ),
        "Relative Error (%)": (
            list(df_5mm["evans_index_rel_error_pct"]) +
            list(df_5mm["bpf_rel_error_pct"]) +
            list(df_5mm["vbr_rel_error_pct"]) +
            list(df_5mm["subcortical_proxy_vol_cm3_rel_error_pct"])
        )
    })

    box_colors = ["#2ca02c", "#1f77b4", "#ff7f0e", "#d62728"]
    sns.boxplot(data=plot_data, x="Biomarker", y="Relative Error (%)", ax=ax2, palette=box_colors, width=0.5, fliersize=3)
    sns.stripplot(data=plot_data, x="Biomarker", y="Relative Error (%)", ax=ax2, color="black", alpha=0.35, jitter=0.2, size=5)

    ax2.set_xlabel("Biomarker Category", fontweight="bold")
    ax2.set_ylabel("Relative Error (%) at 5.0 mm", fontweight="bold")
    ax2.set_title("B. Measurement Resilience at Routine Clinical 5.0 mm", fontweight="bold", pad=12)
    ax2.set_ylim(-5, 110)

    # Add threshold callout line
    ax2.axhline(15.0, color="gray", linestyle="--", alpha=0.7)
    ax2.text(0.1, 17.0, "15% Robustness Ceiling", color="dimgray", fontsize=9, fontstyle="italic")

    plt.tight_layout()
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Run Experiment 1: Controlled Synthetic Degradation Lab.")
    parser.add_argument("--bids_dir", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/data/bids")
    parser.add_argument("--output_csv", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results/tables/experiment1_synthetic_degradation.csv")
    parser.add_argument("--summary_csv", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results/tables/experiment1_degradation_summary.csv")
    parser.add_argument("--plot_path", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results/figures/figure1_synthetic_degradation_curves.png")
    args = parser.parse_args()

    run_experiment(
        bids_dir=Path(args.bids_dir),
        output_csv=Path(args.output_csv),
        summary_csv=Path(args.summary_csv),
        plot_path=Path(args.plot_path)
    )


if __name__ == "__main__":
    main()
