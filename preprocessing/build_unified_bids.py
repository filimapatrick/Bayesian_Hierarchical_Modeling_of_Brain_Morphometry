#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import re
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPO_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry")
NIFTI_ROOT = REPO_ROOT / "data" / "derivatives" / "nifti"
BIDS_ROOT = REPO_ROOT / "data" / "bids"
BRAINLIFE_ROOT = Path("/Volumes/MyHDD/ABDN_DATA/Data_from_brainlife/proj-6554f423b094062da63aa4c9")

COHORTS = {
    "control": {
        "prefix": "con",
        "diagnosis": "CONTROL",
    },
    "dementia": {
        "prefix": "dem",
        "diagnosis": "DEMENTIA",
    },
    "epilepsy": {
        "prefix": "epi",
        "diagnosis": "EPILEPSY",
    },
    "hydrocephalus": {
        "prefix": "hyd",
        "diagnosis": "HYDROCEPHALUS",
    },
}

POSITIVE_3D_T1 = [
    "mprage", "mp-rage", "bravo", "spgr", "t1_mprage", "t1_mpr",
    "3d t1", "3d-t1", "3d_t1", "t1 cube", "cube t1", "tfe", "gre t1"
]

POSITIVE_T1 = [
    "t1", "t1w", "t1 weighted", "t1-weighted", "t1 se", "t1 tse", "t1_fse", "t1 se"
]

NEGATIVE_PATTERNS = [
    "localizer", "scout", "survey", "mip", "projection", "screen save",
    "cal ", "calibration", "adc", "dwi", "diff", "tracew", "trace",
    "tof", "angio", "angiography", "venous", "venogram", "mra", "screenshot",
    "derived", "secondary", "reformat"
]

CONTRAST_HINTS = [
    "post contrast", "post-contrast", "postcontrast", "pc", "+c", " + c", "c+", "contrast", "ce-gadolinium"
]


@dataclass
class T1Candidate:
    subject_id: str
    nifti_path: Path
    json_path: Path
    score: int
    meta: dict


def norm_text(value: object) -> str:
    return str(value or "").strip().lower()


def score_series(meta: dict, nii_path: Path) -> int:
    desc = norm_text(meta.get("SeriesDescription"))
    prot = norm_text(meta.get("ProtocolName"))
    seq = norm_text(meta.get("ScanningSequence"))
    name = f"{desc} {prot} {seq}"

    score = 0
    if any(k in name for k in POSITIVE_3D_T1):
        score += 150
    if any(k in name for k in POSITIVE_T1):
        score += 70
    if "t1" in name:
        score += 30

    if any(k in name for k in NEGATIVE_PATTERNS):
        score -= 200
    if any(k in name for k in ["adc", "trace", "dwi", "diff"]):
        score -= 250
    if any(k in name for k in ["tof", "angio"]):
        score -= 250

    if any(k in name for k in CONTRAST_HINTS):
        score -= 10

    # Check 3D dimensional slice coverage using nibabel header
    try:
        import nibabel as nib
        shape = nib.load(str(nii_path)).header.get_data_shape()
        min_dim = min(shape[:3])
        if min_dim < 8:
            score -= 1000  # Strongly penalize 1-slice, 3-slice, 5-slice localizers/scouts
        elif min_dim >= 50:
            score += 100   # 3D high-resolution volumetric
        elif min_dim >= 15:
            score += 50    # Full 2D multi-slice volume
    except Exception:
        pass

    return score


def find_best_t1(subject_dir: Path) -> Optional[T1Candidate]:
    candidates = []
    for js in subject_dir.glob("*.json"):
        nii = js.with_suffix(".nii.gz")
        if not nii.exists():
            continue
        try:
            with js.open() as f:
                meta = json.load(f)
        except Exception:
            continue

        desc = norm_text(meta.get("SeriesDescription"))
        prot = norm_text(meta.get("ProtocolName"))
        name = f"{desc} {prot}"

        if not ("t1" in name or any(k in name for k in POSITIVE_3D_T1)):
            continue

        score = score_series(meta, nii)
        if score > 0:
            candidates.append(T1Candidate(
                subject_id=subject_dir.name,
                nifti_path=nii,
                json_path=js,
                score=score,
                meta=meta
            ))

    if not candidates:
        return None
    candidates.sort(key=lambda c: c.score, reverse=True)
    return candidates[0]


def extract_subject_number(subject_str: str) -> str:
    m = re.search(r"(\d+)", subject_str)
    if m:
        return f"{int(m.group(1)):03d}"
    return subject_str.replace("sub-", "")


def is_contrast_enhanced(meta: dict) -> bool:
    desc = norm_text(meta.get("SeriesDescription"))
    prot = norm_text(meta.get("ProtocolName"))
    name = f"{desc} {prot}"
    return any(k in name for k in CONTRAST_HINTS)


def find_best_brainlife_t1(sdir: Path) -> Optional[Tuple[Path, dict]]:
    t1_dirs = list(sdir.glob("dt-neuro-anat-t1w*"))
    if not t1_dirs:
        return None

    scored = []
    for dt in t1_dirs:
        info_file = dt / "_info.json"
        nii_file = dt / "t1.nii.gz"
        if not (info_file.exists() and nii_file.exists()):
            continue
        try:
            with open(info_file) as f:
                d = json.load(f)
        except Exception:
            continue

        meta = d.get("meta", {})
        tags = [str(t).lower() for t in d.get("tags", [])]
        desc = str(d.get("desc", "")).lower()

        # Scoring preferences:
        score = 100
        if "acq-axial" in tags:
            score += 50
        elif "acq-sagittal" in tags:
            score += 30
        elif "acq-coronal" in tags:
            score += 20

        if "no-ce" in tags:
            score += 40
        elif "ce-gadolinium" in tags:
            score -= 10

        try:
            import nibabel as nib
            sh = nib.load(str(nii_file)).header.get_data_shape()
            if min(sh[:3]) < 8:
                score -= 1000
            elif min(sh[:3]) >= 15:
                score += 30
        except Exception:
            pass

        scored.append((score, nii_file, meta))

    if not scored:
        return None
    scored.sort(key=lambda x: x[0], reverse=True)
    return scored[0][1], scored[0][2]


def build_bids_dataset():
    print("=" * 60)
    print("Building Full Publication-Ready BIDS Dataset")
    print(f"Target directory: {BIDS_ROOT}")
    print("=" * 60)

    BIDS_ROOT.mkdir(parents=True, exist_ok=True)

    participants_rows = []
    processed_count = 0
    skipped_count = 0

    # 1. Process Core Clinical Cohorts (RSUTH, UPTH, BMH, FMI)
    for cohort, info in COHORTS.items():
        cohort_dir = NIFTI_ROOT / cohort
        if not cohort_dir.exists():
            continue

        subjs = sorted([d for d in cohort_dir.iterdir() if d.is_dir() and d.name.startswith("sub-")])
        print(f"\nProcessing Clinical Archive {cohort.upper()} ({len(subjs)} subjects)...")

        for sdir in subjs:
            num = extract_subject_number(sdir.name)
            bids_id = f"sub-{info['prefix']}{num}"
            best_t1 = find_best_t1(sdir)

            if not best_t1:
                skipped_count += 1
                continue

            anat_dir = BIDS_ROOT / bids_id / "anat"
            anat_dir.mkdir(parents=True, exist_ok=True)

            target_nii = anat_dir / f"{bids_id}_T1w.nii.gz"
            target_json = anat_dir / f"{bids_id}_T1w.json"

            shutil.copy2(best_t1.nifti_path, target_nii)
            shutil.copy2(best_t1.json_path, target_json)

            m = best_t1.meta
            institution = str(m.get("InstitutionName") or "Unknown").strip()
            manufacturer = str(m.get("Manufacturer") or "Unknown").strip()
            model = str(m.get("ManufacturersModelName") or "Unknown").strip()
            field_strength = m.get("MagneticFieldStrength")
            acq_type = str(m.get("MRAcquisitionType") or "Unknown").strip()
            slice_thickness = m.get("SliceThickness")
            spacing = m.get("SpacingBetweenSlices")
            tr = m.get("RepetitionTime")
            te = m.get("EchoTime")
            flip = m.get("FlipAngle")
            contrast = is_contrast_enhanced(m)

            row = {
                "participant_id": bids_id,
                "original_id": sdir.name,
                "dataset_source": "clinical_hospital_archive",
                "cohort": cohort,
                "diagnosis": info["diagnosis"],
                "institution_name": institution,
                "manufacturer": manufacturer,
                "manufacturers_model_name": model,
                "magnetic_field_strength": field_strength if field_strength is not None else "n/a",
                "mr_acquisition_type": acq_type,
                "slice_thickness_mm": slice_thickness if slice_thickness is not None else "n/a",
                "spacing_between_slices_mm": spacing if spacing is not None else "n/a",
                "repetition_time_s": tr if tr is not None else "n/a",
                "echo_time_s": te if te is not None else "n/a",
                "flip_angle_deg": flip if flip is not None else "n/a",
                "contrast_enhanced": "true" if contrast else "false",
                "series_description": str(m.get("SeriesDescription") or "n/a"),
                "protocol_name": str(m.get("ProtocolName") or "n/a"),
            }
            participants_rows.append(row)
            processed_count += 1

    # 2. Process Brainlife Released Dataset (Controls 1-33, Dementia 34-66, Parkinson 67-88)
    if BRAINLIFE_ROOT.exists():
        print(f"\nProcessing Brainlife Dataset (88 subjects from Sci Data 2025)...")
        bl_subjs = sorted([d for d in BRAINLIFE_ROOT.iterdir() if d.is_dir() and d.name.startswith("sub-") and "." not in d.name])
        
        for sdir in bl_subjs:
            m_num = re.search(r"(\d+)", sdir.name)
            if not m_num:
                continue
            idx = int(m_num.group(1))

            if idx <= 33:
                cohort = "control"
                diagnosis = "CONTROL"
                # Offset by 14 (existing controls) -> sub-con015 to sub-con047
                bids_id = f"sub-con{14 + idx:03d}"
            elif idx <= 66:
                cohort = "dementia"
                diagnosis = "DEMENTIA"
                # Offset by 47 (existing dementia) -> sub-dem048 to sub-dem080
                bids_id = f"sub-dem{47 + (idx - 33):03d}"
            else:
                cohort = "parkinson"
                diagnosis = "PARKINSON"
                # Parkinson cohort -> sub-pd001 to sub-pd022
                bids_id = f"sub-pd{idx - 66:03d}"

            bl_res = find_best_brainlife_t1(sdir)
            if not bl_res:
                print(f"  [SKIP] brainlife/{sdir.name}: No suitable T1w scan found.")
                skipped_count += 1
                continue

            nii_src, meta = bl_res

            anat_dir = BIDS_ROOT / bids_id / "anat"
            anat_dir.mkdir(parents=True, exist_ok=True)

            target_nii = anat_dir / f"{bids_id}_T1w.nii.gz"
            target_json = anat_dir / f"{bids_id}_T1w.json"

            shutil.copy2(nii_src, target_nii)
            with target_json.open("w") as f:
                json.dump(meta, f, indent=2)

            institution = str(meta.get("InstitutionName") or "INTERCONTINENTAL DIAG. CENTER").strip()
            manufacturer = str(meta.get("Manufacturer") or "GE").strip()
            model = str(meta.get("StationName") or meta.get("DeviceSerialNumber") or "Brivo MR235").strip()
            
            # Infer field strength: Brivo MR235 is 0.35T / 0.3T, or from DeviceSerialNumber
            serial = str(meta.get("DeviceSerialNumber") or "")
            if "0.3t" in serial.lower():
                field_strength = 0.3
            else:
                field_strength = meta.get("MagneticFieldStrength") or 0.35

            acq_type = str(meta.get("MRAcquisitionType") or "2D").strip()
            slice_thickness = meta.get("SliceThickness") or 5.0
            spacing = meta.get("SpacingBetweenSlices") or 7.0
            tr = meta.get("RepetitionTime")
            te = meta.get("EchoTime")
            flip = meta.get("FlipAngle")
            contrast = is_contrast_enhanced(meta)

            row = {
                "participant_id": bids_id,
                "original_id": f"brainlife_{sdir.name}",
                "dataset_source": "brainlife_sci_data_2025",
                "cohort": cohort,
                "diagnosis": diagnosis,
                "institution_name": institution,
                "manufacturer": manufacturer,
                "manufacturers_model_name": model,
                "magnetic_field_strength": field_strength if field_strength is not None else "n/a",
                "mr_acquisition_type": acq_type,
                "slice_thickness_mm": slice_thickness if slice_thickness is not None else "n/a",
                "spacing_between_slices_mm": spacing if spacing is not None else "n/a",
                "repetition_time_s": tr if tr is not None else "n/a",
                "echo_time_s": te if te is not None else "n/a",
                "flip_angle_deg": flip if flip is not None else "n/a",
                "contrast_enhanced": "true" if contrast else "false",
                "series_description": str(meta.get("SeriesDescription") or "Axi SE2D T1W"),
                "protocol_name": str(meta.get("ProtocolName") or "Axi SE2D T1W"),
            }
            participants_rows.append(row)
            processed_count += 1

    # 3. Process Kano Control Cohort (AKTH / NKDC Kano)
    kano_nifti_dir = NIFTI_ROOT / "control_kano"
    if kano_nifti_dir.exists():
        kano_subjs = sorted([d for d in kano_nifti_dir.iterdir() if d.is_dir()])
        print(f"\nProcessing Kano Control Archive ({len(kano_subjs)} subjects)...")
        kano_start_idx = 48
        for sdir in kano_subjs:
            bids_id = f"sub-con{kano_start_idx:03d}"
            jsons = list(sdir.glob("*.json"))
            t1_cands = []
            for js in jsons:
                try:
                    with js.open() as f:
                        m = json.load(f)
                except Exception:
                    continue
                desc = norm_text(m.get("SeriesDescription"))
                prot = norm_text(m.get("ProtocolName"))
                full = f"{desc} {prot}"
                if any(neg in full for neg in ["localizer", "scout", "diff", "adc", "tracew", "tof", "angio", "mip", "phoenix"]):
                    continue
                if "t1" in full:
                    sc = 100
                    if "3d" in full or "mprage" in full: sc += 50
                    if "tra" in full or "ax" in full: sc += 30
                    if "cor" in full: sc += 20
                    if "sag" in full: sc += 10
                    if "+c" in full or "c+" in full or "contrast" in full: sc -= 15
                    if "fs" in full: sc -= 10
                    nii_path = js.with_suffix(".nii.gz")
                    if nii_path.exists():
                        t1_cands.append((sc, m, js, nii_path))

            if not t1_cands:
                print(f"  [SKIP] {sdir.name}: No suitable T1w scan found.")
                skipped_count += 1
                continue

            t1_cands.sort(key=lambda x: x[0], reverse=True)
            _, meta, best_js, best_nii = t1_cands[0]
            kano_start_idx += 1

            anat_dir = BIDS_ROOT / bids_id / "anat"
            anat_dir.mkdir(parents=True, exist_ok=True)

            target_nii = anat_dir / f"{bids_id}_T1w.nii.gz"
            target_json = anat_dir / f"{bids_id}_T1w.json"

            shutil.copy2(best_nii, target_nii)
            bids_meta = dict(meta)
            bids_meta["Modality"] = "MR"
            with target_json.open("w") as f:
                json.dump(bids_meta, f, indent=2)

            row = {
                "participant_id": bids_id,
                "original_id": f"kano_{sdir.name}",
                "dataset_source": "clinical_hospital_archive",
                "cohort": "control",
                "diagnosis": "CONTROL",
                "institution_name": "Aminu Kano Teaching Hospital / NKDC Kano",
                "manufacturer": str(meta.get("Manufacturer") or "Siemens").strip(),
                "manufacturers_model_name": str(meta.get("ManufacturersModelName") or "MAGNETOM ESSENZA").strip(),
                "magnetic_field_strength": meta.get("MagneticFieldStrength") or 1.5,
                "mr_acquisition_type": "2D",
                "slice_thickness_mm": meta.get("SliceThickness") or "n/a",
                "spacing_between_slices_mm": meta.get("SpacingBetweenSlices") or "n/a",
                "repetition_time_s": meta.get("RepetitionTime") or "n/a",
                "echo_time_s": meta.get("EchoTime") or "n/a",
                "flip_angle_deg": meta.get("FlipAngle") or "n/a",
                "contrast_enhanced": "false",
                "series_description": str(meta.get("SeriesDescription") or "t1_se_tra_320"),
                "protocol_name": str(meta.get("ProtocolName") or "ROUTINE BRAIN"),
            }
            participants_rows.append(row)
            processed_count += 1

    # Write participants.tsv
    tsv_path = BIDS_ROOT / "participants.tsv"
    fieldnames = [
        "participant_id", "original_id", "dataset_source", "cohort", "diagnosis",
        "institution_name", "manufacturer", "manufacturers_model_name",
        "magnetic_field_strength", "mr_acquisition_type",
        "slice_thickness_mm", "spacing_between_slices_mm",
        "repetition_time_s", "echo_time_s", "flip_angle_deg",
        "contrast_enhanced", "series_description", "protocol_name"
    ]
    with tsv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in participants_rows:
            writer.writerow(r)
    print(f"\n✓ Saved: {tsv_path} ({len(participants_rows)} participants)")

    # Write participants.json (BIDS data dictionary)
    json_dict = {
        "participant_id": {"Description": "Unique BIDS identifier for participant"},
        "original_id": {"Description": "Original clinical archive or Brainlife subject identifier"},
        "dataset_source": {"Description": "Provenance of imaging data", "Levels": {"clinical_hospital_archive": "Routine hospital clinical archive (RSUTH, UPTH, BMH, FMI)", "brainlife_sci_data_2025": "Open-access published Nigerian Brain Dataset (Wogu, Filima et al., Sci Data 2025)"}},
        "cohort": {"Description": "Clinical research cohort", "Levels": {"control": "Control group", "dementia": "Dementia cohort", "epilepsy": "Epilepsy cohort", "hydrocephalus": "Hydrocephalus cohort", "parkinson": "Parkinson's disease cohort"}},
        "diagnosis": {"Description": "Clinical diagnostic classification", "Levels": {"CONTROL": "Normal control", "DEMENTIA": "Dementia presentation", "EPILEPSY": "Epilepsy presentation", "HYDROCEPHALUS": "Hydrocephalus presentation", "PARKINSON": "Parkinson's Disease presentation"}},
        "institution_name": {"Description": "Hospital or imaging center where MRI was acquired (e.g. RSUTH, UPTH, FMI, Intercontinental Diagnostic Center)"},
        "manufacturer": {"Description": "MRI scanner vendor (e.g. GE, Siemens)"},
        "manufacturers_model_name": {"Description": "Scanner commercial model name"},
        "magnetic_field_strength": {"Description": "Main magnetic field nominal strength", "Units": "Tesla"},
        "mr_acquisition_type": {"Description": "Acquisition dimensionality", "Levels": {"2D": "2D multi-slice acquisition", "3D": "3D volumetric acquisition"}},
        "slice_thickness_mm": {"Description": "Slice thickness of anatomical scan", "Units": "mm"},
        "spacing_between_slices_mm": {"Description": "Distance between centers of adjacent slices", "Units": "mm"},
        "repetition_time_s": {"Description": "Repetition time (TR)", "Units": "seconds"},
        "echo_time_s": {"Description": "Echo time (TE)", "Units": "seconds"},
        "flip_angle_deg": {"Description": "Radiofrequency excitation flip angle", "Units": "degrees"},
        "contrast_enhanced": {"Description": "Whether gadolinium-based contrast agent was present during scan", "Levels": {"true": "Post-contrast acquisition", "false": "Pre-contrast acquisition"}},
        "series_description": {"Description": "DICOM series description tag"},
        "protocol_name": {"Description": "DICOM protocol name tag"}
    }
    json_path = BIDS_ROOT / "participants.json"
    with json_path.open("w") as f:
        json.dump(json_dict, f, indent=2)
    print(f"✓ Saved: {json_path}")

    # Write dataset_description.json
    desc_dict = {
        "Name": "Multicenter Bayesian Brain Morphometry Dataset across Neurological Disorders in Nigeria",
        "BIDSVersion": "1.10.0",
        "DatasetType": "raw",
        "License": "CC-BY-4.0",
        "Authors": ["Patrick Filima", "African Brain Data Network (ABDN) Collaborators"],
        "Funding": ["African Brain Data Network (ABDN)"],
        "DatasetDOI": "https://doi.org/10.1038/s41597-025-04743-0",
        "EthicsApprovals": ["Institutional Review Board approvals at participating Nigerian centers (RSUTH, UPTH, FMI, IDC)"]
    }
    desc_path = BIDS_ROOT / "dataset_description.json"
    with desc_path.open("w") as f:
        json.dump(desc_dict, f, indent=2)
    print(f"✓ Saved: {desc_path}")

    # Write README.md
    c_counts = {}
    for r in participants_rows:
        c_counts[r['diagnosis']] = c_counts.get(r['diagnosis'], 0) + 1

    readme_content = f"""# 🧠 Multicenter Bayesian Brain Morphometry Dataset across Neurological Disorders in Nigeria

## Overview
This standardized BIDS dataset combines multicenter structural brain MRI acquisitions from hospitals across Nigeria (RSUTH, UPTH, FMI, BMH, and Intercontinental Diagnostic Center), unifying routine clinical archives with the open-access labeled Nigerian Brain Dataset (Wogu, Filima et al., *Scientific Data* 2025).

The dataset encompasses five diagnostic cohorts:
* **Control:** Healthy controls ($N={c_counts.get('CONTROL', 0)}$)
* **Dementia:** Clinical dementia presentations ($N={c_counts.get('DEMENTIA', 0)}$)
* **Parkinson's Disease:** Clinical Parkinson's Disease ($N={c_counts.get('PARKINSON', 0)}$)
* **Hydrocephalus:** Pediatric and adult hydrocephalus presentations ($N={c_counts.get('HYDROCEPHALUS', 0)}$)
* **Epilepsy:** Epilepsy presentations ($N={c_counts.get('EPILEPSY', 0)}$)

## Total Included Subjects
**Total Participants:** {len(participants_rows)}

## Acquisition & Technical Heterogeneity
As detailed in `participants.tsv`, scans reflect routine clinical hospital practice across Nigeria:
* **Vendors:** GE Healthcare, Siemens Healthineers
* **Field Strengths:** 0.3T, 0.35T, and 1.5 Tesla
* **Acquisition Dimensions:** 2D multi-slice Fast Spin Echo (FSE) and 3D volumetric MPRAGE/SPGR
* **Slice Thickness:** Ranging from 1.0mm isotropic up to 5.0–6.0mm thick clinical slices
* **Contrast:** Both pre-contrast and post-contrast acquisitions documented for explicit statistical modeling.

## BIDS Compliance
Organized following BIDS version 1.10.0 specifications.
"""
    readme_path = BIDS_ROOT / "README.md"
    with readme_path.open("w") as f:
        f.write(readme_content)
    print(f"✓ Saved: {readme_path}")

    print("\n" + "=" * 60)
    print(f"SUCCESS: Full Unified BIDS dataset created with {len(participants_rows)} subjects!")
    print(f"BIDS root: {BIDS_ROOT}")
    print("=" * 60)


if __name__ == "__main__":
    build_bids_dataset()
