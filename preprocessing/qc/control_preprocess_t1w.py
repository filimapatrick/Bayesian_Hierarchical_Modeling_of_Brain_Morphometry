#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


REPO_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry")
DEFAULT_SELECTED = Path("/Volumes/MyHDD/ABDN_DATA/selected_t1w_control.tsv")
DEFAULT_QC = Path("/Volumes/MyHDD/ABDN_DATA/qc_t1_flags_control.tsv")
DEFAULT_OUT_ROOT = REPO_ROOT / "data" / "derivatives" / "preprocessed" / "control"


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def load_selected_manifest(path: Path) -> dict[str, Path]:
    if not path.exists():
        raise FileNotFoundError(f"Selected manifest not found: {path}")

    rows: dict[str, Path] = {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"subject_id", "nifti_path"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")

        for row in reader:
            sid = row["subject_id"].strip()
            nifti_path = Path(row["nifti_path"].strip())
            rows[sid] = nifti_path
    return rows


def load_good_subjects(qc_flags_path: Path) -> set[str]:
    if not qc_flags_path.exists():
        raise FileNotFoundError(f"QC flags file not found: {qc_flags_path}")

    good: set[str] = set()
    with qc_flags_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"subject_id", "status"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{qc_flags_path} is missing columns: {sorted(missing)}")

        for row in reader:
            if row["status"].strip().upper() == "GOOD":
                good.add(row["subject_id"].strip())
    return good


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def preprocess_subject(subject_id: str, input_nii: Path, out_root: Path, overwrite: bool = False) -> None:
    subject_out = out_root / subject_id
    subject_out.mkdir(parents=True, exist_ok=True)

    t1w = subject_out / "T1w.nii.gz"
    t1w_json = subject_out / "T1w.json"

    reoriented = subject_out / "T1w_reoriented.nii.gz"
    brain = subject_out / "T1w_brain.nii.gz"
    brain_mask = subject_out / "T1w_brain_mask.nii.gz"
    fast_prefix = subject_out / "T1w_brain_fast"

    provenance = subject_out / "preprocess_provenance.txt"

    if overwrite:
        for p in [
            t1w, t1w_json, reoriented, brain, brain_mask,
            fast_prefix.with_name(f"{fast_prefix.name}_restore.nii.gz"),
            fast_prefix.with_name(f"{fast_prefix.name}_restore_flair.nii.gz"),
            fast_prefix.with_name(f"{fast_prefix.name}_restore_bias.nii.gz"),
            fast_prefix.with_name(f"{fast_prefix.name}_pve_0.nii.gz"),
            fast_prefix.with_name(f"{fast_prefix.name}_pve_1.nii.gz"),
            fast_prefix.with_name(f"{fast_prefix.name}_pve_2.nii.gz"),
            fast_prefix.with_name(f"{fast_prefix.name}_seg.nii.gz"),
            provenance,
        ]:
            if p.exists():
                p.unlink()

    print(f"[RUN] {subject_id}: {input_nii.name}")

    # Copy selected T1 into a stable filename in the preprocessing folder.
    shutil.copy2(input_nii, t1w)
    sidecar = input_nii.with_suffix("").with_suffix(".json")
    if sidecar.exists():
        shutil.copy2(sidecar, t1w_json)

    # Reorient to standard.
    run(["fslreorient2std", str(t1w), str(reoriented)])

    # Brain extraction.
    # -f is intentionally conservative to avoid cutting off anatomy.
    run(["bet", str(reoriented), str(brain), "-m", "-f", "0.30"])

    # FSL BET writes mask as *_mask.nii.gz when given output basename without extension.
    # Depending on platform, the exact filename can vary slightly, so we normalize it.
    bet_mask_candidates = [
        subject_out / "T1w_brain_mask.nii.gz",
        subject_out / "T1w_brain_mask.nii.gz",
        subject_out / "T1w_brain_mask.nii",
        subject_out / "T1w_brain_mask",
    ]
    mask_found = None
    for cand in bet_mask_candidates:
        if cand.exists():
            mask_found = cand
            break

    if mask_found is None:
        # Common BET output is brain_mask.nii.gz based on the output basename.
        alt = subject_out / "T1w_brain_mask.nii.gz"
        if alt.exists():
            mask_found = alt

    if mask_found is None:
        print(f"[WARN] {subject_id}: missing output T1w_brain_mask.nii.gz")
    else:
        if mask_found != brain_mask:
            try:
                shutil.move(str(mask_found), str(brain_mask))
            except Exception:
                pass

    # FAST segmentation on the brain-extracted image.
    run(["fast", "-o", str(fast_prefix), str(brain)])

    # Normalize expected outputs.
    # FAST usually creates:
    #   T1w_brain_fast_restore.nii.gz
    #   T1w_brain_fast_pve_0.nii.gz
    #   T1w_brain_fast_pve_1.nii.gz
    #   T1w_brain_fast_pve_2.nii.gz
    #   T1w_brain_fast_seg.nii.gz
    # The warning below mirrors your earlier workflow if a filename differs.
    expected_fast = [
        subject_out / "T1w_brain_fast_restore.nii.gz",
        subject_out / "T1w_brain_fast_pve_0.nii.gz",
        subject_out / "T1w_brain_fast_pve_1.nii.gz",
        subject_out / "T1w_brain_fast_pve_2.nii.gz",
        subject_out / "T1w_brain_fast_seg.nii.gz",
    ]
    missing_fast = [p.name for p in expected_fast if not p.exists()]
    if missing_fast:
        print(f"[WARN] {subject_id}: missing output {', '.join(missing_fast)}")

    provenance.write_text(
        "\n".join(
            [
                f"subject_id={subject_id}",
                f"input_nifti={input_nii}",
                "steps=fslreorient2std -> bet -> fast",
            ]
        )
        + "\n"
    )

    print(f"[OK] {subject_id}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Preprocess control T1 images with FSL.")
    parser.add_argument(
        "--selected-manifest",
        type=Path,
        default=DEFAULT_SELECTED,
        help="TSV from T1 selection",
    )
    parser.add_argument(
        "--qc-flags",
        type=Path,
        default=DEFAULT_QC,
        help="TSV from visual QC",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=DEFAULT_OUT_ROOT,
        help="Output preprocessing root",
    )
    parser.add_argument(
        "--subject",
        nargs="*",
        default=None,
        help="Optional subject IDs to process",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing outputs",
    )
    args = parser.parse_args()

    selected = load_selected_manifest(args.selected_manifest)
    good_subjects = load_good_subjects(args.qc_flags)

    if args.subject:
        requested = set(args.subject)
        good_subjects = good_subjects.intersection(requested)

    good_subjects = sorted(good_subjects)
    if not good_subjects:
        print("[WARN] No subjects matched QC GOOD + selected manifest.")
        return 0

    args.out_root.mkdir(parents=True, exist_ok=True)

    ok = 0
    failed = 0

    for sid in good_subjects:
        if sid not in selected:
            continue
        try:
            preprocess_subject(
                subject_id=sid,
                input_nii=selected[sid],
                out_root=args.out_root,
                overwrite=args.overwrite,
            )
            ok += 1
        except Exception as e:
            failed += 1
            print(f"[WARN] {sid}: preprocessing failed: {e}")

    print(f"\nDone. OK={ok}, failed={failed}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())