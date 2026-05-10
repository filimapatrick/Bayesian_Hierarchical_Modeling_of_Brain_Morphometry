#!/usr/bin/env python3
"""
Preprocess selected T1-weighted MRI scans for the dementia pipeline.

Inputs:
- selected_t1w.tsv : one chosen T1 per subject
- qc_t1_flags.tsv  : GOOD/BAD labels from visual QC

Outputs:
- data/derivatives/preprocessed/sub-XXX/
    - T1w.nii.gz
    - T1w.json
    - T1w_reoriented.nii.gz
    - T1w_brain.nii.gz
    - T1w_brain_mask.nii.gz
    - T1w_brain_fast_pve_0.nii.gz
    - T1w_brain_fast_pve_1.nii.gz
    - T1w_brain_fast_pve_2.nii.gz
    - T1w_brain_fast_seg.nii.gz
    - T1w_brain_fast_bias.nii.gz
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set


REPO_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry")
DEFAULT_SELECTED_MANIFEST = Path("/Volumes/MyHDD/ABDN_DATA/selected_t1w.tsv")
DEFAULT_QC_FLAGS = Path("/Volumes/MyHDD/ABDN_DATA/qc_t1_flags.tsv")
DEFAULT_OUT_ROOT = REPO_ROOT / "data" / "derivatives" / "preprocessed"
DEFAULT_WORK_ROOT = REPO_ROOT / "data" / "derivatives" / "preprocessed_work"


@dataclass(frozen=True)
class SelectedRow:
    subject_id: str
    nifti_path: Path
    json_path: Path
    series_description: str
    protocol_name: str
    score: int


def check_tool(name: str) -> None:
    if shutil.which(name) is None:
        raise RuntimeError(f"Required tool not found on PATH: {name}")


def load_good_subjects(qc_flags_path: Path) -> Set[str]:
    good: Set[str] = set()
    with qc_flags_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None or "subject_id" not in reader.fieldnames or "status" not in reader.fieldnames:
            raise ValueError("qc_t1_flags.tsv must contain columns: subject_id, status")

        for row in reader:
            if row["status"].strip().upper() == "GOOD":
                good.add(row["subject_id"].strip())
    return good


def load_selected_manifest(selected_manifest: Path) -> List[SelectedRow]:
    rows: List[SelectedRow] = []
    with selected_manifest.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {
            "subject_id",
            "nifti_path",
            "json_path",
            "series_description",
            "protocol_name",
            "score",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"selected_t1w.tsv is missing columns: {sorted(missing)}")

        for row in reader:
            rows.append(
                SelectedRow(
                    subject_id=row["subject_id"].strip(),
                    nifti_path=Path(row["nifti_path"].strip()),
                    json_path=Path(row["json_path"].strip()),
                    series_description=row["series_description"].strip(),
                    protocol_name=row["protocol_name"].strip(),
                    score=int(row["score"]),
                )
            )
    return rows


def run_cmd(cmd: List[str]) -> None:
    subprocess.run(cmd, check=True)


def preprocess_subject(
    row: SelectedRow,
    out_root: Path,
    work_root: Path,
    overwrite: bool = False,
) -> Dict[str, Path]:
    subj_out = out_root / row.subject_id
    subj_work = work_root / row.subject_id

    if overwrite:
        if subj_out.exists():
            shutil.rmtree(subj_out)
        if subj_work.exists():
            shutil.rmtree(subj_work)

    subj_out.mkdir(parents=True, exist_ok=True)
    subj_work.mkdir(parents=True, exist_ok=True)

    # Copy original selected files into the output folder for traceability.
    t1_copy = subj_out / "T1w.nii.gz"
    json_copy = subj_out / "T1w.json"
    shutil.copy2(row.nifti_path, t1_copy)
    shutil.copy2(row.json_path, json_copy)

    reoriented = subj_work / "T1w_reoriented.nii.gz"
    brain = subj_work / "T1w_brain.nii.gz"
    brain_mask = subj_work / "T1w_brain_mask.nii.gz"
    fast_prefix = subj_work / "T1w_brain_fast"

    # 1) Reorient to standard.
    run_cmd(["fslreorient2std", str(t1_copy), str(reoriented)])

    # 2) Brain extraction.
    # -R helps with robust center estimation.
    # -f is left at the BET default unless you need to tune it later.
    run_cmd(["bet", str(reoriented), str(brain), "-m", "-R"])

    # 3) Bias correction + segmentation.
    # FAST writes:
    #   *_pve_0, *_pve_1, *_pve_2, *_seg, *_bias
    run_cmd(["fast", "-t", "1", "-n", "3", "-o", str(fast_prefix), str(brain)])

    # Move useful outputs into the subject folder.
    outputs: Dict[str, Path] = {
        "T1w": t1_copy,
        "T1w_json": json_copy,
        "T1w_reoriented": reoriented,
        "T1w_brain": brain,
        "T1w_brain_mask": brain_mask,
        "T1w_brain_fast_pve_0": fast_prefix.with_name(fast_prefix.name + "_pve_0.nii.gz"),
        "T1w_brain_fast_pve_1": fast_prefix.with_name(fast_prefix.name + "_pve_1.nii.gz"),
        "T1w_brain_fast_pve_2": fast_prefix.with_name(fast_prefix.name + "_pve_2.nii.gz"),
        "T1w_brain_fast_seg": fast_prefix.with_name(fast_prefix.name + "_seg.nii.gz"),
        "T1w_brain_fast_bias": fast_prefix.with_name(fast_prefix.name + "_bias.nii.gz"),
    }

    for key, src in list(outputs.items()):
        if key in {"T1w", "T1w_json"}:
            continue
        if not src.exists():
            print(f"[WARN] {row.subject_id}: missing output {src.name}")
            continue
        dst = subj_out / src.name
        shutil.copy2(src, dst)
        outputs[key] = dst

    # Write a small provenance file.
    prov = subj_out / "preprocess_provenance.txt"
    prov.write_text(
        "\n".join(
            [
                f"subject_id\t{row.subject_id}",
                f"source_nifti\t{row.nifti_path}",
                f"source_json\t{row.json_path}",
                f"series_description\t{row.series_description}",
                f"protocol_name\t{row.protocol_name}",
                f"score\t{row.score}",
                "steps\treorient -> bet -> fast",
            ]
        )
        + "\n"
    )

    outputs["provenance"] = prov
    return outputs


def main() -> int:
    parser = argparse.ArgumentParser(description="Preprocess selected T1-weighted MRI scans.")
    parser.add_argument(
        "--selected-manifest",
        type=Path,
        default=DEFAULT_SELECTED_MANIFEST,
        help="TSV manifest from refine_t1_selection.py",
    )
    parser.add_argument(
        "--qc-flags",
        type=Path,
        default=DEFAULT_QC_FLAGS,
        help="TSV file with GOOD/BAD visual QC labels",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=DEFAULT_OUT_ROOT,
        help="Output directory for preprocessed subjects",
    )
    parser.add_argument(
        "--work-root",
        type=Path,
        default=DEFAULT_WORK_ROOT,
        help="Working directory for intermediate files",
    )
    parser.add_argument(
        "--subject",
        nargs="*",
        default=None,
        help="Optional subject IDs to process, e.g. sub-009 sub-010",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing subject outputs",
    )
    args = parser.parse_args()

    check_tool("fslreorient2std")
    check_tool("bet")
    check_tool("fast")

    good_subjects = load_good_subjects(args.qc_flags)
    selected_rows = load_selected_manifest(args.selected_manifest)

    if args.subject:
        wanted = set(args.subject)
    else:
        wanted = None

    rows_to_process: List[SelectedRow] = []
    for row in selected_rows:
        if row.subject_id not in good_subjects:
            continue
        if wanted is not None and row.subject_id not in wanted:
            continue
        rows_to_process.append(row)

    if not rows_to_process:
        print("[WARN] No subjects matched QC GOOD + selected manifest.")
        return 0

    print(f"Loaded {len(rows_to_process)} GOOD subjects for preprocessing.")

    args.out_root.mkdir(parents=True, exist_ok=True)
    args.work_root.mkdir(parents=True, exist_ok=True)

    ok = 0
    failed = 0

    for row in rows_to_process:
        print(f"[RUN] {row.subject_id}: {row.nifti_path.name}")
        try:
            preprocess_subject(row, args.out_root, args.work_root, overwrite=args.overwrite)
            ok += 1
            print(f"[OK] {row.subject_id}")
        except subprocess.CalledProcessError as exc:
            failed += 1
            print(f"[FAIL] {row.subject_id}: command failed with exit code {exc.returncode}")
        except Exception as exc:
            failed += 1
            print(f"[FAIL] {row.subject_id}: {exc}")

    print(f"\nDone. OK={ok}, failed={failed}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())