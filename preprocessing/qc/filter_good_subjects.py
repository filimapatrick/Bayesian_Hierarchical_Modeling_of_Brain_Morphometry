#!/usr/bin/env python3

"""
Filter GOOD subjects based on QC flags and create a clean dataset.

Input:
- qc_t1_flags.tsv

Output:
- data/derivatives/analysis_ready_good/
"""

from pathlib import Path
import shutil

QC_FILE = Path("/Volumes/MyHDD/ABDN_DATA/qc_t1_flags.tsv")
SOURCE_DIR = Path("/Volumes/MyHDD/bayesian-brain-morphometry/data/derivatives/analysis_ready")
TARGET_DIR = Path("/Volumes/MyHDD/bayesian-brain-morphometry/data/derivatives/analysis_ready_good")


def load_good_subjects():
    good = []
    with QC_FILE.open() as f:
        for line in f:
            if line.startswith("subject_id"):
                continue
            subj, status = line.strip().split("\t")
            if status == "GOOD":
                good.append(subj)
    return good


def copy_subject(subj):
    src = SOURCE_DIR / subj
    dst = TARGET_DIR / subj

    if not src.exists():
        print(f"[WARN] Missing: {subj}")
        return

    if dst.exists():
        print(f"[SKIP] Already exists: {subj}")
        return

    shutil.copytree(src, dst)
    print(f"[OK] Copied: {subj}")


def main():
    TARGET_DIR.mkdir(parents=True, exist_ok=True)

    good_subjects = load_good_subjects()

    print(f"Found {len(good_subjects)} GOOD subjects\n")

    for subj in good_subjects:
        copy_subject(subj)

    print("\nDone.")


if __name__ == "__main__":
    main()