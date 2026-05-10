#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path


def open_fsleyes(image_path: Path) -> None:
    subprocess.run(["fsleyes", str(image_path)], check=False)


def load_manifest(manifest_path: Path):
    with manifest_path.open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"subject_id", "nifti_path"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Manifest is missing columns: {sorted(missing)}")

        rows = []
        for row in reader:
            subject_id = row["subject_id"].strip()
            nifti_path = Path(row["nifti_path"].strip())
            rows.append((subject_id, nifti_path))
    return rows


def load_existing_flags(out_path: Path):
    flags = {}
    if out_path.exists():
        with out_path.open(newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if reader.fieldnames and {"subject_id", "status"} <= set(reader.fieldnames):
                for row in reader:
                    flags[row["subject_id"].strip()] = row["status"].strip().upper()
    return flags


def save_flags(flags: dict[str, str], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["subject_id", "status"])
        for subject_id, status in sorted(flags.items()):
            writer.writerow([subject_id, status])


def main() -> int:
    parser = argparse.ArgumentParser(description="Visual QC for selected T1 images.")
    parser.add_argument(
        "--manifest",
        type=Path,
        required=True,
        help="TSV manifest with subject_id and nifti_path columns",
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output TSV file for QC flags",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Skip subjects already labeled in the output TSV",
    )
    args = parser.parse_args()

    rows = load_manifest(args.manifest)
    if not rows:
        print(f"[WARN] No rows found in {args.manifest}")
        return 0

    flags = load_existing_flags(args.out) if args.resume else {}

    for subject_id, nifti_path in rows:
        if args.resume and subject_id in flags:
            continue

        if not nifti_path.exists():
            print(f"[WARN] {subject_id}: missing file {nifti_path}")
            continue

        print(f"\n=== {subject_id} ===")
        open_fsleyes(nifti_path)

        while True:
            label = input("Label (g=good, b=bad, s=skip, q=quit): ").strip().lower()
            if label in {"g", "b", "s", "q"}:
                break
            print("Invalid input")

        if label == "q":
            break
        if label == "g":
            flags[subject_id] = "GOOD"
        elif label == "b":
            flags[subject_id] = "BAD"
        elif label == "s":
            flags[subject_id] = "SKIP"

        save_flags(flags, args.out)

    save_flags(flags, args.out)
    print(f"\nSaved QC flags to: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
