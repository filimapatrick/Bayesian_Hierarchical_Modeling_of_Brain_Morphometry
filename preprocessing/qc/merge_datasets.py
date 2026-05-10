#!/usr/bin/env python3
from __future__ import annotations

import pandas as pd
from pathlib import Path

DEMENTIA_CSV = Path("/Volumes/MyHDD/bayesian-brain-morphometry/data/derivatives/features/dementia/t1_features.csv")
EPILEPSY_CSV = Path("/Volumes/MyHDD/bayesian-brain-morphometry/data/derivatives/features/epilepsy/t1_features.csv")
HYDROCEPHALUS_CSV = Path("/Volumes/MyHDD/bayesian-brain-morphometry/data/derivatives/features/hydrocephalus/t1_features.csv")
CONTROL_CSV = Path("/Volumes/MyHDD/bayesian-brain-morphometry/data/derivatives/features/control/t1_features.csv")

OUT_CSV = Path("/Volumes/MyHDD/bayesian-brain-morphometry/data/derivatives/features/model_dataset.csv")


def load_and_label(csv_path: Path, label: str) -> pd.DataFrame:
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing file: {csv_path}")
    df = pd.read_csv(csv_path)
    df["diagnosis"] = label
    return df


def main() -> None:
    tables = [
        load_and_label(DEMENTIA_CSV, "DEMENTIA"),
        load_and_label(EPILEPSY_CSV, "EPILEPSY"),
        load_and_label(HYDROCEPHALUS_CSV, "HYDROCEPHALUS"),
        load_and_label(CONTROL_CSV, "CONTROL"),
    ]

    merged = pd.concat(tables, ignore_index=True)

    preferred = [c for c in ["subject_id", "diagnosis"] if c in merged.columns]
    other_cols = [c for c in merged.columns if c not in preferred]
    merged = merged[preferred + other_cols]

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUT_CSV, index=False)

    print(f"Saved merged dataset to: {OUT_CSV}")
    print("\nClass distribution:")
    print(merged["diagnosis"].value_counts())
    print("\nPreview:")
    print(merged.head())


if __name__ == "__main__":
    main()