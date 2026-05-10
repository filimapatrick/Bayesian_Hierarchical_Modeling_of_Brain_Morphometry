#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


REPO_ROOT = Path("/Volumes/MyHDD/bayesian-brain-morphometry")
DEFAULT_DATASET = REPO_ROOT / "data" / "derivatives" / "features" / "updated_model_dataset.csv"
DEFAULT_OUT_DIR = REPO_ROOT / "results" / "pilot_model"


ID_COLS = {"subject_id", "diagnosis", "diagnosis_label", "source_file"}
PREFERRED_TARGET = "diagnosis"


def detect_feature_columns(df: pd.DataFrame) -> list[str]:
    numeric_cols = []
    for col in df.columns:
        if col in ID_COLS:
            continue
        if pd.api.types.is_numeric_dtype(df[col]):
            numeric_cols.append(col)
    return numeric_cols


def make_pipeline() -> Pipeline:
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "clf",
                LogisticRegression(
                    max_iter=5000,
                    class_weight="balanced",
                    random_state=42,
                ),
            ),
        ]
    )


def save_confusion_matrix(cm: np.ndarray, classes: list[str], out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(cm, interpolation="nearest")
    ax.figure.colorbar(im, ax=ax)
    ax.set(
        xticks=np.arange(len(classes)),
        yticks=np.arange(len(classes)),
        xticklabels=classes,
        yticklabels=classes,
        ylabel="True label",
        xlabel="Predicted label",
        title="Pilot model confusion matrix",
    )

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    thresh = cm.max() / 2.0 if cm.size else 0
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(
                j,
                i,
                int(cm[i, j]),
                ha="center",
                va="center",
                color="white" if cm[i, j] > thresh else "black",
            )

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_class_distribution(df: pd.DataFrame, out_path: Path) -> None:
    counts = df[PREFERRED_TARGET].value_counts().sort_index()

    fig, ax = plt.subplots(figsize=(7, 4))
    counts.plot(kind="bar", ax=ax)
    ax.set_title("Class distribution in merged dataset")
    ax.set_xlabel("Diagnosis")
    ax.set_ylabel("Number of subjects")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_feature_summary(df: pd.DataFrame, features: list[str], out_path: Path) -> None:
    summary = df.groupby(PREFERRED_TARGET)[features].mean(numeric_only=True)

    fig, ax = plt.subplots(figsize=(10, max(4, 0.5 * len(summary.index))))
    im = ax.imshow(summary.values, aspect="auto")
    ax.figure.colorbar(im, ax=ax)
    ax.set_xticks(np.arange(len(features)))
    ax.set_xticklabels(features, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(summary.index)))
    ax.set_yticklabels(summary.index)
    ax.set_title("Mean feature values by diagnosis")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Train a simple pilot model on the merged morphometry dataset.")
    parser.add_argument(
        "--dataset",
        type=Path,
        default=DEFAULT_DATASET,
        help="Merged feature CSV",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Directory for outputs",
    )
    args = parser.parse_args()

    if not args.dataset.exists():
        raise FileNotFoundError(f"Dataset not found: {args.dataset}")

    df = pd.read_csv(args.dataset)

    if PREFERRED_TARGET not in df.columns:
        raise ValueError(f"Expected a '{PREFERRED_TARGET}' column in {args.dataset}")

    # Keep only rows with a diagnosis label.
    df = df[df[PREFERRED_TARGET].notna()].copy()
    df[PREFERRED_TARGET] = df[PREFERRED_TARGET].astype(str).str.strip()

    # Drop duplicate subject/diagnosis rows if any.
    if "subject_id" in df.columns:
        df = df.drop_duplicates(subset=["subject_id", PREFERRED_TARGET])

    feature_cols = detect_feature_columns(df)
    if not feature_cols:
        raise ValueError("No numeric feature columns found in the dataset.")

    X = df[feature_cols].copy()
    y = df[PREFERRED_TARGET].copy()

    classes = sorted(y.unique())
    class_counts = y.value_counts()
    min_class_count = int(class_counts.min())

    args.out_dir.mkdir(parents=True, exist_ok=True)
    save_class_distribution(df, args.out_dir / "class_distribution.png")
    save_feature_summary(df, feature_cols, args.out_dir / "feature_means_by_class.png")

    model = make_pipeline()

    report_lines = []
    report_lines.append(f"Dataset: {args.dataset}")
    report_lines.append(f"Rows: {len(df)}")
    report_lines.append(f"Classes: {classes}")
    report_lines.append(f"Feature columns: {feature_cols}")
    report_lines.append("")
    report_lines.append("Class counts:")
    for k, v in class_counts.items():
        report_lines.append(f"  {k}: {v}")

    # Cross-validation only if each class has enough samples.
    if min_class_count >= 2 and len(classes) >= 2:
        n_splits = min(5, min_class_count)
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
        y_pred = cross_val_predict(model, X, y, cv=cv)

        acc = accuracy_score(y, y_pred)
        bal_acc = balanced_accuracy_score(y, y_pred)
        cm = confusion_matrix(y, y_pred, labels=classes)

        save_confusion_matrix(cm, classes, args.out_dir / "confusion_matrix.png")

        report_lines.append("")
        report_lines.append(f"Cross-validated accuracy: {acc:.4f}")
        report_lines.append(f"Cross-validated balanced accuracy: {bal_acc:.4f}")
        report_lines.append("")
        report_lines.append("Classification report:")
        report_lines.append(classification_report(y, y_pred, digits=4, zero_division=0))

        metrics = pd.DataFrame(
            {
                "metric": ["accuracy", "balanced_accuracy", "n_splits"],
                "value": [acc, bal_acc, n_splits],
            }
        )
        metrics.to_csv(args.out_dir / "metrics.csv", index=False)

        preds = pd.DataFrame(
            {
                "subject_id": df["subject_id"] if "subject_id" in df.columns else np.arange(len(df)),
                "true_label": y,
                "pred_label": y_pred,
            }
        )
        preds.to_csv(args.out_dir / "cv_predictions.csv", index=False)

        # Fit a final model on all data and save coefficients for inspection.
        model.fit(X, y)
        clf = model.named_steps["clf"]
        coef = getattr(clf, "coef_", None)
        if coef is not None:
            coef_df = pd.DataFrame(coef, columns=feature_cols)
            coef_df.insert(0, "class", clf.classes_)
            coef_df.to_csv(args.out_dir / "model_coefficients.csv", index=False)

    else:
        # Not enough samples for meaningful CV. Fit a descriptive model only.
        model.fit(X, y)
        clf = model.named_steps["clf"]
        preds = clf.predict(model.named_steps["scaler"].transform(model.named_steps["imputer"].transform(X)))

        report_lines.append("")
        report_lines.append(
            "Not enough subjects per class for stratified cross-validation. "
            "A descriptive fit was run on the full dataset only."
        )
        report_lines.append("")
        report_lines.append("Training-set classification report:")
        report_lines.append(classification_report(y, preds, digits=4, zero_division=0))

        cm = confusion_matrix(y, preds, labels=classes)
        save_confusion_matrix(cm, classes, args.out_dir / "confusion_matrix.png")

        metrics = pd.DataFrame(
            {
                "metric": ["n_classes", "min_class_count", "descriptive_fit_only"],
                "value": [len(classes), min_class_count, True],
            }
        )
        metrics.to_csv(args.out_dir / "metrics.csv", index=False)

        preds_df = pd.DataFrame(
            {
                "subject_id": df["subject_id"] if "subject_id" in df.columns else np.arange(len(df)),
                "true_label": y,
                "pred_label": preds,
            }
        )
        preds_df.to_csv(args.out_dir / "training_predictions.csv", index=False)

        coef = getattr(clf, "coef_", None)
        if coef is not None:
            coef_df = pd.DataFrame(coef, columns=feature_cols)
            coef_df.insert(0, "class", clf.classes_)
            coef_df.to_csv(args.out_dir / "model_coefficients.csv", index=False)

    report_path = args.out_dir / "pilot_model_report.txt"
    report_path.write_text("\n".join(report_lines) + "\n")

    print(f"Saved outputs to: {args.out_dir}")
    print(f"Report: {report_path}")
    print(f"Class counts:\n{class_counts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())