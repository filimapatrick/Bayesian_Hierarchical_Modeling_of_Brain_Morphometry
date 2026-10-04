#!/usr/bin/env bash
# ==============================================================================
# MASTER REPRODUCTION SCRIPT: "What Survives the Blur?"
# Bayesian Hierarchical Modeling of Brain Morphometry Under Clinical Acquisition Heterogeneity
# ==============================================================================
# This orchestrator reproduces all feature extraction, MCMC posterior sampling,
# and experimental laboratory benchmarks (Experiments 1, 2, and 3) end-to-end.
#
# Usage:
#   bash run_study.sh [OPTIONS]
#
# Options:
#   --draws N        Number of MCMC draws per chain (default: 2000)
#   --tune N         Number of MCMC tuning steps per chain (default: 1000)
#   --chains N       Number of parallel Markov chains (default: 4)
#   --skip-features  Skip macro-feature extraction if macro_features.csv already exists
#   --help           Show this help message
# ==============================================================================

set -eo pipefail

# 1. Resolve Script & Project Directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# 2. Resolve Python Environment
if [[ -n "${VIRTUAL_ENV}" ]] && [[ -x "${VIRTUAL_ENV}/bin/python" ]]; then
    PYTHON_EXEC="${VIRTUAL_ENV}/bin/python"
elif [[ -x "/Users/patrick/.venvs/bayesian-brain-morphometry/bin/python" ]]; then
    PYTHON_EXEC="/Users/patrick/.venvs/bayesian-brain-morphometry/bin/python"
elif command -v python3 &>/dev/null; then
    PYTHON_EXEC="$(command -v python3)"
else
    echo "❌ Error: Python 3 executable not found!"
    exit 1
fi

echo "=============================================================================="
echo "🧠 MASTER REPRODUCTION PIPELINE: BAYESIAN BRAIN MORPHOMETRY"
echo "=============================================================================="
echo "Project Directory: ${SCRIPT_DIR}"
echo "Python Executable: ${PYTHON_EXEC}"
echo "Python Version:    $(${PYTHON_EXEC} --version)"
echo "Timestamp:         $(date)"
echo "=============================================================================="

# 3. Parse Command-Line Options
MCMC_DRAWS=2000
MCMC_TUNE=1000
MCMC_CHAINS=4
SKIP_FEATURES=false

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --draws)
            MCMC_DRAWS="$2"
            shift 2
            ;;
        --tune)
            MCMC_TUNE="$2"
            shift 2
            ;;
        --chains)
            MCMC_CHAINS="$2"
            shift 2
            ;;
        --skip-features)
            SKIP_FEATURES=true
            shift 1
            ;;
        --help|-h)
            echo "Usage: bash run_study.sh [--draws N] [--tune N] [--chains N] [--skip-features]"
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: bash run_study.sh [--draws N] [--tune N] [--chains N] [--skip-features]"
            exit 1
            ;;
    esac
done

# 4. Ensure Directory Hierarchy
mkdir -p results/tables results/figures results/traces

# 5. Verify Core Python Dependencies
echo ""
echo "🔍 Verifying Python Dependencies..."
${PYTHON_EXEC} -c '
import numpy
import scipy
import scipy.signal
import scipy.signal.windows
if not hasattr(scipy.signal, "gaussian"):
    scipy.signal.gaussian = scipy.signal.windows.gaussian
import pandas
import nibabel
import matplotlib
import seaborn
import pymc
import arviz
print("  ✓ All required scientific libraries are installed:")
print(f"    - numpy: {numpy.__version__}")
print(f"    - scipy: {scipy.__version__}")
print(f"    - pandas: {pandas.__version__}")
print(f"    - nibabel: {nibabel.__version__}")
print(f"    - matplotlib: {matplotlib.__version__}")
print(f"    - seaborn: {seaborn.__version__}")
print(f"    - pymc: {pymc.__version__}")
print(f"    - arviz: {arviz.__version__}")
' || {
    echo "❌ Missing scientific dependencies. Run: pip install numpy scipy pandas nibabel matplotlib seaborn pymc arviz"
    exit 1
}

START_TIME=$(date +%s)

# ==============================================================================
# STEP 1: MACRO-MORPHOMETRIC FEATURE EXTRACTION
# ==============================================================================
echo ""
echo "=============================================================================="
echo "▶ STEP 1/5: Macro-Morphometric Feature Extraction (N=218 Cohort)"
echo "=============================================================================="
if [[ "${SKIP_FEATURES}" == true ]] && [[ -f "results/tables/macro_features.csv" ]]; then
    echo "⏩ --skip-features specified and results/tables/macro_features.csv exists. Skipping."
else
    ${PYTHON_EXEC} features/extract_features.py \
        --bids_dir "${SCRIPT_DIR}/data/bids" \
        --output_csv "${SCRIPT_DIR}/results/tables/macro_features.csv"
fi

# ==============================================================================
# STEP 2: HIERARCHICAL BAYESIAN MCMC POSTERIOR INFERENCE (PyMC 5.12 NUTS)
# ==============================================================================
echo ""
echo "=============================================================================="
echo "▶ STEP 2/5: Hierarchical Bayesian MCMC Sampling across Biomarkers (PyMC NUTS)"
echo "  Chains: ${MCMC_CHAINS} | Draws: ${MCMC_DRAWS} | Tune: ${MCMC_TUNE}"
echo "=============================================================================="

# 2a. Parenchymal Envelope Fraction (PEF)
echo "--- 2a. Sampling Posterior for PEF (Parenchymal Envelope Fraction) ---"
${PYTHON_EXEC} modeling/inference.py \
    --features_csv results/tables/macro_features.csv \
    --target_metric pef \
    --draws "${MCMC_DRAWS}" \
    --tune "${MCMC_TUNE}" \
    --chains "${MCMC_CHAINS}"

# 2b. Ventricle-to-Brain Ratio (VBR)
echo "--- 2b. Sampling Posterior for VBR (Ventricle-to-Brain Ratio) ---"
${PYTHON_EXEC} modeling/inference.py \
    --features_csv results/tables/macro_features.csv \
    --target_metric vbr \
    --draws "${MCMC_DRAWS}" \
    --tune "${MCMC_TUNE}" \
    --chains "${MCMC_CHAINS}"

# 2c. Evans' Index
echo "--- 2c. Sampling Posterior for Evans' Index (Flagship Hydrocephalus Phenotype) ---"
${PYTHON_EXEC} modeling/inference.py \
    --features_csv results/tables/macro_features.csv \
    --target_metric evans_index \
    --draws "${MCMC_DRAWS}" \
    --tune "${MCMC_TUNE}" \
    --chains "${MCMC_CHAINS}"

# ==============================================================================
# STEP 3: EXPERIMENT 1 - CONTROLLED SYNTHETIC DEGRADATION LAB
# ==============================================================================
echo ""
echo "=============================================================================="
echo "▶ STEP 3/5: Experiment 1 - Controlled Synthetic Degradation Lab"
echo "  Simulating slice blur: 1.0mm -> 3.0mm -> 4.0mm -> 5.0mm -> 6.0mm"
echo "=============================================================================="
${PYTHON_EXEC} features/synthetic_degradation.py \
    --bids_dir "${SCRIPT_DIR}/data/bids" \
    --output_csv "${SCRIPT_DIR}/results/tables/experiment1_synthetic_degradation.csv" \
    --summary_csv "${SCRIPT_DIR}/results/tables/experiment1_degradation_summary.csv" \
    --plot_path "${SCRIPT_DIR}/results/figures/figure1_synthetic_degradation_curves.png"

# ==============================================================================
# STEP 4: EXPERIMENT 2 - CLINICAL FEASIBILITY & FAILURE BOUNDARIES
# ==============================================================================
echo ""
echo "=============================================================================="
echo "▶ STEP 4/5: Experiment 2 - Real-World Feasibility & Logistic Failure Model"
echo "=============================================================================="
${PYTHON_EXEC} features/clinical_failure_boundaries.py \
    --features_csv "${SCRIPT_DIR}/results/tables/macro_features.csv" \
    --output_csv "${SCRIPT_DIR}/results/tables/experiment2_failure_boundaries.csv" \
    --plot_path "${SCRIPT_DIR}/results/figures/figure2_clinical_feasibility_boundaries.png"

# ==============================================================================
# STEP 5: EXPERIMENT 3 - MODEL BENCHMARKING & CONFOUNDING SENSITIVITY
# ==============================================================================
echo ""
echo "=============================================================================="
echo "▶ STEP 5/6: Experiment 3 - Model Benchmarks & Sensitivity Analysis"
echo "=============================================================================="
${PYTHON_EXEC} modeling/sensitivity_analysis.py \
    --features_csv "${SCRIPT_DIR}/results/tables/macro_features.csv" \
    --benchmarks_csv "${SCRIPT_DIR}/results/tables/experiment3_model_benchmarks.csv" \
    --sensitivity_csv "${SCRIPT_DIR}/results/tables/experiment3_sensitivity_analysis.csv" \
    --fig3_path "${SCRIPT_DIR}/results/figures/figure3_posterior_shrinkage_forest.png" \
    --fig4_path "${SCRIPT_DIR}/results/figures/figure4_variance_partitioning_sensitivity.png"

# ==============================================================================
# STEP 6: EXPERIMENT 3 - POSTERIOR PREDICTIVE CHECKS & BOUNDARY EVALUATION
# ==============================================================================
echo ""
echo "=============================================================================="
echo "▶ STEP 6/6: Experiment 3 - Posterior Predictive Checks (PPC) & Model Assessment"
echo "=============================================================================="
${PYTHON_EXEC} modeling/posterior_predictive_check.py \
    --features_csv "${SCRIPT_DIR}/results/tables/macro_features.csv" \
    --traces_dir "${SCRIPT_DIR}/results/traces" \
    --output_csv "${SCRIPT_DIR}/results/tables/posterior_predictive_summary.csv" \
    --output_plot "${SCRIPT_DIR}/results/figures/figure5_posterior_predictive_checks.png"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# ==============================================================================
# ARTIFACT VERIFICATION & SCORECARD
# ==============================================================================
echo ""
echo "=============================================================================="
echo "🎉 REPRODUCTION COMPLETE: VERIFYING STUDY ARTIFACTS"
echo "   Total Elapsed Time: ${DURATION} seconds"
echo "=============================================================================="

${PYTHON_EXEC} -c '
import os
from pathlib import Path

artifacts = [
    # Tables
    ("results/tables/macro_features.csv", "Table: Macro-Morphometric Cohort Features"),
    ("results/tables/freesurfer_qc.csv", "Table: FreeSurfer QC 218 Subject Log"),
    ("results/tables/pipeline_attrition_comparison.csv", "Table: Exp 2 FreeSurfer Attrition Comparison"),
    ("results/tables/posterior_summary_pef.csv", "Table: PEF MCMC Posterior Summary"),
    ("results/tables/posterior_summary_bpf.csv", "Table: BPF (Legacy Mirror) MCMC Posterior Summary"),
    ("results/tables/posterior_summary_vbr.csv", "Table: VBR MCMC Posterior Summary"),
    ("results/tables/posterior_summary_evans_index.csv", "Table: Evans Index MCMC Posterior Summary"),
    ("results/tables/experiment1_degradation_summary.csv", "Table: Exp 1 Synthetic Degradation Summary"),
    ("results/tables/experiment2_failure_boundaries.csv", "Table: Exp 2 Logistic Failure Model"),
    ("results/tables/experiment3_model_benchmarks.csv", "Table: Exp 3 Model Benchmarks"),
    ("results/tables/experiment3_sensitivity_analysis.csv", "Table: Exp 3 Confounding Sensitivity Tests"),
    ("results/tables/posterior_predictive_summary.csv", "Table: Exp 3 PPC Coverage Summary"),
    # Traces
    ("results/traces/mcmc_traces_pef.npz", "Traces: PEF MCMC Markov Chains"),
    ("results/traces/mcmc_traces_bpf.npz", "Traces: BPF (Legacy Mirror) MCMC Markov Chains"),
    ("results/traces/mcmc_traces_vbr.npz", "Traces: VBR MCMC Markov Chains"),
    ("results/traces/mcmc_traces_evans_index.npz", "Traces: Evans Index MCMC Markov Chains"),
    # Figures
    ("results/figures/figure1_synthetic_degradation_curves.png", "Figure 1: Exp 1 Synthetic Degradation Curves"),
    ("results/figures/figure2_clinical_feasibility_boundaries.png", "Figure 2: Exp 2 Clinical Feasibility Boundaries"),
    ("results/figures/figure3_posterior_shrinkage_forest.png", "Figure 3: Exp 3 Posterior Shrinkage Forest Plot"),
    ("results/figures/figure4_variance_partitioning_sensitivity.png", "Figure 4: Exp 3 Variance Partitioning & Sensitivity"),
    ("results/figures/figure5_posterior_predictive_checks.png", "Figure 5: Exp 3 Posterior Predictive Distributions"),
]

all_ok = True
print(f"Status | Size (KB) | Artifact Path & Description")
print("-" * 78)
for rel_path, desc in artifacts:
    p = Path(rel_path)
    if p.exists() and p.stat().st_size > 0:
        size_kb = p.stat().st_size / 1024.0
        print(f"  ✓    | {size_kb:8.1f}  | {rel_path:<40} ({desc})")
    else:
        print(f"  ❌   |    MISSING | {rel_path:<40} ({desc})")
        all_ok = False

print("-" * 78)
if all_ok:
    print("✨ ALL 21 STUDY ARTIFACTS VERIFIED SUCCESSFULLY!")
else:
    print("⚠️ Some artifacts are missing or zero-sized. Review step logs above.")
'

echo "=============================================================================="
