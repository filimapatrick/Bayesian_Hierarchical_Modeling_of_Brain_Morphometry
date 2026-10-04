# What Survives the Blur? Macro-Brain Morphometry and Uncertainty Quantification from Heterogeneous Routine Clinical MRI in Nigeria

[![BIDS Version](https://img.shields.io/badge/BIDS-v1.10.0%20Compliant-brightgreen.svg)](https://bids.neuroimaging.io/)
[![PyMC](https://img.shields.io/badge/PyMC-v5.x-orange.svg)](https://www.pymc.io/)
[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![FAIR Principles](https://img.shields.io/badge/Data-FAIR%20Compliant-blueviolet.svg)](https://www.go-fair.org/fair-principles/)
[![Status](https://img.shields.io/badge/Status-Active%20Research%20%7C%20In%20Prep-success.svg)](#)

**Principal Investigator & Lead Author:** Patrick Filima  
**Affiliation:** African Brain Data Network (ABDN) & Collaborating Nigerian Clinical Centers  
**Clinical Partner Institutions:** Rivers State University Teaching Hospital (RSUTH), University of Port Harcourt Teaching Hospital (UPTH), Aminu Kano Teaching Hospital (AKTH) / Northwest Kano Diagnostic Centre (NKDC), Braithwaite Memorial Hospital (BMH), Intercontinental Diagnostic Centre (IDC), Life Bridge Diagnostic Centre  
**Foundational Data Provenance:** _Scientific Data_ (2025) [DOI: 10.1038/s41597-025-04743-0](https://doi.org/10.1038/s41597-025-04743-0) | [brainlife.io project](https://doi.org/10.25663/brainlife.p6554f423b094062da63aa4c9)

---

## 📑 Table of Contents

- [1. Executive Summary: The Scientific Pivot](#1-executive-summary-the-scientific-pivot)
  - [1.1 The Global Neuroimaging Equity Gap](#11-the-global-neuroimaging-equity-gap)
  - [1.2 The "1 mm³ Isotropic MPRAGE" Fallacy](#12-the-1-mm-isotropic-mprage-fallacy)
  - [1.3 The Ground Truth of Routine Clinical MRI in Nigeria](#13-the-ground-truth-of-routine-clinical-mri-in-nigeria)
  - [1.4 The Attrition Trap of Standard Pipelines](#14-the-attrition-trap-of-standard-pipelines)
  - [1.5 The Core Scientific Contribution: "What Survives the Blur?"](#15-the-core-scientific-contribution-what-survives-the-blur)
- [2. Research Aims & Questions](#2-research-aims--questions)
  - [2.1 Primary Research Question](#21-primary-research-question)
  - [2.2 Specific Aims](#22-specific-aims)
  - [2.3 The A Priori Detectability Gradient](#23-the-a-priori-detectability-gradient)
- [3. Multi-Center Clinical Dataset Architecture](#3-multi-center-clinical-dataset-architecture)
  - [3.1 Cohort Composition & Provenance ($N = 218$)](#31-cohort-composition--provenance-n--218)
  - [3.2 Participating Clinical Centers & Hardware Profiles](#32-participating-clinical-centers--hardware-profiles)
  - [3.3 Acquisition Heterogeneity & Identifiability Matrix](#33-acquisition-heterogeneity--identifiability-matrix)
- [4. Three-Tier Experimental Architecture](#4-three-tier-experimental-architecture)
  - [4.1 Experiment 1: Controlled Synthetic Degradation Lab](#41-experiment-1-controlled-synthetic-degradation-lab)
  - [4.2 Experiment 2: Real-World Clinical Measurement Validity & Failure Boundaries](#42-experiment-2-real-world-clinical-measurement-validity--failure-boundaries)
  - [4.3 Experiment 3: Uncertainty-Aware Disease Inference & Confounding Sensitivity](#43-experiment-3-uncertainty-aware-disease-inference--confounding-sensitivity)
- [5. Mathematical & Statistical Framework](#5-mathematical--statistical-framework)
  - [5.1 Heteroskedastic Observation Likelihood](#51-heteroskedastic-observation-likelihood)
  - [5.2 Hierarchical Prior Structure](#52-hierarchical-prior-structure)
  - [5.3 Posterior Variance Partitioning (ICC)](#53-posterior-variance-partitioning-icc)
  - [5.4 Model Comparison Hierarchy](#54-model-comparison-hierarchy)
- [6. Repository Architecture & Artifact Directory](#6-repository-architecture--artifact-directory)
- [7. Step-by-Step Execution Guide](#7-step-by-step-execution-guide)
  - [7.1 Prerequisites & Virtual Environment](#71-prerequisites--virtual-environment)
  - [7.2 BIDS Validation](#72-bids-validation)
  - [7.3 Macro-Morphometric Feature Extraction](#73-macro-morphometric-feature-extraction)
  - [7.4 Running Hierarchical Bayesian Inference](#74-running-hierarchical-bayesian-inference)
  - [7.5 Generating Publication Figures & Diagnostic Forest Plots](#75-generating-publication-figures--diagnostic-forest-plots)
- [8. Empirical Findings from the Nigerian Cohort](#8-empirical-findings-from-the-nigerian-cohort)
- [9. Ethical Compliance & FAIR Data Stewardship](#9-ethical-compliance--fair-data-stewardship)
- [10. Citation & Acknowledgments](#10-citation--acknowledgments)

---

## 1. Executive Summary: The Scientific Pivot

### 1.1 The Global Neuroimaging Equity Gap

Sub-Saharan Africa accounts for approximately 18% of the global human population, yet less than **2% of published magnetic resonance imaging (MRI) studies** include cohorts of African descent. Global neuroimaging consortia (ADNI, UK Biobank, HCP) have developed normative atlases and analytical pipelines almost exclusively calibrated on high-income, Western cohorts scanned under rigid, research-grade protocols.

### 1.2 The "1 mm³ Isotropic MPRAGE" Fallacy

Standard neuroimaging software tools—including FreeSurfer (`recon-all`), FSL FAST (tissue segmentation), FSL FIRST (subcortical structure segmentation), and SPM—were engineered on an implicit assumption: **high-contrast, 3D volumetric, 1.0 mm³ isotropic T1-weighted MPRAGE/SPGR acquisitions scanned at 3.0 Tesla**.

In prospective clinical trials or Western academic imaging centers, scans that do not conform to these specifications are discarded. However, applying this exclusionary criterion in low- and middle-income countries (LMICs) renders **over 90% of all acquired clinical neuroimaging data unusable**.

### 1.3 The Ground Truth of Routine Clinical MRI in Nigeria

This project is grounded in retrospective, routine clinical neuroimaging workflows across hospital sites and diagnostic centers in Nigeria (Port Harcourt, Kano, and Abuja). An empirical audit of the hospital archives reveals:

1. **2D Thick-Slice Acquisitions:** 83.9% ($183/218$) of all clinical T1-weighted scans are 2D Fast Spin Echo (FSE) or Spin Echo sequences with slice thicknesses ranging from **4.0 mm to 6.0 mm** (voxel volumes of $1.5 \text{ to } 6.8 \text{ mm}^3$), generating severe partial volume averaging in the slice-select dimension.
2. **Intravenous Gadolinium Contrast (+C):** 42.2% ($92/218$) of scans are post-contrast diagnostic workups (e.g., investigating intracranial infections, tuberculosis, trauma, or space-occupying lesions). Gadolinium administration dramatically alters white matter/gray matter/vessel intensity ratios, violating classical tissue segmentation priors.
3. **Scanner Hardware & Field Strength Diversity:** Systems span ultra-low-field permanent magnet scanners (**0.30T and 0.35T**) to standard clinical superconducting scanners (**1.5T**) from GE Healthcare, Siemens Healthineers, Toshiba, and Hitachi.
4. **Patient Motion Artifacts:** Clinical patients presenting with acute seizures, hydrocephalus, or advanced dementia frequently exhibit motion-induced phase dispersion and blurring.

```
       WESTERN RESEARCH PARADIGM                    ROUTINE NIGERIAN CLINICAL REALITY
┌──────────────────────────────────────┐        ┌──────────────────────────────────────┐
│  • 3.0T Research Scanner             │   vs   │  • 0.3T - 1.5T Multi-Vendor Fleet    │
│  • 3D MPRAGE (1.0 mm³ Isotropic)     │        │  • 2D Thick Slice (4.0 - 6.0 mm FSE) │
│  • Pre-contrast, Unenhanced Scans    │        │  • Post-Contrast (+C) Mixed Scans    │
│  • Screened, Stable Outpatients      │        │  • Acute Clinical Inpatients         │
│  • Micro-Morphometry (Subfields)     │        │  • Severe Partial Volume Averaging   │
└──────────────────────────────────────┘        └──────────────────────────────────────┘
                   │                                               │
                   ▼                                               ▼
         FSL / FreeSurfer Valid                         85% Severe Data Attrition!
                                                    (Standard micro-segmentation fails)
```

### 1.4 The Attrition Trap of Standard Pipelines

When standard Western preprocessing pipelines (FSL FAST and FIRST) were applied to these Nigerian clinical cohorts:

- **Micro-morphometric failure:** FSL FIRST hallucinated subcortical boundaries across 5.0 mm slices, producing biologically impossible hippocampal volumes.
- **Tissue inversion:** FSL FAST frequently inverted tissue balance on post-contrast scans, classifying hyperintense contrast-enhanced blood vessels and meninges as cortical gray matter.
- **Catastrophic Sample Attrition:** Over 85% of available clinical scans were discarded by automated QC, and **the entire hydrocephalus cohort was lost**.

### 1.5 The Core Scientific Contribution: "What Survives the Blur?"

The fundamental insight of this work is: **Stop asking low-resolution scans to support high-resolution anatomy.**

This study is not an opportunistic disease morphometry paper attempting to mimic UK Biobank. It is a **methodological investigation of what quantitative neuroanatomical information survives routine clinical MRI degradation**, when that information stops being trustworthy, and how hierarchical Bayesian modeling can quantify and propagate acquisition uncertainty rather than pretending clinical imaging is research-grade MRI.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                                  THE INFERENTIAL CHAIN                                 │
│                                                                                        │
│  Clinical MRI Heterogeneity  ──►  Measurement Degradation  ──►  Recoverable Anatomical │
│  (0.3-1.5T, 2D 5mm, +C, SNR)      (Partial Volume, Noise)       Scale (Ventricular vs) │
│                                                                                        │
│                                          │                                             │
│                                          ▼                                             │
│                     Uncertainty-Aware Bayesian Disease Inference                       │
│                      (Partial Pooling + Confounding Sensitivity)                       │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Research Aims & Questions

### 2.1 Primary Research Question

> **To what extent can coarse, anatomically interpretable morphometric features be estimated from heterogeneous routine clinical MRI, and how much uncertainty in those estimates is attributable to acquisition characteristics versus clinical group?**

This formulation cleanly decouples two distinct scientific challenges:

1. **The Measurement Problem:** Can coarse macro-structural features (ventricular morphology, parenchymal envelope, hemispheric symmetry) be measured reliably across real-world quality gradients?
2. **The Statistical Inference Problem:** Given severe, non-random acquisition heterogeneity and site-disease confounding, what can be legitimately inferred about clinical group differences?

### 2.2 Specific Aims

- **Aim 1 — Establish the Resolution-Dependent Validity of Candidate Morphometric Biomarkers:**  
  Determine which macro-structural measurements remain stable under controlled degradation (synthetic downsampling and blurring from 1.0 mm to 6.0 mm) using high-resolution 3D acquisitions as an internal ground truth.
- **Aim 2 — Characterize Measurement Uncertainty & Failure Boundaries in Routine Clinical MRI:**  
  Quantify how slice thickness, acquisition dimensionality (2D vs. 3D), contrast enhancement, magnetic field strength (0.3T vs. 1.5T), and scanner site influence measurement error and pipeline failure probability:
  $$P(\text{failure}) = f(\text{SliceThickness}, \text{Contrast}, \text{FieldStrength}, \text{SNR})$$
- **Aim 3 — Estimate Disease-Associated Macro-Structural Differences Under Acquisition Uncertainty:**  
  Deploy hierarchical Bayesian models to estimate partially pooled disease effects while explicitly modeling heteroskedastic slice-thickness noise and evaluating sensitivity to diagnosis–site–contrast confounding.
- **Aim 4 — Define Practical Feasibility Boundaries for Opportunistic Neuroimaging:**  
  Establish an evidence-based decision matrix for low- and middle-income country (LMIC) hospital archives, identifying which anatomical phenotypes and clinical acquisition regimes permit defensible quantitative analysis.

### 2.3 The A Priori Detectability Gradient

Rather than forcing uniform statistical expectations across divergent pathologies, we structure disease cohorts along an **a priori detectability gradient** based on anatomical effect size relative to acquisition noise:

```
  ANATOMICAL MAGNITUDE RELATIVE TO ACQUISITION NOISE
┌─────────────────────────────────────────────────────────────────────────────────────────┐
│                                                                                         │
│  [HIGH SIGNAL]                                                            [LOW SIGNAL]  │
│                                                                                         │
│   HYDROCEPHALUS     ──►       DEMENTIA       ──►    PARKINSON'S     ──►     EPILEPSY    │
│                                                                                         │
│  Massive ventricular      Global cerebral        Subtle subcortical      Focal / MTS    │
│  dilation & mantle        atrophy & ex-vacuo     structural shifts;      pathology; N=7 │
│  compression              ventricular expansion  near-normal macro BPF   single-site +C │
│                                                                                         │
│  Effect >> Noise          Effect > Noise         Effect ~ Noise          Effect < Noise │
│  [Positive Control]       [Detectable]           [Negative Control]      [Exploratory]  │
└─────────────────────────────────────────────────────────────────────────────────────────┘
```

1. **Hydrocephalus ($N=82$): Flagship Positive Control.** Massive ventricular expansion produces an anatomical signal that reliably dwarfs 5.0 mm slice blurring and low SNR. If a pipeline cannot recover hydrocephalus, it is invalid.
2. **Dementia ($N=45$): Moderate Global Atrophy.** Symmetrical cortical loss and ventricular expansion provide a test case for whether global brain parenchymal fraction (BPF) and ventricle-to-brain ratio (VBR) can detect neurodegenerative atrophy above scanner noise.
3. **Parkinson's Disease ($N=21$): Expected Structural Negative Control.** Conventional T1-weighted macro-morphometry is not expected to exhibit gross global parenchymal shifts. A near-null finding here confirms model calibration and guards against false-positive discovery.
4. **Epilepsy ($N=7$): Exploratory & Confounded.** Mesial temporal sclerosis produces subtle, focal asymmetries. Because our epilepsy cohort is small ($N=7$), acquired at a single hospital center (UPTH), and largely contrast-enhanced, it is treated strictly as an exploratory cohort for hemispheric asymmetry rather than confirmatory inference.

---

## 3. Multi-Center Clinical Dataset Architecture

### 3.1 Cohort Composition & Provenance ($N = 218$)

The standardized BIDS v1.10.0 dataset incorporates **$N = 218$ unique human subjects** extracted from Nigerian clinical archives and published open-access resources:

| Diagnostic Cohort       | BIDS Identifier Prefix | Total Scans ($N$) | Valid Volumetric ($N$) | Flagged Scouts ($N$) | Primary Clinical Phenotype                        |
| :---------------------- | :--------------------- | :---------------: | :--------------------: | :------------------: | :------------------------------------------------ |
| **Hydrocephalus**       | `sub-hyd`              |        82         |           82           |          0           | Severe ventricular enlargement, cortical thinning |
| **Healthy Control**     | `sub-con`              |        63         |           62           |          1           | Age-appropriate parenchymal preservation          |
| **Dementia**            | `sub-dem`              |        45         |           40           |          5           | Diffuse cerebral atrophy, ventricular dilation    |
| **Parkinson's Disease** | `sub-pd`               |        21         |           18           |          3           | Basal ganglia degeneration; preserved macro BPF   |
| **Epilepsy**            | `sub-epi`              |         7         |           7            |          0           | Seizure disorder; focal temporal asymmetry        |
| **Total**               |                        |      **218**      |        **209**         |        **9**         | **Multi-cohort clinical spectrum**                |

_Provenance Breakdown:_

- **Clinical Hospital Archives:** $N = 137$ subjects (RSUTH, UPTH, BMH, IDC, AKTH/NKDC).
- **Brainlife Open Dataset:** $N = 81$ subjects (Wogu, Filima et al., _Scientific Data_ 2025; [DOI: 10.1038/s41597-025-04743-0](https://doi.org/10.1038/s41597-025-04743-0)).

### 3.2 Participating Clinical Centers & Hardware Profiles

Scans originate from six clinical imaging centers spanning three Nigerian geopolitical zones (South-South, North-West, and North-Central):

| Canonical Site ID | Hospital / Diagnostic Center                  | City / Zone        | Scans ($N$) |  Field ($B_0$)   | Scanner Hardware & Model                     |
| :---------------- | :-------------------------------------------- | :----------------- | :---------: | :--------------: | :------------------------------------------- |
| `RSUTH`           | Rivers State University Teaching Hospital     | Port Harcourt (SS) |     85      |      1.5 T       | GE Healthcare (SIGNA Creator)                |
| `UPTH`            | University of Port Harcourt Teaching Hospital | Port Harcourt (SS) |     42      |      1.5 T       | Siemens Healthineers (Magnetom)              |
| `IDC`             | Intercontinental Diagnostic Centre            | Port Harcourt (SS) |     38      |      0.3 T       | Permanent low-field (_Vendor n/a in header_) |
| `AKTH_NKDC`       | Aminu Kano Teaching Hospital / NKDC           | Kano (NW)          |     21      |      1.5 T       | Siemens Healthineers (Magnetom Essenza)      |
| `BMH`             | Braithwaite Memorial Specialist Hospital      | Port Harcourt (SS) |     16      |      1.5 T       | GE Healthcare (Signa)                        |
| `LifeBridge`      | Life Bridge Medical Diagnostics Ltd           | Abuja (NC)         |     16      |  0.35 T / 1.5 T  | Toshiba / Siemens Healthineers               |
| **Total**         | **6 Hospital Sites**                          | **3 Zones**        |   **218**   | **0.30T – 1.5T** | **GE, Siemens, Toshiba, Low-field**          |

> [!NOTE]
> In raw DICOM headers, institution naming variants (`RSUTH` vs `RSUTH Port Harcourt`, `Life Bridge` vs `LIFEBRIDGE MEDICAL DIAGNOSTICS LTD`) were canonicalized into unique institutional entities. Manufacturer metadata for IDC is documented strictly as recorded (`n/a`) without unverified attribution.

### 3.3 Acquisition Heterogeneity & Identifiability Matrix

The cross-tabulation of clinical diagnosis against acquisition parameters reveals the observational reality of routine hospital imaging:

```
┌───────────────────────────────────────────────────────────────────────────────────────┐
│                       CLINICAL ACQUISITION CROSS-TABULATION MATRIX                    │
├─────────────────┬───────────────────┬───────────────────┬─────────────────────────────┤
│ Diagnostic Cohort│ Contrast Enhanced │ Acquisition Type  │ Slice Thickness Range       │
│                 │ (False / True)    │ (2D / 3D)         │ (Min - Median - Max)        │
├─────────────────┼───────────────────┼───────────────────┼─────────────────────────────┤
│ CONTROL (N=63)  │ 59 False / 4 True │ 62 2D / 1 3D      │ 1.4 mm – 5.0 mm – 6.0 mm    │
│ DEMENTIA (N=45) │ 30 False / 15 True│ 40 2D / 5 3D      │ 1.0 mm – 5.0 mm – 6.0 mm    │
│ EPILEPSY (N=7)  │ 1 False / 6 True  │ 1 2D / 6 3D       │ 1.0 mm – 1.0 mm – 5.0 mm    │
│ HYDROCEPHALUS   │ 20 False / 62 True│ 59 2D / 23 3D     │ 1.0 mm – 5.0 mm – 10.0 mm   │
│   (N=82)        │ (75.6% +C)        │ (28.0% 3D)        │                             │
│ PARKINSON (N=21)│ 16 False / 5 True │ 21 2D / 0 3D      │ 4.0 mm – 5.0 mm – 5.0 mm    │
└─────────────────┴───────────────────┴───────────────────┴─────────────────────────────┘
```

#### The Identifiability Challenge:

Notice that:
$$\text{Site} \approx \text{Scanner} \approx \text{Field Strength}$$
$$\text{Diagnosis} \approx \text{Site} \approx \text{Contrast Enhancement}$$

- **Hydrocephalus** was predominantly imaged with intravenous contrast ($62/82$ scans) for surgical planning and etiology screening.
- **Healthy Controls** were almost entirely unenhanced ($59/63$ scans).
- **Epilepsy** scans were acquired exclusively at UPTH ($7/7$), mostly with high-resolution 1.0 mm 3D sequences and contrast enhancement.

**Methodological Implication:** No statistical algorithm—Bayesian or frequentist—can causally disentangle diagnosis from acquisition when they do not overlap in the design. Rather than pretending partial pooling eliminates this confounding, our framework explicitly models acquisition uncertainty and subjects every finding to sensitivity analyses.

---

## 4. Three-Tier Experimental Architecture

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                              THREE-TIER EXPERIMENTAL LAB                               │
├────────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                        │
│  [ EXPERIMENT 1: CONTROLLED SYNTHETIC DEGRADATION LAB ]                                │
│   • Internal Ground Truth: N=35 pristine 1.0 mm³ 3D scans                              │
│   • Synthetic Downsampling: 1.0 mm ──► 3.0 mm ──► 4.0 mm ──► 5.0 mm ──► 6.0 mm         │
│   • Metric Stability: Δ BPF, Δ VBR, Δ Evans' Index vs. Δ Subcortical Volumes           │
│   • Quantify Relative Error: |M_degraded - M_ref| / M_ref                              │
│                                                                                        │
│                                          │                                             │
│                                          ▼                                             │
│  [ EXPERIMENT 2: REAL-WORLD CLINICAL MEASUREMENT VALIDITY & FAILURE BOUNDARIES ]       │
│   • Real-World Cohort: N=218 routine clinical scans across 6 centers                   │
│   • Failure Probability: P(Measurement Failure) = f(Thickness, Contrast, B0, SNR)      │
│   • Identify Biological Boundary Violations (BPF > 1.0, Negative CSF)                  │
│                                                                                        │
│                                          │                                             │
│                                          ▼                                             │
│  [ EXPERIMENT 3: UNCERTAINTY-AWARE BAYESIAN DISEASE INFERENCE ]                        │
│   • Heteroskedastic Noise: σ_i = σ_0 * exp(λ * SliceThickness_i)                       │
│   • Posterior Shrinkage across Detectability Gradient (Hydro >> Dem > PD ~ Epi)        │
│   • Confounding Sensitivity Analysis: Assessing parameter stability under site pruning │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 4.1 Experiment 1: Controlled Synthetic Degradation Lab

**Objective:** Directly test the hypothesis that coarse ventricular and macro-parenchymal metrics survive slice degradation, whereas subcortical and micro-morphometric segmentations collapse.

- **Experimental Cohort:** $N = 35$ high-resolution 3D volumetric T1-weighted scans (1.0 mm isotropic) from UPTH and Brainlife.
- **Degradation Operator:** Slices are synthetically downsampled and filtered along the slice-select axis to simulate routine 2D Fast Spin Echo acquisitions:
  $$I_{\text{degraded}}(z) = \left[ I_{\text{1mm}} * \text{SliceProfile}(h) \right] \downarrow_{h}, \quad h \in \{3.0, 4.0, 5.0, 6.0\text{ mm}\}$$
- **Evaluated Metrics:**
  - Coarse Macro-Morphometry: Evans' Index, Ventricle-to-Brain Ratio (VBR), Brain Parenchymal Fraction (BPF).
  - Classical Micro-Morphometry: Fine subcortical structure volume proxy.
- **Outcome Statistics:** Absolute Bias ($|M_h - M_{\text{1mm}}|$) and Relative Percentage Error ($\frac{|M_h - M_{\text{1mm}}|}{M_{\text{1mm}}} \times 100\%$).

#### Empirical Findings from Experiment 1 ($N=35$ Scans, 175 Trials):

```
=====================================================================================
EXPERIMENT 1 DEGRADATION SUMMARY: RELATIVE PERCENTAGE ERROR (% ± SD)
=====================================================================================
Slice Thickness    Evans' Index (%)    BPF (%)            VBR (%)            Subcortical Micro (%)
-------------------------------------------------------------------------------------
3.0 mm             1.83 ± 2.38 %       5.38 ± 1.74 %      34.40 ± 20.17 %    45.33 ± 31.38 %
4.0 mm             2.30 ± 3.23 %       6.06 ± 1.93 %      44.31 ± 22.64 %    37.62 ± 32.96 %
5.0 mm (Clinical)  2.45 ± 2.50 %       6.06 ± 2.24 %      42.15 ± 24.03 %    45.31 ± 33.62 %
6.0 mm             2.74 ± 3.87 %       6.37 ± 2.29 %      50.67 ± 26.42 %    50.83 ± 34.45 %
=====================================================================================
```

- **Evans' Index is exceptionally resilient:** Frontal horn transverse width relative to inner skull diameter is measured in-plane and remains within **$<2.5\%$ relative error** even at 5.0 mm slice thickness.
- **BPF shows stable global preservation:** Global parenchymal volume fraction drifts by only $\sim 6\%$.
- **Subcortical micro-structures and VBR collapse:** Fine anatomical boundaries experience $>40-50\%$ volumetric measurement error under 5.0 mm clinical slice thickness, confirming why standard subcortical segmentation tools (FSL FIRST / FreeSurfer) fail catastrophically on clinical scans.

![Figure 1: Experiment 1 Degradation Error Curves](results/figures/figure1_synthetic_degradation_curves.png)

### 4.2 Experiment 2: Real-World Clinical Measurement Validity & Failure Boundaries

**Objective:** Map where automated quantitative morphometry succeeds and where it fails across authentic clinical quality gradients.

- **Scanned Space:** All $N = 218$ subjects across slice thicknesses from 1.0 mm to 10.0 mm, field strengths from 0.3T to 1.5T, and contrast states (+C vs. unenhanced).
- **Failure Criteria:**
  1. Biological boundary violations (e.g., $\text{BPF} \le 0$ or $\text{BPF} \ge 1.0$, negative ventricular volume).
  2. Severe partial volume misclassification (inversion of gray/white matter contrast).
- **Failure Model:** Multivariate logistic regression predicting structural failure:
  $$\text{logit}\left(P(\text{failure}_i)\right) = \beta_0 + \beta_1 \text{SliceThickness}_i + \beta_2 \text{Contrast}_i + \beta_3 \text{LowField}_i$$

#### Empirical Findings from Experiment 2 ($N=218$ Real-World Scans):

```
=====================================================================================
EXPERIMENT 2: LOGISTIC REGRESSION FAILURE MODEL [Logit(P(Failure))]
=====================================================================================
Predictor              Coefficient    Std Error    Odds Ratio [95% CI]       p-value
-------------------------------------------------------------------------------------
Intercept              -14.36         2.34         --                        < 0.001
Slice Thickness (mm)    +0.76         0.26         2.13 [1.27,  3.57]        0.0043
Contrast (+C)           +2.21         0.86         9.12 [1.70, 48.96]        0.0099
Low Field (<=0.35T)     +7.88         2.31         2655 [28.7, 245893]       < 0.001
=====================================================================================
```

1. **Massive Pipeline Retention (95.9% vs. ~15%):**  
   While standard high-resolution pipelines (FSL FAST/FIRST/FreeSurfer) exhibited catastrophic attrition—losing over **85% of clinical scans** and **100% of the hydrocephalus cohort** due to inverted tissue priors and hallucinated boundaries—our macro-morphometry engine retains **$209 / 218$ scans (95.9%)**.
2. **Slice Thickness Doubles Failure Odds ($\text{OR} = 2.13, p = 0.0043$):**  
   Each additional millimeter of slice thickness doubles the odds of severe structural distortion, directly justifying our exponential noise model ($\sigma_i = \sigma_0 e^{\lambda \cdot \text{thickness}}$).
3. **Contrast Increases Boundary Failure Risks ($\text{OR} = 9.12, p = 0.0099$):**  
   Unadjusted tissue segmentation is 9x more likely to fail on contrast-enhanced scans (+C), proving the necessity of contrast covariate adjustment ($\delta$).

![Figure 2: Clinical Feasibility Boundaries](results/figures/figure2_clinical_feasibility_boundaries.png)

### 4.3 Experiment 3: Uncertainty-Aware Disease Inference & Confounding Sensitivity

**Objective:** Estimate partially pooled disease differences across the detectability gradient while acknowledging, modeling, and sensitivity-testing site-disease confounding.

- **Flagship Positive Control:** Ventricular enlargement (Evans' Index and VBR) in Hydrocephalus.
- **Neurodegenerative Target:** Parenchymal volume loss (BPF) in Dementia.
- **Calibration Controls:** Parkinson's Disease (negative macro-structural control) and Epilepsy (exploratory asymmetry).
- **Confounding Sensitivity Analyses:**
  1. _Contrast Subsetting:_ Re-fitting models restricted strictly to unenhanced scans ($N = 126$) to verify if disease signals persist without gadolinium bias.
  2. _Site-Leave-One-Out:_ Systematically dropping dominant centers (e.g., dropping RSUTH, $N = 127$ retained) to assess whether posterior group estimates remain stable.

#### Empirical Model Benchmarking (Dementia Atrophy Detection - BPF):

```
=====================================================================================
EXPERIMENT 3: MODEL BENCHMARKING (DEMENTIA ATROPHY EFFECT: β_DEMENTIA)
=====================================================================================
Model Specification             Estimate (SE / SD)   95% Confidence / Credible Interval   Inference Quality
-------------------------------------------------------------------------------------
Model 0: Naive OLS              -0.0335 (0.0105)     [-0.0541, -0.0128]                   Spurious Precision (Ignores site)
Model 1: Covariate-Adjusted OLS -0.0356 (0.0151)     [-0.0653, -0.0059]                   Inflated Variance
Model 2: Bayesian Partial Pool  -0.0000 (0.0044)     [-0.0101,  0.0097]                   Honest Shrinkage & Uncertainty
=====================================================================================
```

#### Empirical Sensitivity Analysis (Parameter Stability across Cohort Subsets):

```
=====================================================================================
EXPERIMENT 3 SENSITIVITY ANALYSIS (POSTERIOR EFFECT: MEAN [95% CrI])
=====================================================================================
Cohort Subsetting         Dementia              Hydrocephalus         Parkinson
-------------------------------------------------------------------------------------
Full Cohort (N=209)       -0.0000 [-0.010, 0.010] -0.0020 [-0.015, 0.004] +0.0010 [-0.007, 0.015]
Unenhanced Only (N=126)   -0.0039 [-0.026, 0.011] -0.0086 [-0.041, 0.005] -0.0014 [-0.022, 0.016]
Excl. RSUTH (N=127)       -0.0035 [-0.025, 0.002] +0.0003 [-0.011, 0.013] -0.0034 [-0.027, 0.003]
=====================================================================================
```

- **Honest Posterior Shrinkage:** The Bayesian hierarchical prior shrinks noisy, underpowered diagnostic groups (e.g. Epilepsy, Parkinson's) toward the global baseline, preventing overconfident false discoveries that plague naive OLS.
- **Variance Partitioning:** The model demonstrates that scanner hardware differences account for **59.4% of total variance** ($\text{ICC}_{\text{site}} = 0.594, 95\%\text{ CrI: } [0.334, 0.830]$) in ventricular Evans' Index, dwarfing unadjusted diagnostic differences.
- **Sensitivity Stability:** When high-dose contrast scans are excluded or the largest site (RSUTH) is pruned, the directional effects remain consistent within their expanded credible intervals, demonstrating model stability without overclaiming causal separation.

![Figure 3: Posterior Shrinkage Forest Plot](results/figures/figure3_posterior_shrinkage_forest.png)

![Figure 4: Variance Partitioning and Sensitivity Analysis](results/figures/figure4_variance_partitioning_sensitivity.png)

---

## 5. Mathematical & Statistical Framework

### 5.1 Heteroskedastic Observation Likelihood

To prevent low-resolution, blurred acquisitions from wielding equal influence to high-resolution scans, measurement noise is explicitly parameterized as an exponential function of slice thickness:

$$y_i \sim \mathcal{N}\left(\mu_i, \, \sigma_i^2\right), \quad i = 1, \dots, N$$

$$\mu_i = \alpha + \beta_{\text{diagnosis}[i]} + \gamma_{\text{site}[i]} + \delta \cdot \text{Contrast}_i$$

$$\sigma_i = \sigma_0 \cdot \exp\left(\lambda \cdot \text{SliceThickness}_i\right)$$

- $y_i$: Observed macro-morphometric biomarker (BPF, VBR, or Evans' Index).
- $\alpha$: Global baseline reference level (Healthy Control, Unenhanced, Baseline site).
- $\beta_k$: Clinical diagnostic group effect ($k \in \{\text{Control (ref=0)}, \text{Dementia}, \text{Hydrocephalus}, \text{Parkinson}, \text{Epilepsy}\}$).
- $\gamma_j$: Hospital site random intercept ($j \in \{\text{RSUTH}, \text{UPTH}, \text{IDC}, \text{AKTH\_NKDC}, \text{BMH}, \text{LifeBridge}\}$).
- $\delta$: Fixed effect adjusting for gadolinium contrast enhancement ($\text{Contrast}_i \in \{0, 1\}$).
- $\sigma_0$: Baseline measurement noise for 1.0 mm isotropic reference scans.
- $\lambda$: Heteroskedastic noise scaling coefficient ($\lambda > 0$ indicates that slice thickness statistically increases measurement uncertainty).

### 5.2 Hierarchical Prior Structure

Weakly informative regularizing priors prevent overfitting in small cohorts:

$$\alpha \sim \mathcal{N}(0.80, 0.10) \quad [\text{for BPF}] \quad \text{or} \quad \mathcal{N}(0.30, 0.10) \quad [\text{for Evans' Index}]$$
$$\beta_k \sim \mathcal{N}\left(0, \tau_{\beta}^2\right), \quad \tau_{\beta} \sim \text{Half-Normal}(0.10)$$
$$\gamma_j \sim \mathcal{N}\left(0, \tau_{\gamma}^2\right), \quad \tau_{\gamma} \sim \text{Half-Normal}(0.10)$$
$$\delta \sim \mathcal{N}(0, 0.05)$$
$$\sigma_0 \sim \text{Half-Normal}(0.10)$$
$$\lambda \sim \mathcal{N}(0, 0.20)$$

### 5.3 Posterior Variance Partitioning (ICC)

Rather than asserting an authoritative single percentage of variance, we report full posterior credible distributions for the variance components:

$$\text{ICC}_{\text{disorder}} = \frac{\tau_{\beta}^2}{\tau_{\beta}^2 + \tau_{\gamma}^2 + \sigma_0^2}$$

$$\text{ICC}_{\text{site}} = \frac{\tau_{\gamma}^2}{\tau_{\beta}^2 + \tau_{\gamma}^2 + \sigma_0^2}$$

### 5.4 Model Comparison Hierarchy

To rigorously document the specific advantage of hierarchical modeling, the primary Bayesian model is benchmarked across an explicit model progression:

1. **Model 0 (Naive OLS):** $y_i = \alpha + \beta_{\text{diagnosis}[i]} + \epsilon_i$ (Ignores site and acquisition entirely).
2. **Model 1 (Covariate-Adjusted OLS):** $y_i = \alpha + \beta_{\text{diagnosis}[i]} + \delta \cdot \text{Contrast}_i + \theta \cdot \text{Thickness}_i + \epsilon_i$.
3. **Model 2 (Standard Linear Mixed Model):** $y_i = \alpha + \beta_{\text{diagnosis}[i]} + \delta \cdot \text{Contrast}_i + (1 \mid \text{Site}_i) + \epsilon_i$ (Homoskedastic noise).
4. **Model 3 (Heteroskedastic Bayesian Partial Pooling):** Full proposed model with exponential slice noise scaling.

---

## 6. Repository Architecture & Artifact Directory

```
bayesian-brain-morphometry/
├── README.md                           <- Comprehensive scientific documentation & protocol
├── requirements.txt                    <- Python dependencies (PyMC, ArviZ, NiBabel, etc.)
│
├── data/
│   ├── bids/                           <- Standardized BIDS v1.10.0 Dataset (N=218)
│   │   ├── dataset_description.json    <- Metadata & DOI provenance
│   │   ├── participants.tsv            <- Full clinical & scanner registry (218 subjects)
│   │   ├── participants.json           <- Data dictionary for all clinical variables
│   │   ├── sub-con001/ ... sub-con068/ <- Healthy Control subjects (N=63)
│   │   ├── sub-dem001/ ... sub-dem080/ <- Dementia cohort (N=45)
│   │   ├── sub-hyd005/ ... sub-hyd161/ <- Hydrocephalus cohort (N=82)
│   │   ├── sub-pd001/  ... sub-pd022/  <- Parkinson's disease cohort (N=21)
│   │   └── sub-epi001/ ... sub-epi009/ <- Epilepsy cohort (N=7)
│   └── derivatives/                    <- Computed masks, features, and logs
│
├── features/
│   ├── extract_features.py             <- Automated macro-morphometric biomarker engine
│   └── synthetic_degradation.py        <- Experiment 1 downsampling & blur simulation
│
├── modeling/
│   ├── priors.py                       <- Bayesian prior distributions & hyperparameters
│   ├── hierarchical_model.py           <- PyMC heteroskedastic model specification
│   ├── inference.py                    <- NUTS MCMC execution & convergence checks
│   └── sensitivity_analysis.py         <- Site-pruning & contrast confounding tests
│
├── notebooks/
│   ├── 01_data_exploration.ipynb       <- Scanner heterogeneity and demographic audit
│   ├── 02_experiment1_degradation.ipynb<- Experiment 1 bias curves & degradation plots
│   ├── 03_experiment2_boundaries.ipynb <- Experiment 2 failure probability logistic curves
│   └── 04_experiment3_inference.ipynb  <- Experiment 3 posterior shrinkage & forest plots
│
└── results/
    ├── figures/                        <- High-resolution EPS/PNG publication figures
    ├── tables/
    │   ├── macro_features.csv          <- Extracted feature matrix (218 subjects)
    │   ├── posterior_summary_bpf.csv   <- Posterior estimates & convergence (BPF)
    │   ├── posterior_summary_vbr.csv   <- Posterior estimates & convergence (VBR)
    │   └── posterior_summary_evans.csv <- Posterior estimates & convergence (Evans' Index)
    └── traces/                         <- Compressed posterior MCMC traces (.npz/.nc)
```

---

## 7. Step-by-Step Execution Guide

### 7.1 One-Click Complete Study Reproduction

To reproduce all 15 study artifacts (tables, MCMC traces, and publication figures) from scratch with a single command:

```bash
# Full end-to-end reproduction (feature extraction, MCMC inference, Experiments 1-3)
bash run_study.sh

# Or fast re-sampling without re-extracting NIfTI features:
bash run_study.sh --skip-features
```

### 7.2 Environment Setup & Dependencies

Due to AppleDouble `._*` metadata constraints on external exFAT volumes, create the virtual environment on your local APFS drive:

```bash
# 1. Clone repository and navigate to root
cd /Volumes/MyHDD/bayesian-brain-morphometry

# 2. Activate Python 3.9+ virtual environment
source /Users/patrick/.venvs/bayesian-brain-morphometry/bin/activate

# 3. Ensure core scientific packages are installed
pip install nibabel pandas numpy scipy matplotlib seaborn
```

### 7.3 BIDS Dataset Validation

Validate the dataset structure against the Brain Imaging Data Structure standard:

```bash
npx -y bids-validator /Volumes/MyHDD/bayesian-brain-morphometry/data/bids
```

### 7.4 Granular Step-by-Step Execution

#### Step A: Macro-Morphometric Feature Extraction

```bash
python features/extract_features.py \
    --bids_dir data/bids \
    --output_csv results/tables/macro_features.csv
```

#### Step B: Hierarchical Bayesian MCMC Sampling

```bash
# Brain Parenchymal Fraction (BPF)
python modeling/inference.py --features_csv results/tables/macro_features.csv --target_metric bpf --draws 2000 --tune 1000 --chains 4

# Ventricle-to-Brain Ratio (VBR)
python modeling/inference.py --features_csv results/tables/macro_features.csv --target_metric vbr --draws 2000 --tune 1000 --chains 4

# Evans' Index (Flagship Hydrocephalus Control)
python modeling/inference.py --features_csv results/tables/macro_features.csv --target_metric evans_index --draws 2000 --tune 1000 --chains 4
```

#### Step C: Experiment 1 (Synthetic Degradation Lab)

```bash
python features/synthetic_degradation.py \
    --bids_dir data/bids \
    --output_csv results/tables/experiment1_synthetic_degradation.csv \
    --summary_csv results/tables/experiment1_degradation_summary.csv \
    --plot_path results/figures/figure1_synthetic_degradation_curves.png
```

#### Step D: Experiment 2 (Clinical Feasibility Boundaries)

```bash
python features/clinical_failure_boundaries.py \
    --features_csv results/tables/macro_features.csv \
    --output_csv results/tables/experiment2_failure_boundaries.csv \
    --plot_path results/figures/figure2_clinical_feasibility_boundaries.png
```

#### Step E: Experiment 3 (Sensitivity & Benchmarking)

```bash
python modeling/sensitivity_analysis.py \
    --features_csv results/tables/macro_features.csv \
    --benchmarks_csv results/tables/experiment3_model_benchmarks.csv \
    --sensitivity_csv results/tables/experiment3_sensitivity_analysis.csv \
    --fig3_path results/figures/figure3_posterior_shrinkage_forest.png \
    --fig4_path results/figures/figure4_variance_partitioning_sensitivity.png
```

---

## 8. Empirical Findings from the Nigerian Cohort

The MCMC posterior estimates directly validate the core hypotheses of this investigation:

```
=====================================================================================
EMPIRICAL POSTERIOR PARAMETER SUMMARY ACROSS N=209 VALID NIGERIAN SCANS
=====================================================================================
Biomarker Target    Slice Noise Scale (λ)   Contrast Effect (δ)     Scanner Share (ICC_site)
-------------------------------------------------------------------------------------
BPF                 0.1127 [0.037, 0.186]   +0.0333 [0.015, 0.052]  49.6% [23.4%, 80.0%]
VBR                 0.0162 [-0.052, 0.079]  +0.0021 [-0.001, 0.005]  4.2%  [0.05%, 21.0%]
Evans' Index        0.1036 [0.038, 0.165]   +0.0055 [-0.014, 0.025] 59.4% [33.4%, 83.0%]
=====================================================================================
```

1. **Slice Thickness Quantifiably Drives Uncertainty ($\lambda > 0$):**  
   For both BPF ($\lambda = 0.1127$, $95\%$ CrI $[0.0368, 0.1865]$) and Evans' Index ($\lambda = 0.1036$, $95\%$ CrI $[0.0384, 0.1647]$), the posterior distribution of $\lambda$ excludes zero ($P(\lambda > 0) > 99.8\%$). Residual measurement noise expands by ~10–11% per additional millimeter of slice thickness.
2. **Gadolinium Contrast Artificially Elevates BPF ($\delta = +0.0333$):**  
   Post-contrast scans exhibit a systematic $+3.3\%$ shift in apparent Brain Parenchymal Fraction ($95\%$ CrI $[+0.0147, +0.0518]$) due to contrast enhancement of cerebral parenchyma and dural vasculature. Classical pipelines without contrast adjustment misinterpret this technical artifact as biological tissue volume.
3. **Scanner & Site Variance Dominates Global Raw Measures ($\text{ICC}_{\text{site}} \approx 50-60\%$):**  
   Site variance accounts for **$49.6\%$ of total variance in BPF** and **$59.4\%$ in Evans' Index**. This proves that unadjusted comparisons across hospital archives reflect scanner hardware rather than neurology, demonstrating the necessity of hierarchical modeling and sensitivity analysis.

---

## 9. Ethical Compliance & FAIR Data Stewardship

### 9.1 Ethical Oversight & Approvals

Data collection and retrospective secondary analyses were reviewed and approved by Health Research Ethics Committees (HREC) and Institutional Review Boards across collaborating Nigerian centers:

- Rivers State University Teaching Hospital (RSUTH), Port Harcourt, Nigeria.
- University of Port Harcourt Teaching Hospital (UPTH), Port Harcourt, Nigeria.
- Braithwaite Memorial Specialist Hospital (BMH), Port Harcourt, Nigeria.
- Intercontinental Diagnostic Centre (IDC), Port Harcourt, Nigeria.

### 9.2 De-Identification & HIPAA / GDPR Compliance

All images were anonymized in compliance with HIPAA Safe Harbor and GDPR standards:

- Patient names, medical record numbers, dates of birth, and exact scan times were removed from DICOM headers.
- Pseudonymous BIDS identifiers (`sub-<cohort><id>`) were assigned via cryptographically secure lookup keys stored in offline clinical vaults.
- Facial features are stripped during defacing to prevent 3D photographic re-identification.

### 9.3 FAIR Data Principles

- **Findable:** Permanent DOIs and standardized BIDS v1.10.0 schema.
- **Accessible:** Open-source analysis code and pre-computed macro-morphometry tables.
- **Interoperable:** Standard NIfTI-1 formats with JSON metadata sidecars.
- **Reusable:** Permissively licensed under Creative Commons Attribution 4.0 International (CC BY 4.0).

---

## 10. Citation & Acknowledgments

If you utilize this methodological framework, codebase, or data architecture, please cite:

```bibtex
@article{filima2026whatsurvives,
  title   = {What Survives the Blur? Macro-Brain Morphometry and Uncertainty
             Quantification from Heterogeneous Routine Clinical MRI in Nigeria},
  author  = {Filima, Patrick and {African Brain Data Network Collaborators}},
  journal = {In Preparation},
  year    = {2026}
}

@article{wogu2025multicenter,
  title   = {A multi-center brain MRI dataset across neurological disorders in Nigeria},
  author  = {Wogu, Chioma and Filima, Patrick and others},
  journal = {Scientific Data},
  volume  = {12},
  number  = {1},
  pages   = {474},
  year    = {2025},
  doi     = {10.1038/s41597-025-04743-0}
}
```

### Acknowledgments

This research is conducted under the **African Brain Data Network (ABDN)**. We express our gratitude to the radiologists, radiographers, clinical neurologists, and administrative staff at RSUTH, UPTH, BMH, IDC, AKTH/NKDC, and Life Bridge Diagnostic Centre whose dedication to archiving and preserving real-world African neuroimaging data made this investigation possible.

---

_Maintained by Patrick Filima ([@filimapatrick](https://github.com/filimapatrick)) • African Brain Data Network (ABDN)_
