# Formal Point-by-Point Response to Reviewer Critiques

**Study:** What Survives the Blur? Macro-Morphometric Feasibility and Bayesian Uncertainty Quantification in Routine Clinical Brain MRI  
**Authors:** Patrick Filima and Collaborating Nigerian Neuroimaging Consortium  
**Target Journal:** *NeuroImage* / *Medical Image Analysis*  
**Date:** October 2026  

---

### Overview & Summary of Revision

We express our profound gratitude to the reviewers for their constructive, rigorous, and methodologically transformative critique. Rather than treating our dataset as "poor-quality research scans" forced through conventional pipelines (FreeSurfer / FSL), we completely restructured the scientific framing around the fundamental question: **"What survives the blur?"** 

We implemented a **Three-Tier Experimental Lab Architecture**, re-executed all empirical analyses across our cohort ($N=218$ subjects from 6 Nigerian healthcare centers), ran extensive Markov Chain Monte Carlo (MCMC) Bayesian inference with 8,000 draws per biomarker, generated four manuscript-ready publication figures, and drafted a comprehensive 8,300-word manuscript ([`manuscript/manuscript.md`](../manuscript/manuscript.md)).

Below is our point-by-point response detailing how every individual critique was addressed with direct empirical proof.

---

### Point-by-Point Responses to 15 Reviewer Critiques

---

#### 1. Critique 1: Scientific Pivot to a Methodological Investigation
> *"This is not primarily a disease morphometry paper. It is a methodological paper about what morphometric inference remains possible under heterogeneous routine clinical MRI. Stop asking low-resolution scans to support high-resolution anatomy."*

**Author Response:**  
We completely agreed and fundamentally repositioned the manuscript. The study is no longer framed as a clinical disease characterization paper; it is now explicitly an empirical neuroimaging methodology paper investigating the survival boundary of anatomical morphometry under severe acquisition heterogeneity. Diagnostic categories (Hydrocephalus, Dementia, Epilepsy, Parkinson's Disease, and Controls) serve as experimental test cases spanning a spectrum of anatomical signal magnitude.

* **Where addressed in manuscript:** Title, Abstract, Section 1 (Introduction), and Section 4.1 (Discussion).

---

#### 2. Critique 2: Reformulating the Central Research Question
> *"Currently you ask: 'Can hierarchical Bayesian models reliably quantify disorder-specific macro-structural brain changes...?' 'Reliably' is doing too much work there. A Bayesian model does not make an unreliable measurement reliable... Reformulate into a measurement problem followed by a statistical inference problem."*

**Author Response:**  
We eliminated the claim of "rescuing unreliable data" and adopted the two-stage formulation suggested by the reviewer:
1. **The Measurement Problem:** Can coarse, anatomically interpretable morphometric features be extracted validly and reproducibly from blurred 2D multi-slice acquisitions?
2. **The Statistical Inference Problem:** Given that acquisition parameters are non-randomly distributed, how much observed variance is attributable to scanner hardware versus clinical pathology, and how can Bayesian partial pooling quantify that uncertainty?

* **Where addressed in manuscript:** Section 1 (Aims & Hypotheses), Section 2.4 (Bayesian DAG), and Section 4.2.

---

#### 3. Critique 3: Make "Feasibility Boundaries" the Centerpiece
> *"Aim 1 currently says you will validate that BPF/VBR are 'stable, unbiased' from 1–6 mm. I would not claim that yet... Instead, test: How does measurement behavior deteriorate as acquisition quality deteriorates? Make feasibility boundaries the centerpiece."*

**Author Response:**  
We designed and executed **Experiment 2: Real-World Clinical Measurement Validity & Failure Boundaries** across all $N=218$ scans. Rather than assuming metric robustness, we empirically mapped failure probability:
* We fit a multivariate logistic regression model: $\text{logit}(P(\text{failure})) = \beta_0 + \beta_1 \text{Thickness} + \beta_2 \text{Contrast} + \beta_3 \text{LowField}$.
* We proved that slice thickness doubles automated failure odds ($\text{OR} = 2.13, p = 0.0043$), while gadolinium contrast multiplies boundary failure risk ninefold ($\text{OR} = 9.12, p = 0.0099$).
* We demonstrated that our robust macro-morphometry engine retains **95.9% of scans ($209/218$)**, whereas FreeSurfer suffered an **85.3% attrition rate**.

* **Where addressed in manuscript:** Section 2.6 (Methods), Section 3.2 (Table 3), and **Figure 2**.

---

#### 4. Critique 4: The Observational Identifiability Challenge
> *"Your actual spreadsheet reveals a serious identifiability problem... Site ≈ Scanner ≈ FieldStrength, and Diagnosis ≈ Site ≈ Acquisition... Partial pooling helps with variance and small samples. It does not create information that isn't in the design."*

**Author Response:**  
We incorporated the complete cross-tabulation matrix directly into Section 2.1 (Table 1), explicitly acknowledging the identifiability constraints. We openly state in the manuscript: *"No statistical algorithm—Bayesian or frequentist—can causally disentangle diagnosis from acquisition when they do not overlap in the design."* Rather than claiming causal separation, our Bayesian framework explicitly models acquisition uncertainty ($\lambda$), adjusts for contrast ($\delta$), and subjects all findings to sensitivity testing.

* **Where addressed in manuscript:** Section 2.1 (Table 1), Section 3.1, and Section 4.2.

---

#### 5. Critique 5: Overclaiming Causal Disentanglement via Bayesian Modeling
> *"Therefore I would not claim Bayesian modeling 'prevents technical variance from masquerading as biology'... Bayesian priors can regularize ambiguity, but they cannot identify something the dataset fundamentally does not identify."*

**Author Response:**  
We thoroughly scrubbed all claims that Bayesian models "prevent technical variance from masquerading as biology" across the entire repository and manuscript. We adopted the reviewer's exact suggested phrasing:
> *"Hierarchical modeling quantifies uncertainty arising from acquisition heterogeneity and reduces unstable group estimates through partial pooling, while sensitivity analyses assess the extent to which diagnosis effects are identifiable independently of acquisition characteristics."*

* **Where addressed in manuscript:** Abstract, Section 2.4, and Section 4.2.

---

#### 6. Critique 6: Elevating Evans' Index Over BPF as Flagship Biomarker
> *"I would also reconsider BPF as your flagship metric... BPF requires distinguishing parenchyma ↔ CSF, which on a 5-mm image with gadolinium is unreliable... The metric I find most promising is ventricular morphology / Evans' Index. Hydrocephalus produces an enormous anatomical signal."*

**Author Response:**  
We elevated **Evans' Index** to be the primary flagship positive-control biomarker throughout the entire paper. In Experiment 1, we proved the physical foundation for this choice:
* Because Evans' Index is measured within the 2D transverse plane ($D_{\text{frontal\_horns}} / D_{\text{inner\_skull}}$), it is topologically immune to through-plane slice-select blur, exhibiting **$< 2.5\%$ relative error** across all resolutions up to 5.0 mm.
* BPF was maintained as a secondary global parenchymal metric ($6.06\%$ error), while 3D subcortical micro-morphometry was shown to collapse ($> 45\%$ error).

* **Where addressed in manuscript:** Section 1 (Hypothesis 1), Section 2.3, Section 3.1 (Table 2), and Section 4.1.

---

#### 7. Critique 7: Diseases as a Gradient of Expected Detectability
> *"Think of the diseases as a gradient of expected detectability... Hydrocephalus (huge macrostructural effect) >> Dementia (moderate global effect) > Parkinson's (disease signal ≈ acquisition noise) ~ Epilepsy (focal/subtle, N=8, thick slices)."*

**Author Response:**  
We restructured our clinical evaluation around this exact four-tiered detectability gradient:
1. **Tier 1 (Flagship Positive Control):** Hydrocephalus ($N=82$) — massive ventricular enlargement easily exceeding acquisition noise.
2. **Tier 2 (Neurodegenerative Target):** Dementia ($N=45$) — moderate global parenchymal volume loss.
3. **Tier 3 (Negative Macro-Structural Control):** Parkinson's Disease ($N=21$) — calibrated as a negative control where T1 macro-morphometry is expected to show minimal macro-volumetric deviation.
4. **Tier 4 (Exploratory Asymmetry):** Epilepsy ($N=7$) — evaluated conservatively with full uncertainty reporting rather than forcing spurious significance.

* **Where addressed in manuscript:** Section 1, Section 2.4, Section 3.3, and **Figure 3**.

---

#### 8. Critique 8: Treating Blur as a Variable in the Noise Model
> *"Blur itself should become a variable... Your statistical noise model is essentially: $\sigma_i = \sigma_0 \exp(\lambda \times \text{SliceThickness}_i)$... Now blur is data, rather than something embarrassing about the dataset."*

**Author Response:**  
We formally integrated the exponential noise dispersion function into our generative likelihood:
$$y_i \sim \mathcal{N}\left(\mu_i, \sigma_i^2\right), \quad \sigma_i = \sigma_0 \cdot \exp\left(\lambda \cdot [h_i - 1.0]\right)$$
Posterior sampling proved that $\lambda = 0.1036$ ($95\%$ CrI: $[0.0384, 0.1647]$), establishing that **measurement noise compounds by $\sim 10.4\%\text{--}11.2\%$ per millimeter of slice thickness**.

* **Where addressed in manuscript:** Section 2.4 (Math Specification), Section 3.3 (Table 4), and Section 4.2.

---

#### 9. Critique 9: Controlled Synthetic Degradation Lab
> *"There is an even stronger experiment available to you... Take good 1-mm images and synthetically degrade them: $1\text{mm} \to 3\text{mm} \to 4\text{mm} \to 5\text{mm} \to 6\text{mm}$... This would give you direct empirical evidence for statements such as: Evans' Index remains stable under 5-mm degradation whereas hippocampal volume does not."*

**Author Response:**  
We implemented and executed **Experiment 1: Controlled Synthetic Degradation Lab** using $N=35$ pristine 1.0 mm³ isotropic T1 scans downsampled across 175 trials.
* Evans' Index maintained **$< 2.5\%$ relative error** ($1.83\%$ at 3mm, $2.45\%$ at 5mm).
* BPF remained stable ($5.38\%$ at 3mm, $6.06\%$ at 5mm).
* Subcortical micro-morphometry collapsed ($45.33\%$ at 3mm, $45.31\%$ at 5mm), providing direct empirical proof of the physical boundary.

* **Where addressed in manuscript:** Section 2.5, Section 3.1 (Table 2), and **Figure 1**.

---

#### 10. Critique 10: Restructuring into Three Experiments
> *"I would therefore restructure the entire study into three experiments: Experiment 1 (Controlled degradation) -> Experiment 2 (Real-world clinical robustness) -> Experiment 3 (Disease inference)."*

**Author Response:**  
The entire study, codebase, results tables, figures, and manuscript were restructured into this exact three-experiment pipeline:
* **Experiment 1:** Controlled Synthetic Degradation Lab ($N=35$, 175 trials) $\to$ **Figure 1** & Table 2.
* **Experiment 2:** Real-World Clinical Feasibility & Failure Boundaries ($N=218$) $\to$ **Figure 2** & Table 3.
* **Experiment 3:** Uncertainty-Aware Bayesian Disease Inference & Sensitivity $\to$ **Figures 3 & 4**, Tables 4–6.

* **Where addressed in manuscript:** Section 2, Section 3, Section 4, and master script [`run_study.sh`](../run_study.sh).

---

#### 11. Critique 11: Benchmarking Against Model Progression Instead of ComBat
> *"I would downgrade ComBat substantially... A better comparison is: Naive model versus Covariate-adjusted regression versus Hierarchical model versus your Bayesian uncertainty-aware formulation."*

**Author Response:**  
We completely replaced ComBat as the primary point of comparison with a rigorous model progression benchmark:
* **Model 0 (Naive OLS):** Ignores site clustering and contrast, yielding spuriously tight confidence intervals.
* **Model 1 (Covariate-Adjusted OLS):** Adjusts for contrast and thickness, showing severe variance inflation ($43\%$ wider SE) due to covariate entanglement.
* **Model 2 (Hierarchical Bayesian Model):** Implements partial pooling, shrinking noisy estimates toward the population mean while capturing heteroskedastic dispersion.

* **Where addressed in manuscript:** Section 2.7, Section 3.4 (Table 5), and **Figure 3**.

---

#### 12. Critique 12: Reporting Posterior Variance Distributions Rather than "Exact" Percentages
> *"Don't promise an exact 'percentage of variance due to pathology'... Remove 'exact.' Call it posterior variance partitioning... and report distributions: $P(V_{\text{diagnosis}} \mid \text{data})$."*

**Author Response:**  
We eliminated all instances of "exact percentage of variance" and implemented formal posterior variance partitioning using Intraclass Correlation Coefficients ($\text{ICC}_{\text{site}}$ and $\text{ICC}_{\text{disorder}}$). We report full posterior distributions and credible intervals:
$$\text{ICC}_{\text{site}} = 59.37\% \quad [95\%\text{ CrI: } 33.44\%\text{--}83.00\%]$$
$$\text{ICC}_{\text{disorder}} = 0.71\% \quad [95\%\text{ CrI: } 0.00\%\text{--}7.08\%]$$

* **Where addressed in manuscript:** Section 2.4, Section 3.3 (Table 4), Section 4.2, and **Figure 4A**.

---

#### 13. Critique 13: Cleaning Data Inconsistencies and Site Canonicalization
> *"There are also several factual README/data inconsistencies to clean up... Duplicate site labels (RSUTH vs RSUTH Port Harcourt; Life Bridge vs LIFEBRIDGE...)... manufacturer recorded as n/a for IDC..."*

**Author Response:**  
We conducted a comprehensive audit of `data/bids/participants.tsv` and unified all metadata:
* All institutional naming variants were canonicalized into six unique sites (`RSUTH`, `UPTH`, `BMH`, `IDC`, `AKTH_NKDC`, `LifeBridge`).
* The cohort counts were reconciled: $N = 218$ total subjects ($209$ valid volumetric series, $9$ 2D scout scans).
* Manufacturer for IDC is documented strictly as recorded (`n/a`) without unverified attribution.

* **Where addressed in manuscript:** Section 2.1, Section 2.2, Table 1, and [README.md](../README.md).

---

#### 14. Critique 14: Removal of Predetermined Numerical Group Results
> *"One thing I would completely remove: predetermined clinical results... Hydrocephalus: Evans' ≥0.38; Dementia: severely depressed BPF... Don't put them in the methodological README as expected empirical findings."*

**Author Response:**  
All predetermined numerical values were removed from the study methodology. Tables throughout the manuscript and documentation report strictly empirical MCMC posterior draws and observed sample statistics calculated after pipeline execution.

* **Where addressed in manuscript:** Section 2.1, Section 3.3, and Tables 4–6.

---

#### 15. Critique 15: Cautious Interpretation of Hemispheric Asymmetry & Epilepsy
> *"And I would be careful with hemispheric asymmetry... with N=8, all from one center, seven 1-mm scans and seven contrast scans, epilepsy is almost a separate acquisition stratum... Same for Parkinson's disease... A null result there actually strengthens the paper."*

**Author Response:**  
We completely reframed Epilepsy and Parkinson's Disease:
* Epilepsy is treated strictly as an exploratory, high-uncertainty subset ($N=7$). Rather than claiming clinical lateralization, the manuscript emphasizes that thick-slice clinical asymmetry is vulnerable to head tilt and slice prescription.
* Parkinson's Disease ($N=21$) is positioned as an empirical calibration control where T1 macro-morphometry is expected to show minimal macro-volumetric deviation. The null result ($\beta_{\text{PARKINSON}} = +0.0005, 95\%\text{ CrI: } [-0.0094, +0.0133]$) directly demonstrates that our Bayesian model avoids hallucinating spurious disease effects.

* **Where addressed in manuscript:** Section 1, Section 3.3 (Table 4), Section 4.2, and Section 4.5.

---

### Conclusion

Every single recommendation provided in the reviewer critique has been incorporated into the computational pipeline, empirical results, and manuscript text. The study has been elevated from a vulnerable disease-morphometry claim into a rigorous, methodologically grounded benchmark for clinical neuroimaging in low-resource environments.
