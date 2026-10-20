# Analysis decisions and scientific rationale

This page tracks **what was done**, **why**, **what to interpret cautiously**, and **what remains to verify**. It documents existing analytic choices; it should not be mistaken for a post-hoc claim that all settings were preregistered.

## Cohort and modeling context

- 854 pediatric AML patients (400 female, 454 male in the previously analyzed cohort).
- 729,146 assayed CpG features in the combined methylation beta matrix used for NMF after quality control.
- Missing beta values were masked in the original NMF fit rather than filled for NMF
- A 23-component NMF H matrix was retained for Script 16. Two components may have zero patients assigned by maximum H score, but their weights remain part of other patients' profiles.
- The NMF rank-selection analysis used cross-validation, separately from the UMAP workflow.

## Decision register

| Topic | Selected approach | Reason | Interpretation / caveat |
|---|---|---|---|
| Input to UMAP | 854 x 23 transpose of fitted H | Retain full patient component profiles rather than just maximum-assignment labels | Represents NMF space, not independent CpG-space evidence |
| Component preservation | Keep 23 features | Zero-dominant-assignment components may still have nonzero weights | Cannot discard on the basis of assignment count alone |
| H normalization | Analyze raw H and each patient's weights divided by their sum | Compare magnitude-sensitive and composition-like views | NMF component scaling is not uniquely identifiable; normalized profiles do not fully resolve component-specific scale ambiguity |
| Distance metrics | Euclidean for raw H, cosine for normalized H | Contrast weight magnitude and profile direction | Simultaneously changing normalization and metric complicates causal attribution of observed differences |
| UMAP primary settings | 30 neighbors, min_dist 0.30, spread 1, 500 epochs, seed 190 | Consistent reproducible reference visualization | Settings are exploratory rather than proof of biological groups |
| Sensitivity grid | Neighbors 15/30/50; min_dist 0.1/0.3/0.5 | Check robustness to display parameters | Do not select settings to maximize visually appealing clinical separation |
| Neighborhood QC | Mean overlap of 15-nearest neighbors between input H and UMAP | Quantify local neighbor preservation | Does not establish biological validity or global distance fidelity |
| NMF confidence | Original Margin and patient-relative Margin | Explore strength of dominant component assignment | Margin reflects model weight dominance, not a calibrated probability |
| Binary lesion handling | Present / Absent / Unknown | Keep missing information separate from confirmed negative | Verify coding for each source column before recoding |
| Sex/risk palette | Blue/pink for male/female; green/red for low/high | Consistent visuals across covariate-specific figures | Exact risk-label compatibility requires checking the input field |


