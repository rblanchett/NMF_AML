# Analysis decisions and scientific rationale

This page tracks **what was done**, **why**, **what to interpret cautiously**, and **what remains to verify**. It documents existing analytic choices; it should not be mistaken for a post-hoc claim that all settings were preregistered.

## Cohort and modeling context

- 854 pediatric AML patients (400 female, 454 male in the previously analyzed cohort).
- 729,146 assayed CpG features in the combined methylation beta matrix used for NMF after quality control.
- Missing beta values were masked in the original NMF fit rather than filled for NMF; imputation used for distinct PCA exploratory work is not part of this UMAP analysis.
- A 23-component NMF H matrix was retained for Script 16. Two components may have zero patients assigned by maximum H score, but their weights remain part of other patients' profiles.
- The NMF rank-selection analysis used cross-validation, separately from the UMAP workflow. **Verify primary and sensitivity holdout analyses and their exact selection results before describing a definitive rank-selection conclusion in a manuscript.**

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
| Metadata | Join by USI internally | Deterministic alignment | Original IDs must not be exported unintentionally |
| NMF confidence | Original Margin and patient-relative Margin | Explore strength of dominant component assignment | Margin reflects model weight dominance, not a calibrated probability |
| Binary lesion handling | Present / Absent / Unknown | Keep missing information separate from confirmed negative | Verify coding for each source column before recoding |
| Sex/risk palette | Blue/pink for male/female; green/yellow/red for low/standard/high | Consistent visuals across covariate-specific figures | Exact risk-label compatibility requires checking the input field |
| Inference | UMAP visualization only | Exploratory view of patient similarity | Group and clinical associations require their original statistical analyses |

## Important distinction about risk groups

`Final.Risk.Group` has the established Low / Standard / High grouping in the prior clinical analyses. Script 16's public example configuration instead uses `Risk.group.classification..by.protocol.1031`. Verify that column's actual values and coding before claiming it is equivalent to `Final.Risk.Group` or before using the configured palette in final figures.

## Planned supporting figure groups (not confirmed results)

- **S01**: Raw vs normalized H colored by NMF group.
- **S02**: Same embeddings colored by sex.
- **S03**: Risk annotation with source-specific category definitions.
- **S04**: NMF Margin and relative Margin, using continuous scales.
- **S05**: Selected molecular lesion overlays with explicit Unknown category.
- **S06**: Sensitivity comparison across UMAP parameters, with neighbor-overlap QC.

Each caption should explicitly state whether a plot uses raw or normalized H, that the coordinates are derived from the same fitted model, and that labels are annotations rather than supervised predictors of the coordinates.

## Verification checklist before publication

- [ ] Check final cohort totals, sample inclusion criteria, and underlying dataset accessions.
- [ ] Verify rank-selection strategy and selected `k` against saved rank-selection results.
- [ ] Verify H dimensions and patient order from the actual saved analysis objects.
- [ ] Check risk field levels and binary-coding conventions.
- [ ] Confirm versions of R, uwot, Plotly, ggplot2, and relevant dependencies from actual session info.
- [ ] Record the exact QC scores from `All_UMAP_Parameter_QC.csv` (do not invent values).
- [ ] Verify every public image and HTML export against source-specific redistribution terms.
- [ ] Document figure output filenames, plot generation date, and script commit ID.
