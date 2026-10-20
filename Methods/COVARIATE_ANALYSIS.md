
## NMF Group-Specific Clinical and Molecular Associations

### Overview

To identify clinical and molecular characteristics associated with individual NMF groups, each group was compared against **all remaining NMF groups combined** (one-vs-all analysis).

Statistical significance was assessed using the **Benjamini–Hochberg (BH) false discovery rate (FDR) correction**, with an adjusted significance threshold of **q < 0.05**.

### Key Findings

- **13 NMF groups** demonstrated statistically significant enrichment or differences for at least one clinical or molecular covariate.
- For each group, results are presented as:
  1. A list of statistically significant covariates.
  2. A contextual description of the group's clinical and molecular characteristics relative to the overall dataset.
- These analyses were intended to identify distinguishing characteristics of individual methylation-derived NMF groups.

> **Important observation — NMF18:**  
> NMF18 consisted exclusively of female patients. However, sex was the only covariate that reached statistical significance after multiple-testing correction. This association was identified by comparing NMF18 against the remaining cohort, which included both female and male patients. No other evaluated covariates reached the significance threshold.

---

### Statistical Methodology

Statistical tests were selected according to covariate type and the distribution of observed categorical counts.

| Covariate Type | Statistical Test | Selection Criteria |
|----------------|------------------|--------------------|
| **Continuous / Numeric** | Welch's two-sample t-test | Used to compare group means without assuming equal variances. |
| **Binary Categorical (2 × 2)** | Pearson's χ² test | Applied when all contingency table cell counts were ≥ 5. |
| **Binary Categorical (2 × 2)** | Fisher's exact test | Applied when at least one cell count was < 5. |
| **Multilevel Categorical (> 2 levels)** | Pearson's χ² test | Applied when all contingency table cell counts were ≥ 5. |
| **Multilevel Categorical (> 2 levels)** | Monte Carlo simulated χ² test | Applied when at least one cell count was < 5, using **2,000 simulations**. |

### Multiple-Testing Correction

- P-values were adjusted using the **Benjamini–Hochberg procedure**.
- Associations were considered statistically significant at **FDR q < 0.05**.
- Significant findings were summarized independently for each NMF group, with accompanying cohort-level context.

### Interpretation

These one-vs-all analyses characterize the clinical and molecular features that distinguish each NMF group from the remaining study population.

Importantly, statistically significant differences indicate associations with NMF group membership, but do not necessarily establish that a covariate is exclusive to a particular group or causally responsible for the underlying DNA methylation patterns.
.

