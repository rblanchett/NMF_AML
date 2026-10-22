# Pediatric AML DNA methylation: NMF-derived subgroup analysis and interactive UMAP exploration

![Status](https://img.shields.io/badge/analysis-ongoing-informational) ![Language](https://img.shields.io/badge/language-R-blue)

**Code and methodological documentation for exploratory methylation-based stratification of pediatric acute myeloid leukemia (AML).** The project uses nonnegative matrix factorization (NMF) to represent patient methylation profiles, evaluates associations with clinical and molecular annotations, and visualizes patient relationships in a two-dimensional UMAP embedding.

> **Scope:** This repository is being populated in stages. At present, the public release focuses on **Script 16 (H-matrix UMAP)** and the **Section 13/14 explorer exporters**. Earlier NMF, survival, and lesion-association scripts will be added at a later point.  

## Quick navigation

| Resource | Description | Link |
|---|---|---|
| Analysis code | Public-facing R scripts, with configurable file paths | [Scripts/](Scripts/) |
| Analysis decisions | Rationale, alternatives, and methodological cautions | [Methods/ANALYSIS_DECISIONS.md](Methods/ANALYSIS_DECISIONS.md) |
| UMAP methods | Detailed reproducible description of Script 16 | [Methods/UMAP_METHODS.md](Methods/UMAP_METHODS.md) |
| Covariate Analysis Methods | Figure naming conventions, legends, and evidence trail | [Methods/COVARIATE_METHODS.md](Methods/COVARIATE_METHODS.md) |
| Supplementary inventory | Track figure/table sources and publication status | [Supplement/FIGURE_MANIFEST.md](Supplement/FIGURE_MANIFEST.md) |
| Interactive explorer | Public GitHub Pages page | [Interactive UMAP](docs/index.html) |

## Study overview

The analysis was developed using a **combined-sex pediatric AML cohort (854 patients)** and a DNA methylation beta matrix of approximately **729,146 CpGs** after quality control. The existing fitted NMF decomposition uses **23 components**. The H matrix contains a 23-component profile for each patient and supports exploration beyond the discrete group assigned from the largest H score. These are descriptions of the analyzed cohort and stored model, **not a new independent validation dataset**.

## Analysis workflow

```text
Input CpG x patient methylation matrix (prior analysis)
        |
        v
NMF factorization: V ~ W H (existing 23-component model)
        |
        v
H: 23 components x 854 patients
        |
        +--> Patient-by-component raw H --> Euclidean UMAP
        |
        +--> Patient-sum-normalized H --> Cosine UMAP
                    |
                    v
      2D patient coordinates (854 rows per embedding)
                    |
                    v
        Merge clinical/molecular annotations by USI
                    |
          +---------+--------------------+
          |                              |
          v                              v
Publication-quality static figures   Private and reviewed public explorer
(PDF / PNG; local outputs)           (interactive HTML)
```

## Software and reproducibility

The public scripts use R with `uwot`, `ggplot2`, `plotly`, and other packages listed at the beginning of each script. The original run saved an `R_Session_Info.txt` file locally.

## Current scripts

| Script | Purpose | Dependency |
|---|---|---|
| [`16_NMF_H_UMAP_Interactive_Exploration_PUBLIC.R`](scripts/16_NMF_H_UMAP_Interactive_Exploration_PUBLIC.R) | Prepares raw/normalized H, fits UMAP, exports QC and figures | Existing NMF H and clinical input files |
| [`16_Section_13_Unified_Explorer_PUBLIC_APPEND_ONLY.R`](scripts/16_Section_13_Unified_Explorer_PUBLIC_APPEND_ONLY.R) | Append-only unified explorer code; export disabled by default | Original Script 16 runtime objects |
| [`16_Section_14_Public_GitHub_UMAP_Exporter.R`](scripts/16_Section_14_Public_GitHub_UMAP_Exporter.R) | Standalone public-export workflow using saved Script 16 embeddings; removes original USI and survival variables from its selected HTML fields | Existing `Script_16_Complete_Analysis_Objects.rds`, plus source-specific release approval |

These two explorer approaches serve different integration needs; **the standalone Section 14 does not require appending Section 13 or rerunning NMF/UMAP**. Do not execute the two exporters merely to create redundant published artifacts.

## Key methodological decisions

| Decision | Implementation / rationale |
|---|---|
| Keep all components | All 23 H rows retained, including components that are not the dominant assignment for any patient |
| Use H, not the full methylation matrix | Inspect the 23-component NMF representation while retaining each patient's full weight profile |
| Use two representations | Raw H and per-patient sum-normalized H, analyzed separately |
| Distances | Euclidean on raw H, cosine on normalized H in the current Script 16 configuration |
| Reproducibility | Seed 190; one UMAP thread in the current configuration |
| Primary UMAP | `n_neighbors = 30`, `min_dist = 0.30`, `spread = 1`, `n_epochs = 500` |
| Sensitivity | Neighborhood sizes 15/30/50, min_dist 0.1/0.3/0.5; neighborhood-overlap QC (k=15) |
| Missing lesion calls | Unknown remains separate from lesion-absent |
| Visualization covariates | Sex, NMF group, risk, molecular/cytogenetic characteristics, and assignment Margin |
| Survival | Not used to calculate UMAP; excluded from the public Section 14 explorer |

See [Analysis Decisions](methods/ANALYSIS_DECISIONS.md) for interpretation and limitations. This table documents configured methods; it is **not a table of scientific results**.

## Supporting figures

Supporting figures are placed in [`Figures/`](Figures/).  

## Data accessions and study citations:
Pollard JA, Alonzo TA, Gerbing R, Brown P, Fox E, Choi J, Fisher B, Hirsch B, Kahwash S, Getz K, Levine J, Brodersen LE, Loken MR, Raimondi S, Tarlock K, Wood A, Sung L, Kolb EA, Gamis A, Meshinchi S, Aplenc R. Sorafenib in Combination With Standard Chemotherapy for Children With High Allelic Ratio FLT3/ITD+ Acute Myeloid Leukemia: A Report From the Children's Oncology Group Protocol AAML1031. J Clin Oncol. 2022 Jun 20;40(18):2023-2035. doi: 10.1200/JCO.21.01612. Epub 2022 Mar 29. PMID: 35349331; PMCID: PMC9197362.
	
Bertrums EJM, Smith JL, Harmon L, Ries RE et al. Comprehensive molecular and clinical characterization of NUP98 fusions in pediatric acute myeloid leukemia. Haematologica 2023 Aug 1;108(8):2044-2058. PMID: 36815378

Marchi F, Shastri VM, Marrero RJ, Nguyen NHK et al. Epigenomic diagnosis and prognosis of Acute Myeloid Leukemia. Nat Commun 2025 Jul 29;16(1):6961. PMID: 40730747

Choudhury SR, Kaushal A, Biswas P, Padilla C et al. Transcriptional rewiring by enhancer methylation in CBFA2T3-GLIS2-driven pediatric acute megakaryoblastic leukemia. Genes Dis 2026 Jan;13(1):101843. PMID: 41141393

DeBruine ZJ, Melcher K, Triche TJ (2021). "High-performance non-negative matrix factorization for large single-cell data." BioRXiv. doi:10.1101/2021.09.01.458620.

## Citation and contact

- **Manuscript citation:** [TO ADD WHEN AVAILABLE]
- **Code version / archived release DOI:** [TO ADD WHEN AVAILABLE]
- **Maintainer/contact:** [OPTIONAL — enter a public professional contact if desired]
- **License:** [DECIDE BEFORE ADDING A LICENSE FILE; absence of a license does not automatically allow code reuse]

## Repository status

This repository is a documented analysis archive under active development, not yet a fully executable end-to-end NMF reproduction. The current public scripts assume existing NMF model outputs. Additional workflow sections and figure legends will be incorporated as they are reviewed.

