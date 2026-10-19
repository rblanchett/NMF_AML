# Script 16: detailed UMAP methods and execution

## Scope

This describes Script 16 and the optional explorer addenda **as implemented in the current public code**. It is suitable as a methods draft, but numerical results, package versions, and source accession details require verification from the executed analysis and manuscript records.

## NMF representation

Let `V` be the feature-by-patient methylation beta matrix and `V ≈ W H` its previously fitted nonnegative matrix factorization. For the combined-sex pediatric AML cohort, H is a 23-by-854 nonnegative matrix (component-by-patient). The Script 16 UMAP input is its transpose, 854-by-23, containing one vector per patient. The model is **not refitted** in Script 16. Each patient has a discrete group assignment determined elsewhere from the dominant H component, but the UMAP uses **all 23 component weights**.

## H representations

- **Raw:** original H component weights in the fitted model.
- **Patient-sum-normalized:** for patient `i`, `h_norm[i,j] = h[i,j] / sum_j h[i,j]`. Script 16 stops rather than silently dropping patients if a normalization denominator is zero.

Raw weights and normalized weights are saved independently for downstream interpretation. These representations are not interchangeable; NMF factor scales can affect the geometry of H.

## UMAP algorithm and settings

`uwot::umap()` is applied separately to each representation, with the **current configured primary values**:

| Parameter | Raw H | Normalized H |
|---|---:|---:|
| Distance metric | Euclidean | Cosine |
| `n_neighbors` | 30 | 30 |
| `min_dist` | 0.30 | 0.30 |
| `spread` | 1.0 | 1.0 |
| `n_epochs` | 500 | 500 |
| Threads | 1 | 1 |
| Random seed | 190 | 190 |
| Output dimensions | 2 | 2 |

A sensitivity grid additionally compares `n_neighbors = 15, 30, 50` and `min_dist = 0.1, 0.3, 0.5`, subject to the script's configured maximum number of runs. A mean 15-nearest-neighbor overlap statistic compares neighborhood identities in the input 23-dimensional representation with the corresponding two-dimensional UMAP. The QC is diagnostic, not a statistical hypothesis test.

## Annotating the embedding

Clinical and molecular variables are merged by the source internal ID (`USI`) after checking the input alignment. Each patient retains the same coordinates across plots within an H representation. Plots are recolored by annotations without refitting UMAP. Missing annotations are kept as Unknown (where appropriate); missing molecular annotations must not be interpreted as absence.

**NMF assignment Margin** is the difference between the largest and second-largest original H component weights for a patient. The optional *relative Margin* is the difference between the largest and second-largest patient-sum-normalized weights. Neither is a calibrated class membership probability.

Sex is represented with a blue/pink categorical palette; molecular lesions can be colored Present (red), Absent (gray), or Unknown; three-category risk schemes may be displayed as green/yellow/red only after checking actual category labels. Continuous margins use a color gradient.

## Visualization and outputs

`ggplot2` generates PDF/PNG figures, and `plotly::ggplotly()` is used for interactive per-covariate outputs. Script 16 records run information, warnings, parameter choices, session information, QC, and output manifest. The Section 13 addendum is an append-only unified explorer option; the standalone Section 14 exports a reviewed GitHub Pages HTML explorer from **saved Script 16 results**, without recomputing NMF/UMAP.

## Public exporter precautions

Section 14 substitutes generic display labels for USI and blocks configured survival variables from its public field list. It does **not** turn individual-level methylation-derived coordinates and clinical features into provably anonymous data. A publication decision still requires validation of the source data access conditions and the exact HTML payload. Do not commit its local private audit/log files.

## Running Script 16

1. Install the R dependencies documented in the script in your own R environment.
2. Keep input/output files outside the Git repository.
3. Define `NMF_PROJECT_DIR` (a local private project root) in your shell environment or local `.Renviron` file; never commit your real `.Renviron`.
4. Provide the input paths expected by the public `CONFIG` block (`inputs/h_matrix.rds` or `inputs/nmf_model.rds`, `inputs/nmf_assignments.csv`, and `inputs/clinical_annotations.csv`) or adapt them locally.
5. Run `Rscript scripts/16_NMF_H_UMAP_Interactive_Exploration_PUBLIC.R`.
6. Check the log and QC under the private output folder before selecting figures for release.
7. For the standalone public exporter, configure `SCRIPT16_OUTPUT_DIR` and run `Rscript scripts/16_Section_14_Public_GitHub_UMAP_Exporter.R` **only after completing the script's redistribution review**.

## Manuscript-ready methods draft (review before reuse)

> The previously estimated 23-component NMF patient-loading matrix (H) was used to visualize relationships among 854 pediatric AML samples. Each sample was represented by its full vector of 23 nonnegative component weights. Two H representations were evaluated: the original H weights and a patient-normalized representation in which the component weights of each sample were divided by their sum. Uniform Manifold Approximation and Projection (UMAP) was performed independently on the two representations using the `uwot` R package with 30 neighbors, a minimum distance of 0.30, a spread of 1.0, 500 training epochs, and random seed 190; Euclidean distance was used for raw H and cosine distance for normalized H. Sensitivity analyses evaluated alternative neighbor counts (15, 30, and 50) and minimum-distance values (0.1, 0.3, and 0.5). Local neighborhood preservation was assessed using mean 15-nearest-neighbor overlap between each H representation and its two-dimensional embedding. Clinical and molecular metadata were joined to each sample after embedding and used solely for visualization. NMF assignment margins were examined as descriptive indicators of dominance of the leading NMF component. UMAPs were considered exploratory visualizations and were not used to establish statistical group separation or survival associations.

Check actual executed settings and R package versions before copying this paragraph into a manuscript.
