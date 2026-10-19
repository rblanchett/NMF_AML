# Supplementary figure and table manifest

Use this file to trace every scientific figure to the data representation, code version, and finalized legend. This is particularly useful for a supplement with many UMAP covariate overlays.

| Figure ID | Intended content | Raw / normalized | Source script | Plot input / parameter record | Public file | Legend checked? | Public release reviewed? |
|---|---|---|---|---|---|---|---|
| FigS01 | NMF group UMAP | Both | Script 16 | `UMAP_Primary_Objects.rds` and config | Pending | No | No |
| FigS02 | Sex UMAP | Both | Script 16 | Same primary embeddings | Pending | No | No |
| FigS03 | Clinical risk UMAP | Both | Script 16 | Check exact risk scheme | Pending | No | No |
| FigS04 | Margin UMAP | Both | Script 16 | Margin and relative Margin | Pending | No | No |
| FigS05 | Lesion overlays | To decide | Script 16 | Record lesion coding per column | Pending | No | No |
| FigS06 | UMAP sensitivity/QC | Both | Script 16 | `All_UMAP_Parameter_QC.csv` | Pending | No | No |

## Traceability metadata per finalized figure

- Git commit hash and script filename:
- UMAP representation, input rows/features, distance metric:
- UMAP primary or sensitivity parameters:
- Covariate field name and source data dictionary:
- Missing data, excluded records, plotting subset:
- Exact source output path retained **locally** (avoid personal paths in public manifest):
- Final repository file (PNG / PDF):
- Figure legend:
- Reviewer/date and data-sharing decision:
- Corresponding manuscript figure or Supplementary Table number:

## Supplementary table suggestions

The following can be reported as **aggregate summaries** after disclosure review, not as identifiable patient tables: group counts by sex/risk, counts of lesion-positive patients per NMF group, UMAP QC comparisons, and NMF group summary characteristics. Statistical tests and their parameters should be linked to the actual survival/association scripts when those scripts are added.
