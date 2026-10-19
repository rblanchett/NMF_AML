# Supporting figures — public figure gallery

This directory contains **only curated, reviewed, publication-appropriate** figures. Raw generated output folders are not copied wholesale.

## Figure gallery (add links after figures are uploaded)

| Planned ID | Figure | H representation | Suggested output files | Current status |
|---|---|---|---|---|
| FigS01 | NMF group-colored UMAP | Raw + normalized | `FigS01_NMF_Group_raw.png`, `FigS01_NMF_Group_normalized.png` | Not yet added |
| FigS02 | Sex-colored UMAP | Raw + normalized | `FigS02_Sex_raw.png`, `FigS02_Sex_normalized.png` | Not yet added |
| FigS03 | Clinical risk UMAP | Raw + normalized | `FigS03_Risk_raw.png`, `FigS03_Risk_normalized.png` | Not yet added |
| FigS04 | NMF assignment Margin | Raw + normalized | `FigS04_Margin_raw.png`, `FigS04_Margin_normalized.png` | Not yet added |
| FigS05 | Molecular lesion overlays | Specify in caption | `FigS05_Lesions_*.png` | Not yet added |
| FigS06 | Parameter sensitivity / QC | Both | `FigS06_UMAP_QC.png` | Not yet added |

**The figure IDs are suggested editorial labels**, not a claim that these are already your manuscript's final supplementary numbering.

## Figure legend template

**Figure SXX. [Descriptive result-neutral title].** UMAP of pediatric AML samples derived from the [raw / patient-sum-normalized] 23-component NMF H matrix (n = [number of samples actually shown]). Coordinates were computed using `uwot` with [distance metric], `n_neighbors = [ ]`, `min_dist = [ ]`, and random seed [ ]. Points represent individual patients and are colored by [variable with exact definition]. Missing annotations are represented as [Unknown / not displayed, explain]. [Explain panel A / panel B, if relevant]. The UMAP is exploratory; relative distances and visual cluster gaps do not establish statistical significance.

## Publishing workflow

1. Choose a static figure from Script 16's `H_raw/Plots` or `H_normalized/Plots` folder.
2. Verify the sample count, title, risk definitions, plot legends, and annotation coding.
3. Check permissions to redistribute the figure, especially where patient-level combinations could be discerned.
4. Copy a curated PNG **and optionally a PDF** to this directory and rename it with a stable `FigSXX_` prefix.
5. Update this gallery and `supplementary/FIGURE_MANIFEST.md` with output source, parameters, legend, and approval status.
6. Avoid uploading patient-level coordinate CSVs, full plots directory, analysis RDS objects, or logs.

## How to embed a figure in Markdown

Once a reviewed image really exists here, insert a relative link like:

```markdown
![Figure S01 — NMF group UMAP, normalized H](FigS01_NMF_Group_normalized.png)
```

Do not include image links to files that have not yet been uploaded.
