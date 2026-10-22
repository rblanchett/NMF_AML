# NMF Group-Specific Clinical and Molecular Characteristics

## Overview

Categorical clinical, molecular, and cytogenetic variables were compared between each target NMF group and all remaining groups combined. The original analysis used Pearson chi-squared, simulated chi-squared, or Fisher exact tests as reported in the source table; statistical significance was defined as Benjamini–Hochberg-adjusted **q < 0.05**.

**Result:** The supplied final table contains **107 significant NMF-group × variable comparisons across 14 occupied NMF groups**. The 23-component NMF model includes two components without dominant-assignment patients.

### Interpretation of the category columns

The underlying **variable-level** omnibus comparison is significant. The “higher in group” and “lower in group” categories below are categories reported in the source table as having relatively higher or lower prevalence in the target group; they should **not be interpreted as independently FDR-significant category-specific tests** without additional testing. Percentages are calculated relative to the target group or the combined comparison groups, as reported in the supplied results. For `Unknown`, `pending`, `#N/A`, or `#NAME?` categories, interpret as source-data labels or artifacts pending verification, not biological states.

## Group-centric overview

| NMF group | Patients (n) | Significant covariates (n) | Principal overrepresented categories (descriptive) |
|---|---:|---:|---|
| **nmf6** | 37 | 10 | CBFB-MYH11 (`Primary.Fusion`); CBFB-MYH11 (`Cytogenetic.Category.1`); inv(16) (`Primary.Cytogenetic.Code`); Yes (`inv.16.`) |
| **nmf7** | 53 | 14 | Yes (`FLT3.ITD.positive.`); Yes (`WT1.mutation.`); NUP98.NSD1 (`Rare.Fusions`); Normal Karyotype (`Cytogenetic.Category.1`) |
| **nmf8** | 59 | 13 | t(8;21) (`Primary.Cytogenetic.Code`); Yes (`t.8.21.`); RUNX1-RUNX1T1 (`Primary.Fusion`); RUNX1-RUNX1T1 (`Cytogenetic.Category.1`) |
| **nmf9** | 43 | 1 | Other (`Primary.Cytogenetic.Code`) |
| **nmf10** | 33 | 6 | MLL (`Primary.Cytogenetic.Code`); Yes (`MLL`); KMT2A (`Cytogenetic.Category.1`); Yes (`t.9.11..p22.q23.`) |
| **nmf11** | 51 | 11 | Infant (`AgeGroup`); OtherAML (`Cytogenetic.Category.1`); Yes (`M7_AML`); M7 (`FAB`) |
| **nmf14** | 79 | 2 | Male (`Sex`); KMT2A (`Cytogenetic.Category.1`) |
| **nmf15** | 63 | 11 | Yes (`NPM.mutation.`); Normal (`Primary.Cytogenetic.Code`); Normal Karyotype (`Cytogenetic.Category.1`); No (`t.8.21.`) |
| **nmf16** | 11 | 2 | Yes (`MLL`); KMT2A (`Cytogenetic.Category.1`) |
| **nmf17** | 34 | 13 | Infant (`AgeGroup`); OtherAML (`Cytogenetic.Category.1`); Yes (`RAM.phenotype`); Other (`Primary.Cytogenetic.Code`) |
| **nmf18** | 20 | 1 | Female (`Sex`) |
| **nmf20** | 44 | 9 | OtherAML (`Cytogenetic.Category.1`); Other (`Primary.Cytogenetic.Code`); AML with myelodysplasia-related changes (`FAB_or_WHO.Classification`); M1 (`FAB`) |
| **nmf21** | 31 | 8 | Yes (`CEBPA.mutation.`); CEBPA (`SNVs`); Provisional entity: AML with mutated CEBPA (`FAB_or_WHO.Classification`); CEBPA.mutation. (`Mutations.Category`) |
| **nmf22** | 55 | 6 | Yes (`NPM.mutation.`); NPM1 (`SNVs`); Normal Karyotype (`Cytogenetic.Category.1`); Provisional entity: AML with mutated NPM1 (`FAB_or_WHO.Classification`) |

## Detailed significant-variable summaries

Each row lists the **reported** highest-relative-prevalence category in the target group and the reported relatively underrepresented category. `q` is the adjusted p-value for the **entire covariate**. The listed category-specific percentages are descriptive.

### nmf6 — n = 37 (vs 817 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `c.Kit.Mutation.Exon.8` | Yes (29.7% vs 2.6%) | No (70.3% vs 97.4%) | 8.71e-07 |
| `CBL.Mutation` | Del_exons_8_9 (16.2% vs 1.2%) | `Unknown` (source label) (67.6% vs 89.1%) | 0.00777 |
| `Cytogenetic.Category.1` | CBFB-MYH11 (86.5% vs 6.7%) | KMT2A (0.0% vs 25.5%) | 0.00777 |
| `FAB_or_WHO.Classification` | AML with inv(16)(p13q22) or t(16;16)(p13;q22); CBFB-MYH11 (81.1% vs 6.0%) | AML, not otherwise categorized:  Acute monoblastic/acute monocytic leukemia (0.0% vs 15.3%) | 0.00777 |
| `Mutations.Category` | c.Kit.Ex8 (21.6% vs 1.6%) | OtherAML (29.7% vs 56.1%) | 0.00777 |
| `Primary.Cytogenetic.Code` | inv(16) (83.8% vs 6.6%) | MLL (0.0% vs 24.2%) | 0.00777 |
| `RAS.Gene` | NRAS (43.2% vs 21.3%) | Not specified (40.5% vs 72.9%) | 0.00777 |
| `inv.16.` | Yes (83.8% vs 6.6%) | No (13.5% vs 92.5%) | 0.00777 |
| `Primary.Fusion` | CBFB-MYH11 (86.5% vs 6.7%) | Not specified (10.8% vs 28.3%) | 0.014 |
| `Primary.CNV` | CBL_deletion (24.3% vs 3.7%) | No.Relevant.CNV (62.2% vs 74.2%) | 0.029 |

### nmf7 — n = 53 (vs 801 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `WT1.mutation.` | Yes (37.7% vs 4.5%) | No (62.3% vs 95.5%) | 2.24e-10 |
| `Cytogenetic.Category.1` | Normal Karyotype (52.8% vs 21.1%) | KMT2A (5.7% vs 25.6%) | 5.07e-08 |
| `MRD.at.end.of.course.2` | `Unknown` (source label) (37.7% vs 20.7%) | No (37.7% vs 70.4%) | 7.89e-05 |
| `AgeGroup` | AYA (37.7% vs 25.2%) | Infant (1.9% vs 26.1%) | 0.00777 |
| `FAB_or_WHO.Classification` | AML, not otherwise categorized:  AML with maturation (22.6% vs 6.5%) | AML with t(8;21)(q22;q22); RUNX1-RUNX1T1 (0.0% vs 14.7%) | 0.00777 |
| `FLT3.ITD.positive.` | Yes (64.2% vs 8.2%) | No (28.3% vs 87.5%) | 0.00777 |
| `MRD.at.end.of.course.1` | Yes (50.9% vs 24.6%) | No (39.6% vs 70.2%) | 0.00777 |
| `Mutations.Category` | FLT3.ITD (30.2% vs 2.7%) | OtherAML (13.2% vs 57.7%) | 0.00777 |
| `Primary.Cytogenetic.Code` | Normal (52.8% vs 22.2%) | MLL (3.8% vs 24.5%) | 0.00777 |
| `RAM.phenotype` | `Unknown` (source label) (60.4% vs 6.6%) | No (39.6% vs 90.8%) | 0.00777 |
| `Rare.Fusions` | NUP98.NSD1 (35.8% vs 3.1%) | OtherAML (60.4% vs 89.4%) | 0.00777 |
| `SNVs` | FLT3-ITD (34.0% vs 3.2%) | OtherAML (17.0% vs 79.5%) | 0.00777 |
| `Treatment.Arm` | Arm C (26.4% vs 4.2%) | Arm B (26.4% vs 49.8%) | 0.00777 |
| `FAB` | M2 (22.6% vs 6.5%) | `pending` (source label) (22.6% vs 57.1%) | 0.014 |

### nmf8 — n = 59 (vs 795 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `Cytogenetic.Category.1` | RUNX1-RUNX1T1 (96.6% vs 9.1%) | OtherAML (0.0% vs 29.3%) | 1.35e-67 |
| `c.Kit.Mutation.Exon.17` | Yes (35.6% vs 4.3%) | No (64.4% vs 95.7%) | 1.41e-10 |
| `AgeGroup` | Child (79.7% vs 47.2%) | Infant (0.0% vs 26.4%) | 4.26e-05 |
| `FAB_or_WHO.Classification` | AML with t(8;21)(q22;q22); RUNX1-RUNX1T1 (88.1% vs 8.3%) | AML, not otherwise categorized:  Acute monoblastic/acute monocytic leukemia (0.0% vs 15.7%) | 0.00777 |
| `Minus.X` | Yes (20.3% vs 1.3%) | No (79.7% vs 97.7%) | 0.00777 |
| `Minus.Y` | Yes (28.8% vs 3.8%) | No (71.2% vs 95.2%) | 0.00777 |
| `Mutations.Category` | c.Kit.Ex17 (28.8% vs 2.6%) | OtherAML (39.0% vs 56.1%) | 0.00777 |
| `Primary.Cytogenetic.Code` | t(8;21) (98.3% vs 8.6%) | Other (0.0% vs 29.1%) | 0.00777 |
| `Primary.Fusion` | RUNX1-RUNX1T1 (96.6% vs 9.1%) | Not specified (0.0% vs 29.6%) | 0.00777 |
| `t.8.21.` | Yes (98.3% vs 8.6%) | No (1.7% vs 90.4%) | 0.00777 |
| `CSF3R` | c.2252insC/50 (6.8% vs 0.4%) | `Unknown` (source label) (91.5% vs 98.5%) | 0.014 |
| `MLL` | No (100.0% vs 74.0%) | Yes (0.0% vs 25.0%) | 0.0189 |
| `c.Kit.Mutation.Exon.8` | Yes (11.9% vs 3.1%) | No (88.1% vs 96.9%) | 0.0465 |

### nmf9 — n = 43 (vs 811 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `Primary.Cytogenetic.Code` | Other (44.2% vs 26.1%) | t(8;21) (0.0% vs 15.5%) | 0.0426 |

### nmf10 — n = 33 (vs 821 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `Cytogenetic.Category.1` | KMT2A (75.8% vs 22.3%) | RUNX1-RUNX1T1 (0.0% vs 15.7%) | 0.00777 |
| `MLL` | Yes (75.8% vs 21.2%) | No (24.2% vs 77.8%) | 0.00777 |
| `Primary.Cytogenetic.Code` | MLL (75.8% vs 21.1%) | Normal (9.1% vs 24.7%) | 0.00777 |
| `t.9.11..p22.q23.` | Yes (48.5% vs 7.2%) | No (51.5% vs 91.8%) | 0.00777 |
| `FAB_or_WHO.Classification` | AML with t(9;11)(p22;q23); MLLT3-MLL (36.4% vs 6.7%) | AML with t(8;21)(q22;q22); RUNX1-RUNX1T1 (0.0% vs 14.4%) | 0.014 |
| `Additional.Fusions.CNV` | AMBRA1-MLLT10; MLLT10-AMBRA1; ABTB2-KMT2A (3.0% vs 0.0%) | Not specified (78.8% vs 94.0%) | 0.0333 |

### nmf11 — n = 51 (vs 803 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `AgeGroup` | Infant (64.7% vs 22.0%) | Child (27.5% vs 50.8%) | 3.12e-09 |
| `Cyto.Fusion.Molecular.Risk` | High (49.0% vs 31.4%) | Low (11.8% vs 41.1%) | 0.00645 |
| `Cytogenetic.Category.2` | M7 (25.5% vs 3.7%) | OtherAML (54.9% vs 85.3%) | 0.00777 |
| `FAB` | M7 (25.5% vs 3.7%) | `pending` (source label) (39.2% vs 55.9%) | 0.00777 |
| `FAB_or_WHO.Classification` | AML, not otherwise categorized:  Acute megakaryoblastic leukemia (19.6% vs 3.4%) | AML with inv(16)(p13q22) or t(16;16)(p13;q22); CBFB-MYH11 (0.0% vs 9.8%) | 0.00777 |
| `M7_AML` | Yes (25.5% vs 3.7%) | No (74.5% vs 96.1%) | 0.00777 |
| `Rare.Fusions` | NUP98.KDM5A (11.8% vs 1.0%) | OtherAML (80.4% vs 88.0%) | 0.00777 |
| `Was.the.patient.diagnosed.with.non.CNS.extramedullary.disease.` | Yes (29.4% vs 13.7%) | No (68.6% vs 86.3%) | 0.00777 |
| `Primary.CNV` | trisomy21 (17.6% vs 2.5%) | No.Relevant.CNV (60.8% vs 74.5%) | 0.014 |
| `Cytogenetic.Category.1` | OtherAML (49.0% vs 25.9%) | CBFB-MYH11 (0.0% vs 10.8%) | 0.0149 |
| `Primary.Fusion` | NUP98-KDM5A (11.8% vs 1.0%) | CBFB-MYH11 (0.0% vs 10.8%) | 0.029 |

### nmf14 — n = 79 (vs 775 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `Sex` | Male (72.2% vs 51.2%) | Female (27.8% vs 48.8%) | 0.00777 |
| `Cytogenetic.Category.1` | KMT2A (34.2% vs 23.4%) | Normal Karyotype (12.7% vs 24.1%) | 0.0298 |

### nmf15 — n = 63 (vs 791 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `NPM.mutation.` | Yes (33.3% vs 7.3%) | No (66.7% vs 92.7%) | 5.16e-10 |
| `Cytogenetic.Category.1` | Normal Karyotype (39.7% vs 21.7%) | RUNX1-RUNX1T1 (0.0% vs 16.3%) | 0.000926 |
| `FAB_or_WHO.Classification` | Provisional entity: AML with mutated NPM1 (17.5% vs 3.3%) | AML with t(8;21)(q22;q22); RUNX1-RUNX1T1 (0.0% vs 14.9%) | 0.00777 |
| `Mutations.Category` | NPM.mutation. (15.9% vs 3.3%) | FLT3.PM (1.6% vs 7.5%) | 0.00777 |
| `Primary.Fusion` | Not specified (42.9% vs 26.3%) | RUNX1-RUNX1T1 (0.0% vs 16.3%) | 0.00777 |
| `SNVs` | NPM1 (20.6% vs 4.9%) | OtherAML (54.0% vs 77.4%) | 0.00777 |
| `t.6.9.` | Yes (9.5% vs 0.8%) | No (90.5% vs 98.2%) | 0.00777 |
| `Primary.Cytogenetic.Code` | Normal (41.3% vs 22.8%) | t(8;21) (0.0% vs 15.9%) | 0.014 |
| `Rare.Fusions` | DEK.NUP214 (9.5% vs 0.8%) | CBFA2T3.GLIS2 (0.0% vs 2.4%) | 0.0189 |
| `t.8.21.` | No (100.0% vs 83.1%) | Yes (0.0% vs 15.9%) | 0.0383 |
| `FLT3.ITD.positive.` | Yes (23.8% vs 10.7%) | No (68.3% vs 85.1%) | 0.0426 |

### nmf16 — n = 11 (vs 843 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `Cytogenetic.Category.1` | KMT2A (81.8% vs 23.6%) | OtherAML (9.1% vs 27.5%) | 0.029 |
| `MLL` | Yes (81.8% vs 22.5%) | No (18.2% vs 76.5%) | 0.0465 |

### nmf17 — n = 34 (vs 820 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `AgeGroup` | Infant (76.5% vs 22.4%) | Child (17.6% vs 50.7%) | 4.92e-10 |
| `Cytogenetic.Category.1` | OtherAML (76.5% vs 25.2%) | KMT2A (8.8% vs 25.0%) | 0.00777 |
| `Cytogenetic.Code.Other` | +19,+21,+der(22)t(1;22)(q21;q13),+8,+10 (2.9% vs 0.0%) | `Unknown` (source label) (38.2% vs 63.4%) | 0.00777 |
| `FAB` | M7 (29.4% vs 4.0%) | `pending` (source label) (26.5% vs 56.1%) | 0.00777 |
| `FAB_or_WHO.Classification` | AML, not otherwise categorized:  Acute megakaryoblastic leukemia (29.4% vs 3.3%) | AML with t(8;21)(q22;q22); RUNX1-RUNX1T1 (0.0% vs 14.4%) | 0.00777 |
| `M7_AML` | Yes (29.4% vs 4.0%) | No (70.6% vs 95.9%) | 0.00777 |
| `Primary.Cytogenetic.Code` | Other (67.6% vs 25.4%) | t(8;21) (0.0% vs 15.4%) | 0.00777 |
| `Primary.Fusion` | CBFA2T3-GLIS2 (41.2% vs 0.6%) | Not specified (5.9% vs 28.4%) | 0.00777 |
| `RAM.phenotype` | Yes (44.1% vs 0.7%) | No (52.9% vs 89.0%) | 0.00777 |
| `Rare.Fusions` | CBFA2T3.GLIS2 (41.2% vs 0.6%) | OtherAML (50.0% vs 89.1%) | 0.00777 |
| `Cytogenetic.Category.2` | M7 (29.4% vs 4.0%) | OtherAML (67.6% vs 84.1%) | 0.014 |
| `MRD.at.end.of.course.1` | Yes (58.8% vs 24.9%) | No (38.2% vs 69.5%) | 0.0189 |
| `t.8.21.` | No (94.1% vs 83.9%) | Yes (0.0% vs 15.4%) | 0.0333 |

### nmf18 — n = 20 (vs 834 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `Sex` | Female (100.0% vs 45.6%) | Male (0.0% vs 54.4%) | 7.49e-05 |

### nmf20 — n = 44 (vs 810 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `Cytogenetic.Category.1` | OtherAML (61.4% vs 25.4%) | RUNX1-RUNX1T1 (0.0% vs 15.9%) | 0.00777 |
| `Cytogenetic.Category.2` | del5q (9.1% vs 0.6%) | M7 (0.0% vs 5.3%) | 0.00777 |
| `Primary.Fusion` | PICALM-MLLT10 (15.9% vs 0.0%) | RUNX1-RUNX1T1 (0.0% vs 15.9%) | 0.00777 |
| `del5q` | Yes (9.1% vs 0.6%) | No (88.6% vs 98.5%) | 0.00777 |
| `FAB_or_WHO.Classification` | AML with myelodysplasia-related changes (25.0% vs 6.8%) | AML with t(8;21)(q22;q22); RUNX1-RUNX1T1 (0.0% vs 14.6%) | 0.014 |
| `Cytogenetic.Code.Other` | ~4n [other abnormalities could not be ruled-out] (2.3% vs 0.0%) | `Unknown` (source label) (50.0% vs 63.1%) | 0.0247 |
| `FAB` | M1 (22.7% vs 6.0%) | `pending` (source label) (36.4% vs 55.9%) | 0.029 |
| `Primary.Cytogenetic.Code` | Other (52.3% vs 25.7%) | t(8;21) (0.0% vs 15.6%) | 0.029 |
| `Primary.CNV` | del5q (9.1% vs 0.6%) | `#NAME?` (source label) (0.0% vs 4.8%) | 0.0333 |

### nmf21 — n = 31 (vs 823 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `CEBPA.mutation.` | Yes (61.3% vs 2.8%) | No (38.7% vs 97.2%) | 1.25e-17 |
| `CEBPA.Double.Single.Allelic.Status` | Double Mutant (45.2% vs 2.6%) | Not specified (38.7% vs 97.2%) | 0.00777 |
| `FAB_or_WHO.Classification` | Provisional entity: AML with mutated CEBPA (45.2% vs 1.7%) | AML, not otherwise categorized:  Acute monoblastic/acute monocytic leukemia (6.5% vs 14.9%) | 0.00777 |
| `GATA2` | R330L (9.7% vs 0.0%) | `Unknown` (source label) (83.9% vs 99.1%) | 0.00777 |
| `Mutations.Category` | CEBPA.mutation. (45.2% vs 2.3%) | OtherAML (22.6% vs 56.1%) | 0.00777 |
| `SNVs` | CEBPA (48.4% vs 2.6%) | OtherAML (32.3% vs 77.3%) | 0.00777 |
| `Cytogenetic.Category.1` | Normal Karyotype (58.1% vs 21.7%) | KMT2A (6.5% vs 25.0%) | 0.0189 |
| `Primary.Cytogenetic.Code` | Normal (58.1% vs 22.8%) | MLL (6.5% vs 23.8%) | 0.029 |

### nmf22 — n = 55 (vs 799 other patients)

| Significant variable | Higher-prevalence category in target (target vs others) | Lower-prevalence category in target (target vs others) | BH q |
|---|---|---|---:|
| `NPM.mutation.` | Yes (36.4% vs 7.4%) | No (63.6% vs 92.6%) | 6.34e-11 |
| `Mutations.Category` | NPM.mutation. (16.4% vs 3.4%) | OtherAML (41.8% vs 55.8%) | 0.00777 |
| `SNVs` | NPM1 (27.3% vs 4.6%) | OtherAML (58.2% vs 76.8%) | 0.00777 |
| `Cytogenetic.Category.1` | Normal Karyotype (43.6% vs 21.7%) | OtherAML (18.2% vs 27.9%) | 0.00925 |
| `FAB_or_WHO.Classification` | Provisional entity: AML with mutated NPM1 (21.8% vs 3.1%) | AML with inv(16)(p13q22) or t(16;16)(p13;q22); CBFB-MYH11 (1.8% vs 9.8%) | 0.014 |
| `AgeGroup` | AYA (41.8% vs 24.9%) | Infant (7.3% vs 25.8%) | 0.0201 |

## Notes for the interactive UMAP explorer

- Mode 1: preserve the current dropdown showing all configured annotations.
- Mode 2: select one of the 14 significant NMF groups and one of its significant *variables*. To distinguish presence versus absence, the selected category must be specified and mapped against the actual patient-level annotation values.
- Mode 3: select a characteristic, then display which NMF groups have a significant **variable-level** association and the category-specific descriptive directions. Group highlighting should not imply independently significant category-level tests.
- For four-way coloring, classify patients as (inside group, characteristic present), (outside group, present), (inside group, absent), (outside group, absent). Preserve missing/unknown annotations separately; they must not silently count as absent.

## Source and audit notes

- Source: user-provided final significant-characteristics TSV (107 rows).
- The table contains some data labels such as `Unknown`, `pending`, `#N/A`, and `#NAME?`. These are carried through rather than silently repaired.
- This summary does not perform new hypothesis tests and does not independently verify category-specific enrichment significance.
- Repeated representations of a biological event across multiple covariates (e.g., fusion and cytogenetic code) are retained as separate source variables.
