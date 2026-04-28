# Lesley TB — Mother–Infant Paired BCG Stimulation Study

Immunological analysis of a **paired mother–infant cohort** profiling responses to **BCG vaccination** and other *Mycobacterium tuberculosis* (Mtb) antigens. Samples are stratified by **maternal IGRA status** and **infant IGRA status**, with infants sampled longitudinally at **12 weeks** and **44 weeks** of life.

The repository contains two complementary projects sharing the same cohort:

1. **Flow cytometry + plasma biomarkers** (`Scripts/`) — the primary analysis: intracellular cytokine staining (ICS), proliferation gating, and a multiplex plasma cytokine/biomarker panel.
2. **TB CITE-seq** (`TB_CITEseq/`) — single-cell RNA + surface protein (ADT) + αβ TCR (VDJ) profiling on a subset of the same samples, with isotype-corrected DSB-normalized ADTs.

---

## Cohort and stimulations

| | |
|---|---|
| **Subjects** | Mothers and their infants (paired via `Maternal PID` / `FamilyID`) |
| **Maternal timepoints** | Entry, Dx |
| **Infant timepoints** | 12 weeks, 44 weeks (longitudinally paired) |
| **In vitro stimulations** | **MED** (medium control, baseline), **BCG** (vaccine), **DosR** (Mtb dormancy regulon antigens), **E6C10** (ESAT-6 / CFP-10, Mtb-specific RD1 antigens), **GAG** (HIV Gag protein) |
| **IGRA stratification** | Maternal IGRA⁺/⁻ and Infant IGRA⁺/⁻ |

The stimulation panel covers complementary Mtb biology: BCG models vaccine-induced responses, E6C10 reads out Mtb-specific (RD1) effector immunity, and DosR captures the dormancy/latency-associated repertoire. GAG provides an HIV-relevant comparator antigen.

---

## Project A — Flow cytometry + plasma biomarkers

### Inputs

```
Raw_Data/
├── TB_Flow_Cytokine DATA_102725.xlsx        # ICS / proliferation flow data (FlowJo export)
└── TB Flow and Plasma Biomarker data.xlsx   # Combined workbook including Plasma Biomarker sheet
```

### Pipeline (`Scripts/`)

#### `1.Data_Cleaning_Norm.R` — Data cleaning + normalization

Reads the FlowJo export, harmonizes typed-missing variants (`NA`, `N/A`, `--`, etc.), pivots to long form keyed by `GatePath`, and re-parses the gating hierarchy. Then performs two normalization passes:

- **MED-baseline subtraction.** Within each `PID × Time_Point × GatePath`, the medium-control mean is subtracted from each stimulated value to give `Value_BaselineAdj` (`ΔStim = Stim − MED`).
- **Logit-scale batch normalization.** Using a healthy-control sample at MED (`HC-SP @ MED`) as the per-batch reference, raw absolute percentages are shifted on the **logit scale** (additive) with a 2.0-cap on absolute shift (~7.4× odds), then back-transformed to percent. ΔStim values are kept on the linear percentage scale (because they can be negative). This avoids the multiplicative blow-up that scaling raw frequencies produces near zero.

Below-detection (`< 0.01%`) values are flagged separately so downstream analyses can decide whether to treat them as zero or as missing. Outputs are written as both `.rds` and `.csv` to `saved_R_dat/`:

```
TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.{rds,csv}
```

#### `2. Paired_Test_ProlifvsNonProlif.R` — Paired and group comparisons

Restricts to TB stimulations (`BCG`, `DosR`, `E6C10`, `GAG`) and `Freq. of Parent` cytokine gates, builds composite phenotype names from the T-level gate columns (e.g. `CCR6+ CXCR3-`, `CD45RO+ CD27+`), and runs four pre-defined comparisons:

| | Comparison | Test |
|---|---|---|
| **Comp 1** | Infants paired 12wks vs 44wks (by stim) | Paired Wilcoxon (within infant) |
| **Comp 2** | Infants by maternal IGRA⁺/⁻ (per timepoint × stim) | Mann-Whitney |
| **Comp 3** | Infants by infant IGRA⁺/⁻ (per timepoint × stim) | Mann-Whitney |
| **Comp 4** | Mothers across groups (by stim) | Kruskal–Wallis + pairwise Wilcoxon |

A separate **Proliferating vs Non-Proliferating** paired Wilcoxon (with rank-biserial effect size) is run within each subset across all readouts; results live one-folder-per-readout under `Results/Prolif_vs_NonProlif/` (e.g. `CCR6+ CXCR3+/`, `CD107a/`, `CD45RO+ CD27-/`, …) along with a top-level `_ALL_readouts_summary.csv`.

#### `3.Plasma_Biomarker.R` — Plasma cytokine panel

Reads the **Plasma Biomarker** sheet of the Excel workbook, parses qualitative flags (`UND`, `OVER`, `NA`) out of analyte values, log-transforms and scales within timepoint, and exports both wide and long forms:

```
saved_R_dat/Plasma_Biomarker_long_scaled.{rds,csv}
saved_R_dat/Plasma_Biomarker_wide_scaled.{rds,csv}
```

It then runs the same four comparison structure as the flow analysis, with results in `Plasma_Results/`:

```
Plasma_Results/
├── Infants_12vs44_paired/   stats_*.csv + plots/
├── Infants_byInfantIGRA/    stats_*.csv + plots/
├── Infants_byMaternalIGRA/  stats_*.csv + plots/
└── Mothers_byMaternalIGRA/  stats_*.csv + plots/
```

#### `4.Lasso_maternal_Predictors.R` — Maternal predictors of infant responses

Builds a master analysis-ready table joining maternal flow + plasma features to the infant 12-week proliferation outcome (CD8⁺), then fits **LASSO logistic regression** (`glmnet`) with **ROC evaluation** (`pROC`) separately per stimulation antigen / data layer:

```
Model_Results/
├── DosR/      coefficients_DosR.csv, summary_DosR.csv, roc_DosR.png
├── E6C10/     coefficients_E6C10.csv, summary_E6C10.csv, roc_E6C10.png
├── GAG/       coefficients_GAG.csv, summary_GAG.csv, roc_GAG.png
├── Plasma/    coefficients_Plasma.csv, summary_Plasma.csv, roc_Plasma.png
└── summary_all_layers.csv
```

The question this answers: *which maternal flow-cytometric or plasma-biomarker features predict infant CD8⁺ proliferative response to a given antigen at 12 weeks?*

### Per-readout results (`Results/Prolif_vs_NonProlif/`)

One folder per phenotypic gate, with a paired-differences table and a Wilcoxon stats CSV. Phenotypes covered include: T-helper subset markers (`CCR6± CXCR3±`, `CCR4+`, `CD161+`), memory subsets (`CD45RO± CD27±`), T-follicular helper marker (`CXCR5+`), and degranulation marker (`CD107a`).

---

## Project B — TB CITE-seq (`TB_CITEseq/`)

A multimodal single-cell pipeline on Cell Ranger `multi` outputs (10x Genomics 5' chemistry, with Antibody Capture and VDJ-T libraries).

### Pipeline (`TB_CITEseq/Scripts/`)

| | Script | What it does |
|---|---|---|
| 1 | `1.dsb_norm.R` | Reads raw + filtered Cell Ranger h5 matrices, applies **DSB normalization** to the ADT assay (empty droplets + isotype controls), saves merged Seurat as `.qs2` |
| 2 | `2.Basic_QC.R` | RNA QC visualizations and filtering (counts, features, mito, ribo) |
| 3 | `3.Cell_Cycle_Doublet_Removal.R` | `CellCycleScoring` + per-sample doublet detection |
| 4 | `4.Isotype_Correction.R` | Per-marker thresholding using `Isotype.xlsx` panel-defined isotype controls |
| 5 | `5.Integration_WNN.R` | Cross-sample integration + Weighted Nearest Neighbor (RNA + ADT) joint embedding |
| 6 | `6.Pre_Annotation_Plots.R` | UMAP overviews and pre-annotation marker plots (`UMAP_overview_res0.6.png`, etc.) |
| 7 | `7.VDJ.R` | scRepertoire-based αβ TCR clonotype analysis |

### Outputs

```
TB_CITEseq/
├── Isotype.xlsx                # Panel definition: marker → isotype-control mapping
├── Scripts/                    # 7 numbered R scripts (above)
├── QC/                         # RNA QC outputs from script 2
├── Integration/                # WNN integration outputs from script 5
├── Pre_Annotation_Plots/       # UMAPs / marker plots from script 6
└── VDJ/                        # scRepertoire clonotype outputs from script 7
```

---

## Repository structure

```
Lesley_TB_Mother_Infant_Paired_BCG_Stim/
├── Raw_Data/
│   ├── TB_Flow_Cytokine DATA_102725.xlsx
│   └── TB Flow and Plasma Biomarker data.xlsx
├── Scripts/                                # Project A pipeline (flow + plasma)
│   ├── 1.Data_Cleaning_Norm.R
│   ├── 2. Paired_Test_ProlifvsNonProlif.R
│   ├── 3.Plasma_Biomarker.R
│   └── 4.Lasso_maternal_Predictors.R
├── saved_R_dat/                            # Cached preprocessed objects (RDS + CSV)
│   ├── TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.{rds,csv}
│   ├── Plasma_Biomarker_long_scaled.{rds,csv}
│   ├── Plasma_Biomarker_wide_scaled.{rds,csv}
│   └── master_analysis_ready.{rds,csv}
├── Results/
│   ├── Prolif_Only/                        # Comp1–4 proliferation-only summaries
│   │   ├── Comparisons/                    # Per-comparison stats CSVs
│   │   └── Plots/                          # Slopes, violins, box plots, summary heatmaps
│   └── Prolif_vs_NonProlif/                # One subfolder per phenotypic gate
│       └── _ALL_readouts_summary.csv
├── Plasma_Results/                         # Plasma biomarker comparisons (4 groups)
├── Model_Results/                          # LASSO coefficients + ROC per layer
│   ├── DosR/, E6C10/, GAG/, Plasma/
│   └── summary_all_layers.csv
└── TB_CITEseq/                             # Project B (CITE-seq + VDJ)
    ├── Isotype.xlsx
    ├── Scripts/                            # 7 numbered R scripts
    ├── QC/, Integration/, Pre_Annotation_Plots/, VDJ/
    └── saved_R_data/                       # Intermediate Seurat .qs2 objects
```

---

## Dependencies

### Project A (flow + plasma + LASSO)
- [`readxl`](https://CRAN.R-project.org/package=readxl), [`writexl`](https://CRAN.R-project.org/package=writexl), [`readr`](https://CRAN.R-project.org/package=readr) — I/O
- [`dplyr`](https://CRAN.R-project.org/package=dplyr), [`tidyr`](https://CRAN.R-project.org/package=tidyr), [`stringr`](https://CRAN.R-project.org/package=stringr), [`purrr`](https://CRAN.R-project.org/package=purrr) — wrangling
- [`ggplot2`](https://ggplot2.tidyverse.org/), [`patchwork`](https://patchwork.data-imaginist.com/), [`ggtext`](https://CRAN.R-project.org/package=ggtext), [`ragg`](https://CRAN.R-project.org/package=ragg), `glue` — plotting
- [`broom`](https://CRAN.R-project.org/package=broom) — tidy stats
- [`glmnet`](https://CRAN.R-project.org/package=glmnet), [`pROC`](https://CRAN.R-project.org/package=pROC) — LASSO + ROC
- [`digest`](https://CRAN.R-project.org/package=digest) — reproducible IDs

### Project B (CITE-seq)
- [`Seurat`](https://satijalab.org/seurat/) v5
- [`dsb`](https://CRAN.R-project.org/package=dsb) — ADT denoising
- [`scRepertoire`](https://www.borch.dev/uploads/screpertoire/) — TCR analysis
- [`scCustomize`](https://samuel-marsh.github.io/scCustomize/) — single-cell visualization
- [`qs2`](https://CRAN.R-project.org/package=qs2) — fast object serialization
- `tidyverse`, `Matrix`, `hdf5r`, `data.table`, `cowplot`, `biomaRt`

---

## Reproducing the analysis

1. Clone the repo and open `Lesley_TB_Mother_Infant_Paired_BCG_Stim.Rproj` in RStudio.
2. Update the path variables near the top of each script (the scripts currently switch between Linux paths under `/home/akshay-iyer/Documents/...` and Windows paths under `C:/Users/ammas/Documents/...`).
3. Run scripts 1–4 of `Scripts/` in order (Project A — flow + plasma).
4. For Project B, point `in.path` in `TB_CITEseq/Scripts/1.dsb_norm.R` at the Cell Ranger `multi` output directory and run scripts 1–7 in order.
5. The pipeline caches intermediate objects to `saved_R_dat/` (Project A) and `TB_CITEseq/saved_R_data/` (Project B) so later steps pick up where earlier ones left off.
