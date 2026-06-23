# BCG-Stim Infants 12 vs 44 Weeks — Parse Bio scRNA-seq

Single-cell RNA-seq (gene expression only) on **Parse Biosciences** of peripheral blood from infants sampled longitudinally and stimulated *in vitro* with **BCG** or media. The design enables three orthogonal comparisons in the same dataset: vaccine response (BCG vs media), maturation effect (44 weeks vs 12 weeks), and the impact of TB exposure status (IGRA⁺ vs IGRA⁻) on each.

This is the GEX companion to the flow / plasma / CITE-seq analysis in the [Lesley TB Mother–Infant Paired BCG Stim](../Lesley_TB_Mother_Infant_Paired_BCG_Stim) repository — same infants, same stim conditions, same timepoints, single-cell transcriptional readout.

---

## Cohort and design

| | |
|---|---|
| **Subjects** | 10 infants (`6044931`, `6044981`, `6061161`, `6062151`, `6063411`, `6063421`, `6063481`, `6065761`, `6067021`, `6071261`) |
| **Timepoints** | 12 weeks, 44 weeks (longitudinally paired within infant) |
| **In vitro stimulations** | BCG (vaccine), media (unstimulated control) |
| **Stratification** | Infant IGRA⁺ vs IGRA⁻ |
| **Total libraries** | 40 (10 × 2 timepoints × 2 stims, fully balanced) |
| **Platform** | Parse Biosciences combinatorial-barcoding scRNA-seq |
| **Modality** | Gene expression only (no ADT, no VDJ) |

---

## Repository structure

```
BCG_Stim_Infants_12_44_wks_PARSE_GEX/
├── Scripts/
│   ├── 1A.Create_Seurat_Parse_raw.R            # Read raw Parse outputs → merged Seurat
│   ├── 1B.Create_Seurat_Parse_filtered.R       # Read Parse filtered outputs (alt path)
│   ├── 2.Basic_QC.R                            # nFeature/nCount/mito/ribo/heme/platelet/complexity
│   ├── 3.Cell_Cycle_Doublet_Removal.R          # CellCycleScoring + scDblFinder
│   ├── 4.Integration.R                         # CCA + Harmony + MNN, side-by-side
│   └── 5.VLN_DGE.R                             # Three-way DGE design + violins
│
├── QC_Plots/
│   ├── UNFILTERED_QC/                          # Pre-Parse-filter exploration
│   │   ├── 1_Exploration/                      # Cells per orig.ident / sample_id / stim / timepoint
│   │   ├── 2_FeatureScatter_QC/                # Gene-vs-Mito, UMI-vs-Gene
│   │   └── 3_Histogram_Plots/                  # nGenes, UMI, Mito/Ribo/Heme/Platelet ratios, complexity
│   ├── FILTERED_QC/                            # Same panel post-Parse cell calling
│   ├── FINAL_QC/                               # Post-pipeline filtering: cells per stim/timepoint/orig.ident
│   ├── Cell_Cycle/                             # Phase scoring, ridge plots, PCA & UMAP by phase
│   └── Doublets/                               # 40 per-sample scDblFinder diagnostic plots
│
├── Integration_Plots/
│   ├── CCA_Clustering_Grid.png                 # Three integration methods compared side-by-side
│   ├── Harmony_Clustering_Grid.png
│   └── MNN_Clustering_Grid.png
│
├── DGE/
│   ├── Stim_vs_Unstim_byTimepoint/             # BCG vs media, per cluster, per timepoint
│   │   └── DGE_cluster<N>_week{12,44}_BCG_vs_media.csv
│   ├── Week44_vs_Week12_byStim/                # Maturation effect, per cluster, per stim
│   └── IGRApos_vs_IGRAneg_byStimTimepoint/     # TB-exposure effect, per cluster, per stim×TP
│       ├── DGE_cluster<N>__IGRApos_{BCG,Med}_{12,44}wk_vs_IGRAneg_{BCG,Med}_{12,44}wk.csv
│       └── IGRApos_vs_IGRAneg__ALL_clusters__ALL_results_LONG.csv
│
└── Heatmaps/
    └── IGRA_Stim_Timepoint_AvgExpr/
        ├── MarkerSelection/                    # Top 100 genes per lineage + source-cluster mapping
        │   ├── HeatGenes__CD4__n100.csv
        │   ├── HeatGenes__CD4__n100__gene_to_source_cluster.csv
        │   ├── HeatGenes__CD8_GD__n100.csv
        │   ├── HeatGenes__CD8_GD__n100__gene_to_source_cluster.csv
        │   ├── HeatGenes__Monocytes__n100.csv
        │   └── HeatGenes__Monocytes__n100__gene_to_source_cluster.csv
        └── Plots_SE/                           # Heatmaps + dot plots per lineage with SE-style labeling
            ├── Heatmap_SE__{CD4,CD8_GD,Monocytes}.png
            └── DotPlot2_SE__{CD4,CD8_GD,Monocytes}.png
```

---

## Pipeline overview

The analysis is implemented in `Scripts/`, run in numerical order.

### 1 – Object construction (`1A`, `1B`)
Per-sample Parse Bio output folders are loaded with `ReadParseBio()` and merged into a single Seurat object. `1A` reads from the `DGE_unfiltered/` outputs (raw Parse counts, before Parse's own cell-calling), and `1B` reads from `DGE_filtered/`. Per-cell metadata from Parse is preserved; `orig.ident`, `sample_id`, `stim`, `timepoint`, and IGRA status are propagated for downstream stratification.

### 2 – Basic QC (`2.Basic_QC.R`)
QC over `nFeature_RNA`, `nCount_RNA`, percent mitochondrial / ribosomal / hemoglobin / platelet, and the `log10GenesPerUMI` complexity score. Three views are kept side-by-side so the effect of Parse's own filter is auditable:

- `QC_Plots/UNFILTERED_QC/` — pre-Parse-filter
- `QC_Plots/FILTERED_QC/` — post-Parse-filter
- `QC_Plots/FINAL_QC/` — post-pipeline custom filter (cells-per-stim, cells-per-timepoint, etc.)

### 3 – Cell cycle and doublets (`3.Cell_Cycle_Doublet_Removal.R`)
Seurat `CellCycleScoring` with the standard S/G2M gene panels, plus per-sample doublet detection via **scDblFinder** (one diagnostic plot per library in `QC_Plots/Doublets/`). The cleaned, doublet-removed object is saved for downstream steps.

### 4 – Integration (`4.Integration.R`)
The cleaned object is split by sample and integrated three different ways with `IntegrateLayers()`:

- **CCA** (Canonical Correlation Analysis)
- **Harmony**
- **MNN** (fast Mutual Nearest Neighbors)

Side-by-side clustering grids in `Integration_Plots/` allow the integration choice to be made transparently.

### 5 – DGE & visualization (`5.VLN_DGE.R`)
Three differential expression designs are run **per cluster**, exporting one CSV per cluster × condition combination plus a long-format consolidated table for the IGRA contrast:

| Design | Question | Output folder |
|---|---|---|
| **BCG vs media** at each timepoint | What does the vaccine do to each cell type at 12 weeks? At 44 weeks? | `DGE/Stim_vs_Unstim_byTimepoint/` |
| **Week 44 vs Week 12** at each stim | How does each cell type mature across the 12-to-44 week window, with and without BCG? | `DGE/Week44_vs_Week12_byStim/` |
| **IGRA⁺ vs IGRA⁻** at each stim × timepoint | Does TB-exposure status alter cell-type-specific transcriptional state, in baseline (Med) or in BCG response, at each timepoint? | `DGE/IGRApos_vs_IGRAneg_byStimTimepoint/` |

**Average-expression heatmaps** are produced per lineage (CD4, CD8 / γδ, Monocytes), using the top 100 markers per lineage selected from cluster-level differential expression. Each heatmap displays mean expression across the **IGRA × Stim × Timepoint** combinations, providing a compact view of how each compartment's signature shifts across the design. Source-cluster traceability is preserved via `*_gene_to_source_cluster.csv`. Companion dot plots (`DotPlot2_SE__*.png`) accompany each heatmap.

---

## Dependencies

- [`Seurat`](https://satijalab.org/seurat/) v5 — multimodal single-cell framework
- [`scDblFinder`](https://bioconductor.org/packages/scDblFinder/) — doublet detection
- [`harmony`](https://github.com/immunogenomics/harmony) — batch correction
- [`SeuratWrappers`](https://github.com/satijalab/seurat-wrappers) — fast MNN integration
- [`scCustomize`](https://samuel-marsh.github.io/scCustomize/) — single-cell visualization
- [`Azimuth`](https://azimuth.hubmapconsortium.org/) — reference-based annotation utilities
- [`ComplexHeatmap`](https://bioconductor.org/packages/ComplexHeatmap/), `circlize` — heatmaps
- [`ggplot2`](https://ggplot2.tidyverse.org/), `patchwork`, `cowplot`, `ggrepel`
- `tidyverse`, `Matrix`, `data.table`

Parse Biosciences `ReadParseBio()` is provided by recent Seurat versions for reading Parse's combinatorial-barcoding output.

---

## Key questions the design addresses

| Question | Where to look |
|---|---|
| Per-cell-type response to BCG at 12 weeks | `DGE/Stim_vs_Unstim_byTimepoint/DGE_cluster*_week12_BCG_vs_media.csv` |
| How does the BCG response change between 12 and 44 weeks? | Compare 12-week vs 44-week files in `Stim_vs_Unstim_byTimepoint/`, plus `Week44_vs_Week12_byStim/DGE_cluster*_BCG.csv` |
| Baseline maturation of each cell type from 12 → 44 weeks | `DGE/Week44_vs_Week12_byStim/DGE_cluster*_media.csv` |
| Effect of TB exposure on baseline transcriptome | `DGE/IGRApos_vs_IGRAneg_byStimTimepoint/DGE_cluster*__IGRApos_Med_*wk_vs_IGRAneg_Med_*wk.csv` |
| Effect of TB exposure on the BCG response | `DGE/IGRApos_vs_IGRAneg_byStimTimepoint/DGE_cluster*__IGRApos_BCG_*wk_vs_IGRAneg_BCG_*wk.csv` |
| Combined view across CD4 / CD8&γδ / Monocyte lineages | `Heatmaps/IGRA_Stim_Timepoint_AvgExpr/Plots_SE/` |
| Pan-cluster summary of every IGRA contrast | `DGE/IGRApos_vs_IGRAneg_byStimTimepoint/IGRApos_vs_IGRAneg__ALL_clusters__ALL_results_LONG.csv` |

---

## Reproducing the analysis

1. Clone the repo and open `BCG_Stim_Infants_12_44_wks_PARSE_GEX.Rproj` in RStudio.
2. Update path variables near the top of each script — the scripts currently reference absolute paths under `/home/akshay-iyer/Documents/...`.
3. Place the per-sample Parse output folders (`sample_<id>/DGE_unfiltered/` and/or `DGE_filtered/`) under the directory referenced in `1A`/`1B`.
4. Run scripts 1–5 in order. Intermediate Seurat objects are saved between steps so later scripts pick up where earlier ones left off.
