# Quality Control Visualizations
##Ref -> https://www.singlecellcourse.org/scrna-seq-analysis-with-bioconductor.html
##Ref -> http://bioconductor.org/books/3.15/OSCA.intro/getting-scrna-seq-datasets.html
#Ref -> https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html
# Ref -> https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html

##############################################################
# QC COMMENTS (Applicable to Both RAW and PARSEF Pipelines)
#
# The following QC steps are applied to all Seurat objects 
# (RAW and PARSEF). Each metric is assessed against expected 
# biological and technical thresholds:
#
# 1. nGenes (nFeature_RNA):
#    - High-quality data should have a single main peak.
#    - Shoulders or bimodal shapes can indicate failed cells 
#      or heterogeneous biology (e.g. quiescent populations).
#    - Cutoff example: > 250 features retained.
#
# 2. UMI Counts (nCount_RNA):
#    - Low UMI counts (< 500) indicate low sequencing depth.
#    - Between 500–1000 is borderline but may still be usable.
#
# 3. Complexity Score (log10GenesPerUMI):
#    - Ratio of genes detected per UMI (novelty).
#    - Low complexity (< 0.8) indicates over-sequencing a 
#      narrow set of transcripts (e.g. poor cell diversity).
#
# 4. Percent Mitochondrial Reads (percent_mito):
#    - High % mitochondrial reads (> 15–20%) can indicate 
#      stressed/dying cells.
#    - May vary by tissue type; context-dependent.
#
# 5. Percent Ribosomal Reads (percent_ribo):
#    - Very low ribosomal RNA (< 5%) can indicate poor quality.
#    - Very high ribosomal RNA may indicate unremoved debris 
#      or technical bias.
# Ribosomal content (percent_ribo) is expected to be low (~1%) in Parse data due to rRNA depletion.
# No hard cutoff applied. Cells with percent_ribo > 5% will be reviewed manually,
# but filtering will only occur if high-ribo clusters are identified as artifacts.
# 10x Genomics scRNA-seq ribosomal content:
# Typical cells show ~5–15% ribosomal gene expression.
# High percent_ribo (>20%) may indicate stressed, damaged, or highly proliferative cells.
# Unlike Parse (where ~1% is normal), 10x baseline ribo content is higher.
# 6. Percent Hemoglobin Genes (percent_hb):
#    - High levels (> 20%) indicate strong RBC contamination.
#
# 7. Percent Platelet Genes (percent_plat):
#    - High levels (> 2%) indicate platelet contamination.
#
# All cutoffs should be adjusted based on dataset-specific
# distributions, protocol differences, and biological context.
##############################################################
##############################################################
# BCG Parse: Quality Control Visualizations
##############################################################

# Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(scCustomize)
library(qs2) ### For more efficient save/loading
library(scales)

# Base directory for BCG Parse project
base_dir <- "/home/akshay-iyer/Documents/BCG_Stim_Infants_12_44_wks_PARSE_GEX"

qc_dir_raw <- file.path(base_dir, "QC_Plots", "UNFILTERED_QC")
dir.create(qc_dir_raw, recursive = TRUE, showWarnings = FALSE)

# Load UNFILTERED Seurat object
seu_raw <- readRDS(file.path(base_dir, "saved_R_data", "BCG_PARSE_unfiltered_merged.rds"))
Idents(seu_raw) <- "orig.ident"

md_raw <- seu_raw@meta.data

# Create combined factor for nicer colors
md_raw$tp_stim <- interaction(md_raw$timepoint, md_raw$stim, drop = TRUE)

# QC axis helper functions
pretty_log10_counts <- scale_x_log10(
  breaks = c(100, 250, 500, 1000, 3000, 6000, 10000, 40000),
  labels = label_number(accuracy = 1, big.mark = ",")
)

pretty_log10_percent <- scale_x_log10(
  breaks = c(0.1, 0.5, 1, 2, 5, 10, 20, 50),
  labels = label_number(accuracy = 0.1)
)

qc_theme <- function() {
  theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

##############################################################
# 1. EXPLORATION PLOTS
##############################################################

exploratory_dir_raw <- file.path(qc_dir_raw, "1_Exploration")
dir.create(exploratory_dir_raw, showWarnings = FALSE)

# Cells per orig.ident
png(file.path(exploratory_dir_raw, "Cells_per_orig.ident.png"), width = 2200, height = 1400)
md_raw %>%
  ggplot(aes(x = orig.ident, fill = timepoint)) +
  geom_bar() +
  qc_theme() +
  labs(
    title = "Number of Cells per Library (orig.ident)",
    x = "Library (sample_id_timepoint_stim)", y = "Cell Count"
  )
dev.off()

# Cells per sample_id
png(file.path(exploratory_dir_raw, "Cells_per_sample_id.png"), width = 1600, height = 1200)
md_raw %>%
  ggplot(aes(x = sample_id, fill = sample_id)) +
  geom_bar() +
  qc_theme() +
  labs(title = "Number of Cells per Sample ID", x = "Sample ID", y = "Cell Count")
dev.off()

# Cells per timepoint
png(file.path(exploratory_dir_raw, "Cells_per_timepoint.png"), width = 1400, height = 1000)
md_raw %>%
  ggplot(aes(x = timepoint, fill = timepoint)) +
  geom_bar() +
  qc_theme() +
  labs(title = "Number of Cells per Timepoint", x = "Timepoint", y = "Cell Count")
dev.off()

# Cells per stim
png(file.path(exploratory_dir_raw, "Cells_per_stim.png"), width = 1400, height = 1000)
md_raw %>%
  ggplot(aes(x = stim, fill = stim)) +
  geom_bar() +
  qc_theme() +
  labs(title = "Number of Cells per Stimulation", x = "Stimulation", y = "Cell Count")
dev.off()

##############################################################
# 2. FEATURE SCATTER QC PLOTS
##############################################################

fs_dir_raw <- file.path(qc_dir_raw, "2_FeatureScatter_QC")
dir.create(fs_dir_raw, showWarnings = FALSE)

# UMI vs Gene coloured by mito
png(file.path(fs_dir_raw, "UMI_vs_Gene_by_Mito.png"), width = 1800, height = 1200)
QC_Plot_UMIvsGene(
  seurat_object = seu_raw,
  meta_gradient_name = "percent_mito",
  low_cutoff_gene = 200,
  high_cutoff_gene = 6000,
  low_cutoff_UMI = 500,
  high_cutoff_UMI = 40000,
  meta_gradient_low_cutoff = 20,
  combination = TRUE
)
dev.off()

# Gene vs percent_mito
png(file.path(fs_dir_raw, "Gene_vs_Mito.png"), width = 1800, height = 1200)
QC_Plot_GenevsFeature(
  seurat_object = seu_raw,
  feature1 = "percent_mito",
  low_cutoff_gene = 200,
  high_cutoff_gene = 6000,
  high_cutoff_feature = 20
)
dev.off()

##############################################################
# 3. HISTOGRAM / DENSITY QC PLOTS (with clean axes)
##############################################################

hist_dir_raw <- file.path(qc_dir_raw, "3_Histogram_Plots")
dir.create(hist_dir_raw, showWarnings = FALSE)

md_raw$log10GenesPerUMI <- log10(md_raw$nFeature_RNA / md_raw$nCount_RNA)

### UMI COUNT DISTRIBUTION
png(file.path(hist_dir_raw, "UMI_Count.png"), width = 1800, height = 1200)
md_raw %>%
  ggplot(aes(x = nCount_RNA, color = tp_stim, fill = tp_stim)) +
  geom_density(alpha = 0.25) +
  pretty_log10_counts +
  qc_theme() +
  geom_vline(xintercept = 500, linetype = "dashed") +
  labs(
    title = "UMI Count Distribution (Unfiltered)",
    x = "UMI Count per Cell", color = "Timepoint_Stimulation", fill = "Timepoint_Stimulation"
  )
dev.off()

### nGENES DISTRIBUTION
png(file.path(hist_dir_raw, "nGenes.png"), width = 1800, height = 1200)
md_raw %>%
  ggplot(aes(x = nFeature_RNA, color = tp_stim, fill = tp_stim)) +
  geom_density(alpha = 0.25) +
  pretty_log10_counts +
  qc_theme() +
  geom_vline(xintercept = 250, linetype = "dashed") +
  labs(
    title = "Detected Gene Count Distribution (Unfiltered)",
    x = "Detected Genes per Cell", color = "Timepoint_Stimulation", fill = "Timepoint_Stimulation"
  )
dev.off()

### COMPLEXITY SCORE
png(file.path(hist_dir_raw, "Complexity_Score.png"), width = 1800, height = 1200)
md_raw %>%
  ggplot(aes(x = log10GenesPerUMI, color = tp_stim, fill = tp_stim)) +
  geom_density(alpha = 0.25) +
  qc_theme() +
  geom_vline(xintercept = 0.8, linetype = "dashed") +
  labs(
    title = "Complexity Score (log10 Genes / UMIs)",
    x = "log10(Genes per UMI)"
  )
dev.off()

### MITOCHONDRIAL RATIO
png(file.path(hist_dir_raw, "Mito_Ratio.png"), width = 1800, height = 1200)
md_raw %>%
  ggplot(aes(x = percent_mito, color = tp_stim, fill = tp_stim)) +
  geom_density(alpha = 0.25) +
  pretty_log10_percent +
  qc_theme() +
  geom_vline(xintercept = 15, linetype = "dashed") +
  labs(
    title = "Mitochondrial Percentage Distribution (Unfiltered)",
    x = "percent_mito"
  )
dev.off()

### RIBOSOMAL RATIO
png(file.path(hist_dir_raw, "Ribo_Ratio.png"), width = 1800, height = 1200)
md_raw %>%
  ggplot(aes(x = percent_ribo, color = tp_stim, fill = tp_stim)) +
  geom_density(alpha = 0.25) +
  pretty_log10_percent +
  qc_theme() +
  geom_vline(xintercept = 5, linetype = "dashed") +
  labs(
    title = "Ribosomal Percentage Distribution (Unfiltered)",
    x = "percent_ribo"
  )
dev.off()

### HEMOGLOBIN RATIO
png(file.path(hist_dir_raw, "Heme_Ratio.png"), width = 1800, height = 1200)
md_raw %>%
  ggplot(aes(x = percent_hb, color = tp_stim, fill = tp_stim)) +
  geom_density(alpha = 0.25) +
  pretty_log10_percent +
  qc_theme() +
  geom_vline(xintercept = 20, linetype = "dashed") +
  labs(
    title = "Hemoglobin Percentage Distribution (Unfiltered)",
    x = "percent_hb"
  )
dev.off()

### PLATELET RATIO
png(file.path(hist_dir_raw, "Platelet_Ratio.png"), width = 1800, height = 1200)
md_raw %>%
  ggplot(aes(x = percent_plat, color = tp_stim, fill = tp_stim)) +
  geom_density(alpha = 0.25) +
  pretty_log10_percent +
  qc_theme() +
  geom_vline(xintercept = 2, linetype = "dashed") +
  labs(
    title = "Platelet Marker Percentage Distribution (Unfiltered)",
    x = "percent_plat"
  )
dev.off()

##############################################################
# QC for FILTERED MERGED SEURAT OBJECT (pre-QC filtering)
##############################################################

qc_dir_filtered <- file.path(base_dir, "QC_Plots", "FILTERED_QC")
dir.create(qc_dir_filtered, recursive = TRUE, showWarnings = FALSE)

seu_filtered <- qs_read(file.path(base_dir, "saved_R_data", "BCG_PARSE_filtered_merged.qs2"))
Idents(seu_filtered) <- "orig.ident"

########################################
# 1. Exploration (FILTERED, pre-QC)
########################################

exploratory_dir_filtered <- file.path(qc_dir_filtered, "1_Exploration")
dir.create(exploratory_dir_filtered, showWarnings = FALSE)

md_filt <- seu_filtered@meta.data

# Cells per orig.ident
png(file.path(exploratory_dir_filtered, "Cells_per_orig.ident.png"), width = 2200, height = 1400)
md_filt %>%
  ggplot(aes(x = orig.ident, fill = timepoint)) +
  geom_bar() +
  qc_theme() +
  labs(title = "Number of Cells per Library (Filtered, Pre-QC)",
       x = "Library (sample_id_timepoint_stim)", y = "Cell Count")
dev.off()

# Cells per sample_id
png(file.path(exploratory_dir_filtered, "Cells_per_sample_id.png"), width = 1600, height = 1200)
md_filt %>%
  ggplot(aes(x = sample_id, fill = sample_id)) +
  geom_bar() +
  qc_theme() +
  labs(title = "Number of Cells per Sample ID (Filtered, Pre-QC)",
       x = "Sample ID", y = "Cell Count")
dev.off()

# Cells per timepoint
png(file.path(exploratory_dir_filtered, "Cells_per_timepoint.png"), width = 1200, height = 1000)
md_filt %>%
  ggplot(aes(x = timepoint, fill = timepoint)) +
  geom_bar() +
  qc_theme() +
  labs(title = "Number of Cells per Timepoint (Filtered, Pre-QC)",
       x = "Timepoint", y = "Cell Count")
dev.off()

# Cells per stim
png(file.path(exploratory_dir_filtered, "Cells_per_stim.png"), width = 1200, height = 1000)
md_filt %>%
  ggplot(aes(x = stim, fill = stim)) +
  geom_bar() +
  qc_theme() +
  labs(title = "Number of Cells per Stimulation (Filtered, Pre-QC)",
       x = "Stimulation", y = "Cell Count")
dev.off()

########################################
# 2. FeatureScatter QC (FILTERED, pre-QC)
########################################

fs_dir_filtered <- file.path(qc_dir_filtered, "2_FeatureScatter_QC")
dir.create(fs_dir_filtered, showWarnings = FALSE)

png(file.path(fs_dir_filtered, "UMI_vs_Gene_colored_by_Mito.png"), width = 1800, height = 1200)
QC_Plot_UMIvsGene(
  seurat_object = seu_filtered,
  meta_gradient_name = "percent_mito",
  low_cutoff_gene = 200,
  high_cutoff_gene = 6000,
  low_cutoff_UMI = 500,
  high_cutoff_UMI = 40000,
  meta_gradient_low_cutoff = 20,
  combination = TRUE
)
dev.off()

png(file.path(fs_dir_filtered, "Gene_vs_Mito.png"), width = 1800, height = 1200)
QC_Plot_GenevsFeature(
  seurat_object = seu_filtered,
  feature1 = "percent_mito",
  low_cutoff_gene = 200,
  high_cutoff_gene = 6000,
  high_cutoff_feature = 20
)
dev.off()

########################################
# 3. Histogram QC (FILTERED, pre-QC)
########################################

hist_dir_filtered <- file.path(qc_dir_filtered, "3_Histogram_Plots")
dir.create(hist_dir_filtered, showWarnings = FALSE)

md_filt$log10GenesPerUMI <- log10(md_filt$nFeature_RNA / md_filt$nCount_RNA)
md_filt$tp_stim <- interaction(md_filt$timepoint, md_filt$stim, drop = TRUE)

png(file = file.path(hist_dir_filtered, "UMI_Count.png"), width = 1800, height = 1200)
md_filt %>% 
  ggplot(aes(color = tp_stim, x = nCount_RNA, fill = tp_stim)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + qc_theme() +
  geom_vline(xintercept = 500) +
  labs(title = "UMI Count Distribution (Filtered, Pre-QC)",
       x = "nCount_RNA (log10)", color = "Timepoint_Stim", fill = "Timepoint_Stim")
dev.off()

png(file = file.path(hist_dir_filtered, "nGenes.png"), width = 1800, height = 1200)
md_filt %>% 
  ggplot(aes(color = tp_stim, x = nFeature_RNA, fill = tp_stim)) + 
  geom_density(alpha = 0.2) + scale_x_log10() + qc_theme() +
  geom_vline(xintercept = 250) +
  labs(title = "Gene Count Distribution (Filtered, Pre-QC)",
       x = "nFeature_RNA (log10)", color = "Timepoint_Stim", fill = "Timepoint_Stim")
dev.off()

png(file = file.path(hist_dir_filtered, "Complexity_Score.png"), width = 1800, height = 1200)
md_filt %>%
  ggplot(aes(x = log10GenesPerUMI, color = tp_stim, fill = tp_stim)) +
  geom_density(alpha = 0.2) + qc_theme() +
  geom_vline(xintercept = 0.8) +
  labs(title = "Complexity Score (Genes per UMI, Filtered, Pre-QC)",
       x = "log10(Genes per UMI)", color = "Timepoint_Stim", fill = "Timepoint_Stim")
dev.off()

png(file = file.path(hist_dir_filtered, "Mito_Ratio.png"), width = 1800, height = 1200)
md_filt %>% 
  ggplot(aes(color = tp_stim, x = percent_mito, fill = tp_stim)) + 
  geom_density(alpha = 0.2) + qc_theme() +
  geom_vline(xintercept = 15) +
  labs(title = "Mitochondrial Ratio Distribution (Filtered, Pre-QC)",
       x = "percent_mito", color = "Timepoint_Stim", fill = "Timepoint_Stim")
dev.off()

png(file = file.path(hist_dir_filtered, "Ribo_Ratio.png"), width = 1800, height = 1200)
md_filt %>% 
  ggplot(aes(color = tp_stim, x = percent_ribo, fill = tp_stim)) + 
  geom_density(alpha = 0.2) + qc_theme() +
  geom_vline(xintercept = 5) +
  labs(title = "Ribosomal Ratio Distribution (Filtered, Pre-QC)",
       x = "percent_ribo", color = "Timepoint_Stim", fill = "Timepoint_Stim")
dev.off()

png(file = file.path(hist_dir_filtered, "Heme_Ratio.png"), width = 1800, height = 1200)
md_filt %>% 
  ggplot(aes(color = tp_stim, x = percent_hb, fill = tp_stim)) + 
  geom_density(alpha = 0.2) + qc_theme() +
  geom_vline(xintercept = 20) +
  labs(title = "Hemoglobin Ratio Distribution (Filtered, Pre-QC)",
       x = "percent_hb", color = "Timepoint_Stim", fill = "Timepoint_Stim")
dev.off()

png(file = file.path(hist_dir_filtered, "Platelet_Ratio.png"), width = 1800, height = 1200)
md_filt %>% 
  ggplot(aes(color = tp_stim, x = percent_plat, fill = tp_stim)) + 
  geom_density(alpha = 0.2) + qc_theme() +
  geom_vline(xintercept = 2) +
  labs(title = "Platelet Ratio Distribution (Filtered, Pre-QC)",
       x = "percent_plat", color = "Timepoint_Stim", fill = "Timepoint_Stim")
dev.off()


##############################################################
# FINAL QC FILTERING + POST-QC VISUALIZATION
##############################################################

qc_final_dir <- file.path(base_dir, "QC_Plots", "FINAL_QC")
dir.create(qc_final_dir, recursive = TRUE, showWarnings = FALSE)

# Pre-QC violin plots by orig.ident
feats_parse <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_hb", "percent_plat")

png(file.path(qc_final_dir, "PreQC_Features_by_orig.ident.png"), width = 2200, height = 1400)
VlnPlot(seu_filtered, group.by = "orig.ident", features = feats_parse, pt.size = 0.1, ncol = 3) +
  NoLegend()
dev.off()

# Also violin plots grouped by stim and timepoint (cleaner)
png(file.path(qc_final_dir, "PreQC_Features_by_stim.png"), width = 1600, height = 1200)
VlnPlot(seu_filtered, group.by = "stim", features = feats_parse, pt.size = 0.1, ncol = 3) +
  NoLegend()
dev.off()

png(file.path(qc_final_dir, "PreQC_Features_by_timepoint.png"), width = 1600, height = 1200)
VlnPlot(seu_filtered, group.by = "timepoint", features = feats_parse, pt.size = 0.1, ncol = 3) +
  NoLegend()
dev.off()

# Parse-specific cutoffs
seu_filtered <- JoinLayers(seu_filtered)

filtered_final <- subset(
  seu_filtered,
  subset = (nCount_RNA >= 500) &
    (nFeature_RNA >= 600) &
    (percent_mito < 15) &
    (percent_hb < 10) &
    (percent_plat < 2)
)

# Post-QC: cells per orig.ident / timepoint / stim

md_post <- filtered_final@meta.data

png(file.path(qc_final_dir, "Cells_per_orig.ident_PostQC.png"), width = 2200, height = 1400)
md_post %>%
  ggplot(aes(x = orig.ident, fill = timepoint)) +
  geom_bar() +
  qc_theme() +
  labs(title = "Number of Cells per Library (Post-QC)",
       x = "Library (sample_id_timepoint_stim)", y = "Cell Count")
dev.off()

png(file.path(qc_final_dir, "Cells_per_timepoint_PostQC.png"), width = 1200, height = 1000)
md_post %>%
  ggplot(aes(x = timepoint, fill = timepoint)) +
  geom_bar() +
  qc_theme() +
  labs(title = "Number of Cells per Timepoint (Post-QC)",
       x = "Timepoint", y = "Cell Count")
dev.off()

png(file.path(qc_final_dir, "Cells_per_stim_PostQC.png"), width = 1200, height = 1000)
md_post %>%
  ggplot(aes(x = stim, fill = stim)) +
  geom_bar() +
  qc_theme() +
  labs(title = "Number of Cells per Stimulation (Post-QC)",
       x = "Stimulation", y = "Cell Count")
dev.off()

# Post-QC violin plots by stim/timepoint
png(file.path(qc_final_dir, "PostQC_Features_by_stim.png"), width = 1600, height = 1200)
VlnPlot(filtered_final, group.by = "stim", features = feats_parse, pt.size = 0.1, ncol = 3) +
  NoLegend()
dev.off()

png(file.path(qc_final_dir, "PostQC_Features_by_timepoint.png"), width = 1600, height = 1200)
VlnPlot(filtered_final, group.by = "timepoint", features = feats_parse, pt.size = 0.1, ncol = 3) +
  NoLegend()
dev.off()

# Save post-QC object
qs_save(filtered_final, file.path(base_dir, "saved_R_data", "BCG_PARSE_filtered_postQC.qs2"))
