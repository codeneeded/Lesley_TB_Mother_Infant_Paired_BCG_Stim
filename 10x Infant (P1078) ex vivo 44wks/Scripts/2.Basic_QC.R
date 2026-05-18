###############################################################################
#  TB CITE-seq Pipeline — Script 2: QC Visualizations + Filtering
#  Adapted from OPIS CITE-seq pipeline
###############################################################################

library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(data.table)
library(ggplot2)
library(biomaRt)
library(scCustomize)
library(patchwork)
library(qs2)

###############################################################################
#                         PATHS & GLOBAL VARIABLES
###############################################################################

project_root <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/TB_CITEseq"
qc_root      <- file.path(project_root, "QC")
out.path     <- file.path(project_root, "saved_R_data")

dir.create(qc_root, recursive = TRUE, showWarnings = FALSE)
setwd(qc_root)

# Load merged Seurat object
merged_seurat <- qs_read(file.path(out.path, "TB_Seuratv5_CITEseq_dsbnorm_merged.qs2"))

###############################################################################
#                         GLOBAL PLOT THEME + HELPERS
###############################################################################

# Clean, publication-ready base theme
theme_qc <- theme_minimal(base_size = 16) +
  theme(
    text             = element_text(family = "sans"),
    plot.title       = element_text(size = 20, face = "bold", hjust = 0.5,
                                    margin = margin(b = 12)),
    plot.subtitle    = element_text(size = 13, color = "grey40", hjust = 0.5,
                                    margin = margin(b = 10)),
    axis.title       = element_text(size = 14, face = "bold"),
    axis.text        = element_text(size = 12, color = "grey20"),
    axis.text.x      = element_text(angle = 40, hjust = 1, vjust = 1),
    legend.title     = element_text(size = 12, face = "bold"),
    legend.text      = element_text(size = 11),
    legend.position  = "right",
    panel.grid.major = element_line(color = "grey92", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(size = 13, face = "bold"),
    plot.margin      = margin(15, 15, 15, 15)
  )

# Color palettes
igra_pal <- c("IGRA-" = "#4393C3", "IGRA+" = "#D6604D")

group_pal <- c(
  "Infant:IGRA-_Mother:IGRA-" = "#4393C3",
  "Infant:IGRA-_Mother:IGRA+" = "#92C5DE",
  "Infant:IGRA+_Mother:IGRA-" = "#F4A582",
  "Infant:IGRA+_Mother:IGRA+" = "#D6604D"
)

# --- Reusable bar plot function ---
bar_plot <- function(data, x_var, fill_var, title, subtitle = NULL,
                     palette = NULL, show_count = TRUE) {
  p <- ggplot(data, aes(x = .data[[x_var]], fill = .data[[fill_var]])) +
    geom_bar(color = "white", linewidth = 0.3, width = 0.75) +
    theme_qc +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Number of Cells") +
    guides(fill = "none")
  
  if (!is.null(palette)) p <- p + scale_fill_manual(values = palette)
  
  if (show_count) {
    p <- p + geom_text(stat = "count", aes(label = after_stat(count)),
                       vjust = -0.5, size = 4.5, fontface = "bold", color = "grey30")
  }
  p
}

# --- Reusable density plot function ---
density_plot <- function(data, x_var, color_var, title, xlab_text,
                         vline_x = NULL, log_x = TRUE) {
  p <- ggplot(data, aes(x = .data[[x_var]], color = .data[[color_var]],
                        fill = .data[[color_var]])) +
    geom_density(alpha = 0.15, linewidth = 0.7) +
    theme_qc +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    labs(title = title, x = xlab_text, y = "Cell Density",
         color = "Sample", fill = "Sample")
  
  if (log_x) p <- p + scale_x_log10(labels = label_comma())
  
  if (!is.null(vline_x)) {
    p <- p + geom_vline(xintercept = vline_x, linetype = "dashed",
                        color = "grey30", linewidth = 0.6)
  }
  p
}

###############################################################################
#                           PRE-QC VISUALIZATIONS
###############################################################################

preqc_dir <- file.path(qc_root, "Pre-QC")
dir.create(preqc_dir, recursive = TRUE, showWarnings = FALSE)
setwd(preqc_dir)

metadata <- merged_seurat@meta.data

# --- Cells per Sample ---
png(file = "Cells_per_sample.png", width = 1800, height = 1200, res = 150)
bar_plot(metadata, "orig.ident", "orig.ident",
         title    = "Cells per Sample",
         subtitle = "Pre-QC")
dev.off()

# --- Cells per Infant IGRA Status ---
png(file = "Cells_per_IGRA_infant.png", width = 1200, height = 1000, res = 150)
bar_plot(metadata, "IGRA_infant", "IGRA_infant",
         title    = "Cells per Infant IGRA Status",
         subtitle = "Pre-QC",
         palette  = igra_pal) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
dev.off()

# --- Cells per Maternal IGRA Status ---
png(file = "Cells_per_IGRA_mother.png", width = 1200, height = 1000, res = 150)
bar_plot(metadata, "IGRA_mother", "IGRA_mother",
         title    = "Cells per Maternal IGRA Status",
         subtitle = "Pre-QC",
         palette  = igra_pal) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
dev.off()

# --- Cells per Combined IGRA Group ---
png(file = "Cells_per_IGRA_group.png", width = 1600, height = 1100, res = 150)
bar_plot(metadata, "IGRA_group", "IGRA_group",
         title    = "Cells per IGRA Group",
         subtitle = "Infant \u00d7 Maternal IGRA Status \u2014 Pre-QC",
         palette  = group_pal)
dev.off()

# --- QC Features Violin Plots ---
feats.1 <- c("nCount_RNA", "nFeature_RNA", "nCount_ADT", "nFeature_ADT",
             "percent_mito", "percent_ribo", "percent_hb", "percent_plat")

png(file = "Pre-QC_features_grouped.png", width = 2000, height = 1400, res = 150)
VlnPlot(merged_seurat, group.by = "orig.ident", features = feats.1,
        pt.size = 0.1, ncol = 2) +
  NoLegend()
dev.off()

# --- Density Plots ---
png(file = "UMI_Count.png", width = 1600, height = 1000, res = 150)
density_plot(metadata, "nCount_RNA", "orig.ident",
             title = "UMI Count Distribution", xlab_text = "nCount_RNA (log10)",
             vline_x = 500)
dev.off()

png(file = "nGenes.png", width = 1600, height = 1000, res = 150)
density_plot(metadata, "nFeature_RNA", "orig.ident",
             title = "Gene Count Distribution", xlab_text = "nFeature_RNA (log10)",
             vline_x = 5000)
dev.off()

png(file = "Complexity_Score.png", width = 1600, height = 1000, res = 150)
density_plot(metadata, "log10GenesPerUMI", "orig.ident",
             title = "Complexity (Genes per UMI)", xlab_text = "log10(Genes/UMI)",
             vline_x = 0.8, log_x = FALSE)
dev.off()

png(file = "Mito_Ratio.png", width = 1600, height = 1000, res = 150)
density_plot(metadata, "percent_mito", "orig.ident",
             title = "Mitochondrial %", xlab_text = "% Mitochondrial (log10)",
             vline_x = 15)
dev.off()

png(file = "Ribo_Ratio.png", width = 1600, height = 1000, res = 150)
density_plot(metadata, "percent_ribo", "orig.ident",
             title = "Ribosomal %", xlab_text = "% Ribosomal (log10)",
             vline_x = 5)
dev.off()

png(file = "Heme_Ratio.png", width = 1600, height = 1000, res = 150)
density_plot(metadata, "percent_hb", "orig.ident",
             title = "Hemoglobin %", xlab_text = "% Hemoglobin (log10)",
             vline_x = 20)
dev.off()

png(file = "Platelet_Ratio.png", width = 1600, height = 1000, res = 150)
density_plot(metadata, "percent_plat", "orig.ident",
             title = "Platelet %", xlab_text = "% Platelet (log10)",
             vline_x = 2)
dev.off()

# --- scCustomize QC Plots ---
p1 <- QC_Plots_Genes(seurat_object = merged_seurat, low_cutoff = 600, high_cutoff = 5500)
p2 <- QC_Plots_UMIs(seurat_object = merged_seurat, low_cutoff = 1200, high_cutoff = 45000)
p3 <- QC_Plots_Mito(seurat_object = merged_seurat, high_cutoff = 20)
p4 <- QC_Plots_Complexity(seurat_object = merged_seurat, high_cutoff = 0.8)

png(file = "Grouped_Cutoff.png", width = 2400, height = 1200, res = 150)
wrap_plots(p1, p2, p3, p4, ncol = 4)
dev.off()

# --- Scatter QC Plots ---
png(file = "UMIvsGene.png", width = 1800, height = 1200, res = 150)
QC_Plot_UMIvsGene(
  seurat_object    = merged_seurat,
  low_cutoff_gene  = 600,
  high_cutoff_gene = 5500,
  low_cutoff_UMI   = 500,
  high_cutoff_UMI  = 50000,
  group.by = "orig.ident"
)
dev.off()

png(file = "MitovsGene.png", width = 1800, height = 1200, res = 150)
QC_Plot_GenevsFeature(
  seurat_object       = merged_seurat,
  feature1            = "percent_mito",
  low_cutoff_gene     = 600,
  high_cutoff_gene    = 5500,
  high_cutoff_feature = 20,
  group.by = "orig.ident"
)
dev.off()

png(file = "MitovsGene_gradient.png", width = 1800, height = 1200, res = 150)
QC_Plot_UMIvsGene(
  seurat_object      = merged_seurat,
  meta_gradient_name = "percent_mito",
  low_cutoff_gene    = 600,
  high_cutoff_gene   = 5500,
  high_cutoff_UMI    = 45000
)
dev.off()

###############################################################################
#                              FILTERING
#  ** REVIEW Pre-QC plots before running this section **
#  Thresholds below are carried from OPIS — adjust based on your TB QC plots.
###############################################################################

merged_seurat <- JoinLayers(merged_seurat)
DefaultAssay(merged_seurat) <- "RNA"

filtered_seurat <- subset(
  x = merged_seurat,
  subset = (nCount_RNA       >= 500)  &
    (nFeature_RNA     >= 500)  &
    (log10GenesPerUMI >  0.80) &
    (percent_mito     <  15)   &
    (percent_ribo     >  5)    &
    (percent_hb       <  20)   &
    (percent_plat     <  2)
)

# Keep genes expressed in at least 10 cells
genes_in_10_cells <- rowSums(filtered_seurat@assays$RNA@layers$counts > 0) >= 10
filtered_seurat   <- subset(filtered_seurat, features = names(genes_in_10_cells[genes_in_10_cells]))

###############################################################################
#                        POST-QC VISUALIZATIONS
###############################################################################

postqc_dir <- file.path(qc_root, "Post-QC")
dir.create(postqc_dir, recursive = TRUE, showWarnings = FALSE)
setwd(postqc_dir)

post_meta <- filtered_seurat@meta.data

# --- Cells per Sample Post-QC ---
png(file = "Post-QC_Cells_per_sample.png", width = 1800, height = 1200, res = 150)
bar_plot(post_meta, "orig.ident", "orig.ident",
         title    = "Cells per Sample",
         subtitle = "Post-QC")
dev.off()

# --- Cells per Infant IGRA Status Post-QC ---
png(file = "Post-QC_Cells_per_IGRA_infant.png", width = 1200, height = 1000, res = 150)
bar_plot(post_meta, "IGRA_infant", "IGRA_infant",
         title    = "Cells per Infant IGRA Status",
         subtitle = "Post-QC",
         palette  = igra_pal) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
dev.off()

# --- Cells per Maternal IGRA Status Post-QC ---
png(file = "Post-QC_Cells_per_IGRA_mother.png", width = 1200, height = 1000, res = 150)
bar_plot(post_meta, "IGRA_mother", "IGRA_mother",
         title    = "Cells per Maternal IGRA Status",
         subtitle = "Post-QC",
         palette  = igra_pal) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
dev.off()

# --- Cells per IGRA Group Post-QC ---
png(file = "Post-QC_Cells_per_IGRA_group.png", width = 1600, height = 1100, res = 150)
bar_plot(post_meta, "IGRA_group", "IGRA_group",
         title    = "Cells per IGRA Group",
         subtitle = "Infant \u00d7 Maternal IGRA Status \u2014 Post-QC",
         palette  = group_pal)
dev.off()

# --- Post-QC Feature Violins ---
png(file = "Post-QC_features_grouped.png", width = 2000, height = 1400, res = 150)
VlnPlot(filtered_seurat, group.by = "orig.ident", features = feats.1,
        pt.size = 0.1, ncol = 2) +
  NoLegend()
dev.off()

###############################################################################
#                              SAVE
###############################################################################

qs_save(filtered_seurat,
        file = file.path(out.path, "TB_Seuratv5_filtered_seurat.qs2"))

message("Done! Filtered Seurat object saved.")
message("Pre-QC cells:  ", ncol(merged_seurat))
message("Post-QC cells: ", ncol(filtered_seurat))

