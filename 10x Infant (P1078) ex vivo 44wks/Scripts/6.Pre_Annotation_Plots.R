###############################################################################
#  TB CITE-seq Pipeline — Script 6: Pre-Annotation Plots
#  Feature plots, violin plots, heatmaps, average expression export
###############################################################################

library(Seurat)
library(SeuratExtend)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(grid)
library(qs2)

###############################################################################
#                         PATHS & SETUP
###############################################################################

base_dir  <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/TB_CITEseq"
load.path <- file.path(base_dir, "saved_R_data")

# Main output folder
annot_dir <- file.path(base_dir, "Pre_Annotation_Plots")
dir.create(annot_dir, recursive = TRUE, showWarnings = FALSE)

# Sub-folders
out_feat_adt  <- file.path(annot_dir, "FeaturePlots_ADT")
out_feat_rna  <- file.path(annot_dir, "FeaturePlots_RNA")
out_vln_adt   <- file.path(annot_dir, "ViolinPlots_ADT")
out_vln_rna   <- file.path(annot_dir, "ViolinPlots_RNA")
out_heat      <- file.path(annot_dir, "Heatmaps")
out_avgexpr   <- file.path(annot_dir, "AvgExpression")

for (d in c(out_feat_adt, out_feat_rna, out_vln_adt, out_vln_rna, out_heat, out_avgexpr)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

###############################################################################
#                         LOAD OBJECT
###############################################################################

TB_ALL <- qs_read(file.path(load.path, "TB_ALL_WNN.qs2"))

# Set active cluster identity to res 0.6
Idents(TB_ALL) <- "wsnn_res.0.6"
table(TB_ALL$Sample_ID,TB_ALL$IGRA_mother)
# Ensure numeric ordering of cluster levels
lvls <- as.character(sort(as.numeric(levels(Idents(TB_ALL)))))
Idents(TB_ALL) <- factor(Idents(TB_ALL), levels = lvls)
TB_ALL$seurat_clusters <- Idents(TB_ALL)

message("Clusters: ", paste(levels(Idents(TB_ALL)), collapse = ", "))
message("Total cells: ", ncol(TB_ALL))

###############################################################################
#  MARKER LISTS
###############################################################################

# ── ADT: ALL proteins in assay ──
DefaultAssay(TB_ALL) <- "ADT"
all_adt_markers <- sort(rownames(TB_ALL))
isotype_genes   <- c("Mouse-IgG1", "Mouse-IgG2a", "Mouse-IgG2b",
                     "Rat-IgG2b", "Rat-IgG1", "Rat-IgG2a", "Hamster-IgG")
adt_available   <- setdiff(all_adt_markers, isotype_genes)

cat("ADT markers (excl. isotypes):", length(adt_available), "\n")

# ── RNA: Grouped marker lists ──
rna_markers <- list(
  # --- T cell core ---
  `CD8 Identity` = c(
    "CD8A", "CD8B", "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"
  ),
  `CD4 T cell` = c(
    "CD4", "IL7R", "FOXP3", "IL2RA", "CTLA4", "IKZF2",
    "CCR4", "RORC", "GATA3", "TBX21", "BCL6",
    "CXCR5", "MAF", "ICOS", "BATF"
  ),
  `Treg` = c(
    "FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "TIGIT"
  ),
  # --- Naive / Quiescence ---
  `Naive - Quiescence` = c(
    "CCR7", "SELL", "TCF7", "LEF1", "BACH2", "IL7R", "BCL2",
    "S1PR1", "KLF2", "FOXP1", "SATB1"
  ),
  # --- Effector / Cytotoxicity ---
  `Effector - Cytotoxicity` = c(
    "GZMB", "GZMA", "GZMH", "GZMK", "GZMM",
    "PRF1", "GNLY", "NKG7", "KLRG1", "CX3CR1",
    "FGFBP2", "FCGR3A"
  ),
  `Effector TFs` = c(
    "TBX21", "EOMES", "RUNX3", "PRDM1", "ZEB2", "ID2"
  ),
  # --- Exhaustion / Inhibitory ---
  `Exhaustion` = c(
    "TOX", "TOX2", "PDCD1", "HAVCR2", "TIGIT", "LAG3",
    "CTLA4", "ENTPD1", "CD160", "CD244"
  ),
  # --- Memory / Stemness ---
  `Memory - Stemness` = c(
    "IL7R", "CD27", "CD28", "BACH2", "TCF7", "BCL2",
    "BCL6", "CXCR5", "CCR7"
  ),
  # --- Activation / Proliferation ---
  `Activation - Proliferation` = c(
    "MKI67", "HLA-DRA", "HLA-DRB1", "CD38", "FAS",
    "ICOS", "TNFRSF9", "CXCR6", "CD69"
  ),
  # --- Tissue Residency ---
  `Tissue Residency` = c(
    "ITGAE", "ITGA1", "CD69", "CXCR6", "ZNF683", "PRDM1"
  ),
  # --- GammaDelta TCR ---
  `GammaDelta TCR` = c(
    "TRDV1", "TRDV2", "TRGV9", "TRDC", "TRGC1", "TRGC2"
  ),
  # --- MAIT ---
  `MAIT` = c(
    "TRAV1-2", "SLC4A10", "KLRB1", "ZBTB16", "RORC"
  ),
  # --- NK-like / Innate ---
  `NK-like - Innate` = c(
    "TYROBP", "KLRD1", "KIR2DL3", "KIR3DL1", "NCAM1",
    "KLRF1", "KLRC1", "KLRC2"
  ),
  # --- NK cell ---
  `NK cell` = c(
    "NCAM1", "KLRF1", "KLRC1", "KLRB1", "NCR1", "NCR3",
    "FCGR3A", "TYROBP", "FCER1G", "SH2D1B",
    "SPON2", "CLIC3", "MYOM2"
  ),
  # --- Type I IFN ---
  `Type I IFN` = c(
    "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI44L", "MX1",
    "OAS1", "STAT1", "IRF7"
  ),
  # --- Stress / Heat Shock ---
  `Stress - Heat Shock` = c(
    "HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1", "IFI27"
  ),
  # --- Chemokines / Cytokines ---
  `Chemokines - Cytokines` = c(
    "CCL3", "CCL4", "CCL4L2", "CCL5", "XCL1", "XCL2",
    "IFNG", "TNF", "IL2", "CSF2"
  ),
  # --- Terminal Differentiation ---
  `Terminal Differentiation` = c(
    "B3GAT1", "KLRG1", "CX3CR1", "ZEB2", "GZMB"
  ),
  # --- B cell ---
  `B cell` = c(
    "CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "BANK1",
    "BLK", "BLNK", "TCL1A", "IGHM", "IGHD", "IGKC", "IGLC2"
  ),
  # --- Plasma cell ---
  `Plasma cell` = c(
    "SDC1", "XBP1", "MZB1", "JCHAIN", "IGHG1", "IGHA1",
    "PRDM1", "IRF4"
  ),
  # --- Monocyte / Macrophage ---
  `Monocyte - Macrophage` = c(
    "CD14", "LYZ", "S100A8", "S100A9", "S100A12",
    "FCGR3A", "CSF1R", "CD68", "MARCO",
    "VCAN", "FCN1", "MNDA"
  ),
  # --- cDC ---
  `cDC` = c(
    "FCER1A", "CLEC10A", "CD1C", "CLEC9A", "XCR1",
    "BATF3", "IRF8", "IRF4", "LAMP3"
  ),
  # --- pDC ---
  `pDC` = c(
    "LILRA4", "CLEC4C", "IRF7", "TCF4", "GZMB",
    "IL3RA", "JCHAIN"
  ),
  # --- Platelet ---
  `Platelet` = c(
    "PPBP", "PF4", "GP9", "ITGA2B", "TUBB1", "TREML1"
  ),
  # --- Erythrocyte ---
  `Erythrocyte` = c(
    "HBB", "HBA1", "HBA2", "ALAS2", "SLC4A1", "GYPA"
  ),
  # --- HSC / Progenitor ---
  `HSC - Progenitor` = c(
    "CD34", "SPINK2", "AVP", "CRHBP", "CYTL1", "SOX4", "HOPX"
  ),
  # --- ILC ---
  `ILC` = c(
    "IL7R", "KIT", "GATA3", "RORC", "IL1RL1",
    "AREG", "IL22", "NCR2"
  )
)

# Flatten and check availability
rna_markers_flat <- unique(unlist(rna_markers))

DefaultAssay(TB_ALL) <- "RNA"
rna_available <- rna_markers_flat[rna_markers_flat %in% rownames(TB_ALL)]
rna_missing   <- setdiff(rna_markers_flat, rna_available)

cat("RNA markers: ", length(rna_available), " available, ",
    length(rna_missing), " missing\n")
if (length(rna_missing) > 0) cat("Missing:", paste(rna_missing, collapse = ", "), "\n")

###############################################################################
#  AVERAGE EXPRESSION EXPORT
###############################################################################

message("Exporting average expression tables...")

# Helper: min-max scale rows to 0-1
scale_01 <- function(mat) {
  scaled <- t(apply(mat, 1, function(x) {
    rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    if (rng == 0 || is.na(rng)) return(rep(0.5, length(x)))
    (x - min(x, na.rm = TRUE)) / rng
  }))
  colnames(scaled) <- colnames(mat)
  # Replace any remaining NaN/NA with 0
  scaled[is.na(scaled)] <- 0
  scaled
}

# Helper: sort columns numerically (strips 'g' prefix from AverageExpression output)
sort_cluster_cols <- function(mat) {
  cnames <- colnames(mat)
  nums <- as.numeric(gsub("^g", "", cnames))
  ord <- order(nums)
  mat[, ord, drop = FALSE]
}

# --- RNA ---
DefaultAssay(TB_ALL) <- "RNA"

avg_rna_raw <- AverageExpression(
  TB_ALL, assays = "RNA", features = rna_available,
  group.by = "seurat_clusters", slot = "data"
)$RNA

avg_rna_full <- AverageExpression(
  TB_ALL, assays = "RNA", group.by = "seurat_clusters", slot = "data"
)$RNA

avg_rna_scaled <- scale_01(as.matrix(avg_rna_raw))

write.csv(avg_rna_raw,
          file.path(out_avgexpr, "AvgExpr_RNA_markers_per_cluster.csv"), row.names = TRUE)
write.csv(round(avg_rna_scaled, 4),
          file.path(out_avgexpr, "AvgExpr_RNA_markers_per_cluster_scaled.csv"), row.names = TRUE)
write.csv(avg_rna_full,
          file.path(out_avgexpr, "AvgExpr_RNA_AllGenes_per_cluster.csv"), row.names = TRUE)

# --- ADT ---
DefaultAssay(TB_ALL) <- "ADT"

avg_adt_raw <- AverageExpression(
  TB_ALL, assays = "ADT", features = adt_available,
  group.by = "seurat_clusters", slot = "data"
)$ADT

avg_adt_scaled <- scale_01(as.matrix(log10(avg_adt_raw + 1)))

write.csv(avg_adt_raw,
          file.path(out_avgexpr, "AvgExpr_ADT_per_cluster.csv"), row.names = TRUE)
write.csv(round(avg_adt_scaled, 4),
          file.path(out_avgexpr, "AvgExpr_ADT_per_cluster_scaled.csv"), row.names = TRUE)

# --- Percent expression ---
# Join layers first (Seurat v5 requires this for GetAssayData across split layers)
TB_ALL[["RNA"]] <- JoinLayers(TB_ALL[["RNA"]])

DefaultAssay(TB_ALL) <- "RNA"
rna_data_mat <- GetAssayData(TB_ALL, slot = "data")
pct_rna <- as.data.frame(
  t(sapply(levels(Idents(TB_ALL)), function(cl) {
    cells <- WhichCells(TB_ALL, idents = cl)
    mat   <- rna_data_mat[rna_available, cells, drop = FALSE]
    rowMeans(mat > 0) * 100
  }))
)
write.csv(round(pct_rna, 2),
          file.path(out_avgexpr, "PctExpr_RNA_per_cluster.csv"), row.names = TRUE)

DefaultAssay(TB_ALL) <- "ADT"
TB_ALL[["ADT"]] <- JoinLayers(TB_ALL[["ADT"]])
adt_data_mat <- GetAssayData(TB_ALL, slot = "data")
pct_adt <- as.data.frame(
  t(sapply(levels(Idents(TB_ALL)), function(cl) {
    cells <- WhichCells(TB_ALL, idents = cl)
    mat   <- adt_data_mat[adt_available, cells, drop = FALSE]
    rowMeans(mat > 0) * 100
  }))
)
write.csv(round(pct_adt, 2),
          file.path(out_avgexpr, "PctExpr_ADT_per_cluster.csv"), row.names = TRUE)

message("Average expression tables saved to: ", out_avgexpr)

###############################################################################
#  FEATURE PLOTS — ADT (all proteins on WNN UMAP)
###############################################################################

message("Generating ADT feature plots...")
DefaultAssay(TB_ALL) <- "ADT"

for (feat in adt_available) {
  fp <- DimPlot2(
    TB_ALL, features = feat, reduction = "wnn.umap"
  ) + ggtitle(paste("ADT |", feat))
  
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", feat)
  ggsave(
    file.path(out_feat_adt, paste0("FeatADT_", safe_name, ".png")),
    fp, dpi = 300, width = 8, height = 6, bg = "white"
  )
}
message("ADT feature plots done: ", length(adt_available))

###############################################################################
#  FEATURE PLOTS — RNA (all marker genes on WNN UMAP)
###############################################################################

message("Generating RNA feature plots...")
DefaultAssay(TB_ALL) <- "RNA"

for (gene in rna_available) {
  if (!gene %in% rownames(TB_ALL[["RNA"]])) next
  
  fp <- DimPlot2(
    TB_ALL, features = gene, reduction = "wnn.umap"
  ) + ggtitle(paste("RNA |", gene))
  
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", gene)
  ggsave(
    file.path(out_feat_rna, paste0("FeatRNA_", safe_name, ".png")),
    fp, dpi = 300, width = 8, height = 6, bg = "white"
  )
}
message("RNA feature plots done: ", length(rna_available))

###############################################################################
#  VIOLIN PLOTS — ADT (all proteins)
###############################################################################

message("Generating ADT violin plots...")
DefaultAssay(TB_ALL) <- "ADT"

for (feat in adt_available) {
  p <- VlnPlot2(
    TB_ALL, features = feat, group.by = "seurat_clusters",
    assay = "ADT", cols = "light", show.mean = TRUE,
    mean_colors = c("red", "blue")
  )
  
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", feat)
  ggsave(
    file.path(out_vln_adt, paste0("VlnADT_", safe_name, ".png")),
    p, width = 14, height = 5, dpi = 300, bg = "white"
  )
}
message("ADT violin plots done: ", length(adt_available))

###############################################################################
#  VIOLIN PLOTS — RNA (all marker genes)
###############################################################################

message("Generating RNA violin plots...")
DefaultAssay(TB_ALL) <- "RNA"

for (gene in rna_available) {
  if (!gene %in% rownames(TB_ALL[["RNA"]])) next
  
  p <- VlnPlot2(
    TB_ALL, features = gene, group.by = "seurat_clusters",
    assay = "RNA", cols = "light", show.mean = TRUE,
    mean_colors = c("red", "blue")
  )
  
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", gene)
  ggsave(
    file.path(out_vln_rna, paste0("VlnRNA_", safe_name, ".png")),
    p, width = 14, height = 5, dpi = 300, bg = "white"
  )
}
message("RNA violin plots done: ", length(rna_available))

###############################################################################
#  HEATMAPS — Per cell-type category (readable, split by lineage)
#
#  Each heatmap shows ALL clusters (columns) × markers for ONE lineage (rows)
#  Separate heatmaps for ADT and RNA per category
###############################################################################

message("Generating per-lineage heatmaps...")

# Define heatmap groups — curated subsets of markers per lineage
heatmap_groups <- list(
  `CD8_T_cell` = list(
    RNA = c("CD8A", "CD8B", "CD3D", "CD3E", "GZMB", "GZMA", "GZMH", "GZMK",
            "PRF1", "GNLY", "NKG7", "KLRG1", "CX3CR1", "FGFBP2",
            "TBX21", "EOMES", "RUNX3", "PRDM1", "ZEB2", "ID2",
            "TOX", "PDCD1", "HAVCR2", "TIGIT", "LAG3", "CTLA4", "ENTPD1",
            "CCR7", "SELL", "TCF7", "LEF1", "IL7R", "BCL2",
            "MKI67", "CD69", "CD38", "HLA-DRA", "B3GAT1"),
    ADT = c("CD8A", "CD3D", "KLRG1", "CX3CR1", "ENTPD1", "PDCD1",
            "TIGIT", "LAG3", "CTLA4", "SELL", "IL7R", "CD27", "CD28",
            "CD38", "FAS", "CD69", "TNFRSF9", "CD45RA", "CD45RO",
            "ICOS", "ITGAE", "ITGA1", "NCAM1", "HLA-DRA")
  ),
  `CD4_T_cell` = list(
    RNA = c("CD4", "IL7R", "FOXP3", "IL2RA", "CTLA4", "IKZF2",
            "CCR4", "RORC", "GATA3", "TBX21", "BCL6", "CXCR5",
            "MAF", "ICOS", "BATF", "CCR7", "SELL", "TCF7", "LEF1"),
    ADT = c("CD4", "IL7R", "IL2RA", "CTLA4", "CCR4", "CCR6",
            "CXCR5", "CXCR3", "ICOS", "SELL", "CD27", "CD28",
            "CD45RA", "CD45RO", "CD69", "HLA-DRA", "TIGIT", "PDCD1")
  ),
  `NK_cell` = list(
    RNA = c("NCAM1", "KLRF1", "KLRC1", "KLRB1", "NCR1", "FCGR3A",
            "TYROBP", "FCER1G", "KLRD1", "KIR2DL3", "KIR3DL1",
            "GZMB", "GZMA", "PRF1", "GNLY", "NKG7",
            "SPON2", "CLIC3", "MYOM2", "SH2D1B"),
    ADT = c("NCAM1", "FCGR3A", "KLRD1", "KLRB1", "KLRK1", "NCR1",
            "KIR2DL1", "KIR2DL3", "KIR3DL1", "CD38", "CD69",
            "ITGAM", "ITGAX", "CD7", "SIGLEC7")
  ),
  `B_cell_Plasma` = list(
    RNA = c("CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "BANK1",
            "TCL1A", "IGHM", "IGHD", "IGKC",
            "SDC1", "XBP1", "MZB1", "JCHAIN", "PRDM1", "IRF4"),
    ADT = c("CD19", "MS4A1", "CD79B", "CR2", "CD22", "FCER2",
            "CD24", "CD27", "IGHD", "IGHM", "IGKC",
            "CD38", "SELL", "HLA-DRA")
  ),
  `Monocyte_DC` = list(
    RNA = c("CD14", "LYZ", "S100A8", "S100A9", "S100A12",
            "FCGR3A", "CSF1R", "CD68", "VCAN", "FCN1", "MNDA",
            "FCER1A", "CLEC10A", "CD1C", "CLEC9A",
            "LILRA4", "CLEC4C", "IRF7", "TCF4", "IL3RA"),
    ADT = c("CD14", "FCGR3A", "ITGAM", "ITGAX", "CD1C",
            "CLEC12A", "THBD", "SIGLEC1", "CD163",
            "HLA-DRA", "HLA-A", "CLEC4C", "IL3RA",
            "CD33", "FCGR1A", "FCGR2A", "CD40")
  ),
  `GammaDelta_MAIT` = list(
    RNA = c("TRDV1", "TRDV2", "TRGV9", "TRDC", "TRGC1", "TRGC2",
            "TRAV1-2", "SLC4A10", "KLRB1", "ZBTB16", "RORC",
            "GZMB", "GZMA", "GNLY", "NKG7", "CD8A", "CD4"),
    ADT = c("CD3D", "CD8A", "CD4", "NCAM1", "KLRB1", "KLRD1",
            "CD38", "CD69", "CD27", "SELL", "CD45RA", "CD45RO",
            "ITGAE", "ITGB7")
  ),
  `Naive_Memory_Exhaustion` = list(
    RNA = c("CCR7", "SELL", "TCF7", "LEF1", "BACH2", "IL7R", "BCL2",
            "CD27", "CD28", "GZMK", "GZMB", "PRF1",
            "TOX", "PDCD1", "HAVCR2", "TIGIT", "LAG3", "CTLA4",
            "MKI67", "B3GAT1", "CX3CR1", "ZEB2"),
    ADT = c("SELL", "IL7R", "CD27", "CD28", "CD45RA", "CD45RO",
            "PDCD1", "TIGIT", "LAG3", "CTLA4", "ENTPD1",
            "KLRG1", "CX3CR1", "CD38", "FAS", "CD69",
            "TNFRSF9", "ICOS", "HLA-DRA")
  ),
  `IFN_Stress` = list(
    RNA = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI44L", "MX1",
            "OAS1", "STAT1", "IRF7",
            "HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1", "IFI27"),
    ADT = c("HLA-A", "HLA-DRA", "HLA-E", "CD38", "FAS", "CD69")
  ),
  `Platelet_RBC_HSC` = list(
    RNA = c("PPBP", "PF4", "GP9", "ITGA2B",
            "HBB", "HBA1", "HBA2", "ALAS2",
            "CD34", "SPINK2", "AVP", "SOX4", "HOPX"),
    ADT = c("ITGA2B", "GP1BB", "SELP", "PECAM1", "CD34", "CD38")
  )
)

# Heatmap color palette
heatmap_colors <- colorRampPalette(c("#F7FCF5", "#C7E9C0", "#74C476", "#31A354", "#006D2C"))(100)

# Generate per-lineage heatmaps
for (group_name in names(heatmap_groups)) {
  group <- heatmap_groups[[group_name]]
  
  # --- RNA heatmap ---
  DefaultAssay(TB_ALL) <- "RNA"
  rna_feats <- intersect(group$RNA, rownames(TB_ALL))
  
  if (length(rna_feats) > 0) {
    avg_rna <- AverageExpression(
      TB_ALL, assays = "RNA", features = rna_feats,
      group.by = "seurat_clusters", slot = "data"
    )$RNA
    
    avg_rna_sc <- scale_01(as.matrix(avg_rna))
    
    # Order columns numerically (strips 'g' prefix)
    avg_rna_sc <- sort_cluster_cols(avg_rna_sc)
    
    png(file.path(out_heat, paste0("Heatmap_RNA_", group_name, ".png")),
        width = max(14, ncol(avg_rna_sc) * 0.7 + 4),
        height = max(6, nrow(avg_rna_sc) * 0.35 + 2),
        units = "in", res = 300, bg = "white")
    
    pheatmap(
      avg_rna_sc,
      cluster_rows  = TRUE,
      cluster_cols  = FALSE,
      scale         = "none",
      color         = heatmap_colors,
      border_color  = "white",
      cellwidth     = 30,
      cellheight    = 14,
      fontsize      = 10,
      fontsize_row  = 9,
      fontsize_col  = 10,
      angle_col     = 0,
      main          = paste0("RNA — ", gsub("_", " ", group_name))
    )
    dev.off()
  }
  
  # --- ADT heatmap ---
  DefaultAssay(TB_ALL) <- "ADT"
  adt_feats <- intersect(group$ADT, rownames(TB_ALL))
  
  if (length(adt_feats) > 0) {
    avg_adt <- AverageExpression(
      TB_ALL, assays = "ADT", features = adt_feats,
      group.by = "seurat_clusters", slot = "data"
    )$ADT
    
    avg_adt_sc <- scale_01(as.matrix(log10(avg_adt + 1)))
    
    # Order columns numerically (strips 'g' prefix)
    avg_adt_sc <- sort_cluster_cols(avg_adt_sc)
    
    png(file.path(out_heat, paste0("Heatmap_ADT_", group_name, ".png")),
        width = max(14, ncol(avg_adt_sc) * 0.7 + 4),
        height = max(6, nrow(avg_adt_sc) * 0.35 + 2),
        units = "in", res = 300, bg = "white")
    
    pheatmap(
      avg_adt_sc,
      cluster_rows  = TRUE,
      cluster_cols  = FALSE,
      scale         = "none",
      color         = heatmap_colors,
      border_color  = "white",
      cellwidth     = 30,
      cellheight    = 14,
      fontsize      = 10,
      fontsize_row  = 9,
      fontsize_col  = 10,
      angle_col     = 0,
      main          = paste0("ADT — ", gsub("_", " ", group_name))
    )
    dev.off()
  }
  
  message("  ", group_name, " heatmaps done")
}

###############################################################################
#  OVERVIEW UMAP — for reference
###############################################################################

p_overview <- DimPlot2(
  TB_ALL,
  reduction  = "wnn.umap",
  group.by   = "seurat_clusters",
  cols       = "light",
  label      = TRUE,
  box        = TRUE,
  repel      = TRUE,
  label.size = 6,
  pt.size    = 0.4,
  raster     = FALSE,
  theme      = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("TB WNN — wsnn_res.0.6")

ggsave(
  file.path(annot_dir, "UMAP_overview_res0.6.png"),
  p_overview, dpi = 300, width = 10, height = 8, bg = "white"
)

# Azimuth reference
p_azimuth <- DimPlot2(
  TB_ALL,
  reduction  = "wnn.umap",
  group.by   = "predicted.celltype.l2",
  label      = TRUE,
  box        = TRUE,
  repel      = TRUE,
  label.size = 4,
  pt.size    = 0.4,
  raster     = FALSE,
  theme      = list(NoAxes(), theme_umap_arrows())
) + ggtitle("TB WNN — Azimuth L2")

ggsave(
  file.path(annot_dir, "UMAP_Azimuth_L2.png"),
  p_azimuth, dpi = 300, width = 12, height = 8, bg = "white"
)

message("\n=== All pre-annotation plots complete ===")
message("Upload the CSVs from ", out_avgexpr, " for annotation help.")
###############################################################################
#  CLUSTER ANNOTATION (wsnn_res.0.6) — APPLIED
#  ADT-first calls with RNA corroboration; see annotation block above for evidence
###############################################################################

message("\n=== Applying cluster annotations ===")

# --- 1. Fix mislabeled IGRA status for sample 6061161 ---
# Sample 6061161's mother is IGRA-negative, not positive
message("Fixing IGRA_mother for sample 6061161 (IGRA+ -> IGRA-)...")

cat("Before fix:\n")
print(table(TB_ALL$Sample_ID, TB_ALL$IGRA_mother))

TB_ALL$IGRA_mother[TB_ALL$Sample_ID == "6061161"] <- "IGRA-"

# If IGRA_mother is a factor, drop unused levels and re-level
if (is.factor(TB_ALL$IGRA_mother)) {
  TB_ALL$IGRA_mother <- droplevels(TB_ALL$IGRA_mother)
}

cat("\nAfter fix:\n")
print(table(TB_ALL$Sample_ID, TB_ALL$IGRA_mother))

# --- 2. Backup original cluster identity ---
TB_ALL$wsnn_res.0.6_orig <- Idents(TB_ALL)

# --- 3. Fine-grained annotation (per cluster) ---
cluster_annotations <- c(
  "0"  = "Naive CD4 T",
  "1"  = "Naive B (follicular)",
  "2"  = "Naive CD8 T",
  "3"  = "CD56dim NK",
  "4"  = "Naive CD4 T",
  "5"  = "Naive CD4 T (IL7R-hi)",
  "6"  = "TEMRA CD8",
  "7"  = "CM CD4 (Th2/CCR4+)",
  "8"  = "Transitional B (T1/T2)",
  "9"  = "Activated Naive CD4",
  "10" = "Naive CD4 T (low-signal)",
  "11" = "Memory B (CD27- / MZ-like)",
  "12" = "Treg",
  "13" = "Classical Monocyte",
  "14" = "Platelet doublets",
  "15" = "Vd1 gd T",
  "16" = "Activated B / GC-pre",
  "17" = "MAIT + Vg9Vd2 (mixed)",
  "18" = "GZMK+ Tem CD8",
  "19" = "CD56bright NK",
  "20" = "Naive CD4 T (SOX4-hi)",
  "21" = "DC mix (pDC+cDC2+prog)",
  "22" = "ILC (ILC3-like)",
  "23" = "Doublets (T+B+RBC)",
  "24" = "Low-quality (ADT artifact)",
  "25" = "Proliferating Plasmablast"
)

# Fine annotation — strip names so Seurat doesn't try to match them to barcodes
TB_ALL$celltype_fine <- unname(cluster_annotations[as.character(Idents(TB_ALL))])
TB_ALL$celltype_fine <- factor(TB_ALL$celltype_fine, levels = unique(cluster_annotations))


# --- 4. Collapsed annotation (merging Naive CD4 family + cleaner top-level labels) ---
collapsed_map <- c(
  "Naive CD4 T"                 = "Naive CD4 T",
  "Naive CD4 T (IL7R-hi)"       = "Naive CD4 T",
  "Naive CD4 T (low-signal)"    = "Naive CD4 T",
  "Naive CD4 T (SOX4-hi)"       = "Naive CD4 T",
  "Activated Naive CD4"         = "Activated CD4 T",
  "CM CD4 (Th2/CCR4+)"          = "CM CD4 T",
  "Treg"                        = "Treg",
  "Naive CD8 T"                 = "Naive CD8 T",
  "TEMRA CD8"                   = "TEMRA CD8 T",
  "GZMK+ Tem CD8"               = "Tem CD8 T",
  "Vd1 gd T"                    = "gd T",
  "MAIT + Vg9Vd2 (mixed)"       = "Innate-like T (MAIT/gd)",
  "CD56dim NK"                  = "CD56dim NK",
  "CD56bright NK"               = "CD56bright NK",
  "ILC (ILC3-like)"             = "ILC",
  "Naive B (follicular)"        = "Naive B",
  "Transitional B (T1/T2)"      = "Transitional B",
  "Memory B (CD27- / MZ-like)"  = "Memory B",
  "Activated B / GC-pre"        = "Activated B",
  "Proliferating Plasmablast"   = "Plasmablast",
  "Classical Monocyte"          = "Classical Mono",
  "DC mix (pDC+cDC2+prog)"      = "DC",
  "Platelet doublets"           = "Platelet",
  "Doublets (T+B+RBC)"          = "Doublets",
  "Low-quality (ADT artifact)"  = "Low-quality"
)
# Collapsed annotation — same fix
TB_ALL$celltype_main <- unname(collapsed_map[as.character(TB_ALL$celltype_fine)])
TB_ALL$celltype_main <- factor(TB_ALL$celltype_main, levels = unique(collapsed_map))

# --- 5. Mark cells to drop before downstream (doublets / QC failures) ---
drop_clusters <- c("14", "23", "24")
TB_ALL$keep <- !(as.character(Idents(TB_ALL)) %in% drop_clusters)

cat("\n--- Cells flagged for removal ---\n")
cat("Dropping clusters:", paste(drop_clusters, collapse = ", "), "\n")
cat("Cells to drop:", sum(!TB_ALL$keep),
    "(", round(100 * sum(!TB_ALL$keep) / ncol(TB_ALL), 2), "%)\n")
cat("Cells to keep:", sum(TB_ALL$keep), "\n\n")

# Summary tables
cat("--- celltype_fine counts ---\n")
print(table(TB_ALL$celltype_fine))
cat("\n--- celltype_main counts ---\n")
print(table(TB_ALL$celltype_main))

###############################################################################
#  ANNOTATED UMAPs
###############################################################################

annot_umap_dir <- file.path(annot_dir, "Annotated_UMAPs")
dir.create(annot_umap_dir, recursive = TRUE, showWarnings = FALSE)

# Fine annotation
p_fine <- DimPlot2(
  TB_ALL,
  reduction  = "wnn.umap",
  group.by   = "celltype_fine",
  cols       = "light",
  label      = TRUE,
  box        = TRUE,
  repel      = TRUE,
  label.size = 3.5,
  pt.size    = 0.4,
  raster     = FALSE,
  theme      = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("TB WNN — Fine annotation (per cluster)")

ggsave(
  file.path(annot_umap_dir, "UMAP_celltype_fine.png"),
  p_fine, dpi = 300, width = 12, height = 9, bg = "white"
)

# Main (collapsed) annotation
p_main <- DimPlot2(
  TB_ALL,
  reduction  = "wnn.umap",
  group.by   = "celltype_main",
  cols       = "light",
  label      = TRUE,
  box        = TRUE,
  repel      = TRUE,
  label.size = 4.5,
  pt.size    = 0.4,
  raster     = FALSE,
  theme      = list(NoAxes(), theme_umap_arrows())
) + ggtitle("TB WNN — Main annotation (collapsed)")

ggsave(
  file.path(annot_umap_dir, "UMAP_celltype_main.png"),
  p_main, dpi = 300, width = 13, height = 9, bg = "white"
)

# QC-pass cells highlighted
p_keep <- DimPlot2(
  TB_ALL,
  reduction  = "wnn.umap",
  group.by   = "keep",
  cols       = c("TRUE" = "#2E86AB", "FALSE" = "#E63946"),
  label      = FALSE,
  pt.size    = 0.4,
  raster     = FALSE,
  theme      = list(NoAxes(), theme_umap_arrows())
) + ggtitle("Cells flagged for removal (clusters 14, 23, 24)")

ggsave(
  file.path(annot_umap_dir, "UMAP_keep_flag.png"),
  p_keep, dpi = 300, width = 11, height = 8, bg = "white"
)

###############################################################################
#  EXPORT ANNOTATION TABLE
###############################################################################

annot_table <- data.frame(
  cluster        = names(cluster_annotations),
  celltype_fine  = unname(cluster_annotations),
  celltype_main  = unname(collapsed_map[cluster_annotations]),
  n_cells        = as.integer(table(TB_ALL$wsnn_res.0.6_orig)[names(cluster_annotations)]),
  flagged_drop   = names(cluster_annotations) %in% drop_clusters
)

write.csv(annot_table,
          file.path(annot_dir, "Cluster_annotation_table.csv"),
          row.names = FALSE)

cat("\nAnnotation table saved to: ",
    file.path(annot_dir, "Cluster_annotation_table.csv"), "\n")

###############################################################################
#  SAVE ANNOTATED OBJECTS
###############################################################################

save.path <- file.path(base_dir, "saved_R_data")
dir.create(save.path, recursive = TRUE, showWarnings = FALSE)

# Full annotated object (all cells, including flagged ones)
message("\nSaving full annotated object (TB_ALL_Annotated.qs2)...")
qs_save(TB_ALL, file.path(save.path, "TB_ALL_Annotated.qs2"))

# Clean object — flagged doublet/QC-fail clusters removed
message("Subsetting to drop QC-fail clusters and saving TB_ALL_Annotated_Clean.qs2...")
TB_clean <- subset(TB_ALL, subset = keep)

# Drop unused factor levels in cleaned object
TB_clean$celltype_fine <- droplevels(TB_clean$celltype_fine)
TB_clean$celltype_main <- droplevels(TB_clean$celltype_main)
TB_clean$wsnn_res.0.6_orig <- droplevels(TB_clean$wsnn_res.0.6_orig)

qs_save(TB_clean, file.path(save.path, "TB_ALL_Annotated_Clean.qs2"))

cat("\nSaved objects:\n")
cat("  Full:  ", file.path(save.path, "TB_ALL_Annotated.qs2"), "  (", ncol(TB_ALL), " cells)\n", sep = "")
cat("  Clean: ", file.path(save.path, "TB_ALL_Annotated_Clean.qs2"), " (", ncol(TB_clean), " cells)\n", sep = "")

message("\n=== Annotation pipeline complete ===")
