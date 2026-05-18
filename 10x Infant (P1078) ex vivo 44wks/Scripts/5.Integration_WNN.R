###############################################################################
#  TB CITE-seq Pipeline — Script 5: RNA + ADT Integration + WNN + Clustree
#  Adapted from OPIS CITE-seq pipeline
###############################################################################

## --------------------------
## 0) Libraries
## --------------------------
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(SeuratWrappers)
library(Azimuth)
library(ggrepel)
library(patchwork)
library(SeuratExtend)    # DimPlot2
library(batchelor)
library(harmony)
library(qs2)
library(clustree)

## --------------------------
## 1) Paths & Folders
## --------------------------
base_dir   <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/TB_CITEseq"
load.path  <- file.path(base_dir, "saved_R_data")

integration_dir <- file.path(base_dir, "Integration")
plot_dir        <- file.path(base_dir, "Integration", "Plots")
clustree_dir    <- file.path(base_dir, "Integration", "Clustree")

dir.create(integration_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,        recursive = TRUE, showWarnings = FALSE)
dir.create(clustree_dir,    recursive = TRUE, showWarnings = FALSE)

## --------------------------
## 2) Load Isotype-Corrected Object
## --------------------------
seurat_isotype <- qs_read(file.path(load.path, "TB_Seurat_isotype.qs2"))

## --------------------------
## 3) Prep — Join RNA, define TB_ALL
## --------------------------
seurat_isotype[["RNA"]] <- JoinLayers(seurat_isotype[["RNA"]])
TB_ALL <- seurat_isotype
rm(seurat_isotype); gc()

# Split RNA layers per sample for integration
TB_ALL[["RNA"]] <- split(TB_ALL[["RNA"]], f = TB_ALL$Sample_ID)

## --------------------------
## 4) RNA — Unintegrated + Integrated (CCA + FastMNN)
## --------------------------
DefaultAssay(TB_ALL) <- "RNA"

TB_ALL <- NormalizeData(TB_ALL)
TB_ALL <- FindVariableFeatures(TB_ALL)
TB_ALL <- ScaleData(TB_ALL)
TB_ALL <- RunPCA(TB_ALL)

# Unintegrated
TB_ALL <- FindNeighbors(TB_ALL, dims = 1:30, reduction = "pca")
TB_ALL <- FindClusters(TB_ALL, resolution = 2, cluster.name = "unintegrated_clusters")
TB_ALL <- RunUMAP(
  TB_ALL,
  dims           = 1:30,
  reduction      = "pca",
  reduction.name = "umap.unintegrated"
)

# CCA integration
TB_ALL <- IntegrateLayers(
  TB_ALL,
  method         = CCAIntegration,
  orig.reduction = "pca",
  assay          = "RNA",
  new.reduction  = "integrated.cca.rna"
)

# FastMNN integration
TB_ALL <- IntegrateLayers(
  TB_ALL,
  method        = FastMNNIntegration,
  assay         = "RNA",
  new.reduction = "integrated.mnn.rna"
)

# CCA — UMAP + clusters
TB_ALL <- FindNeighbors(TB_ALL, reduction = "integrated.cca.rna", dims = 1:30)
TB_ALL <- FindClusters(TB_ALL, resolution = 2, cluster.name = "cca_clusters_rna")
TB_ALL <- RunUMAP(
  TB_ALL,
  reduction      = "integrated.cca.rna",
  dims           = 1:30,
  reduction.name = "umap.cca.rna"
)

# FastMNN — UMAP + clusters
TB_ALL <- FindNeighbors(TB_ALL, reduction = "integrated.mnn.rna", dims = 1:30)
TB_ALL <- FindClusters(TB_ALL, resolution = 2, cluster.name = "mnn_clusters_rna")
TB_ALL <- RunUMAP(
  TB_ALL,
  reduction      = "integrated.mnn.rna",
  dims           = 1:30,
  reduction.name = "umap.mnn.rna"
)

## --------------------------
## 5) RNA Integration Plots
## --------------------------
p_rna_cca_celltype <- DimPlot2(
  TB_ALL, reduction = "umap.cca.rna", group.by = "predicted.celltype.l2",
  label = TRUE, box = TRUE, repel = TRUE, label.size = 3, pt.size = 0.4,
  theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("RNA CCA — Azimuth L2")

p_rna_cca_cluster <- DimPlot2(
  TB_ALL, reduction = "umap.cca.rna", group.by = "cca_clusters_rna",
  label = TRUE, box = TRUE, repel = TRUE, label.size = 3, pt.size = 0.4,
  theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("RNA CCA — Clusters")

p_rna_cca_sample <- DimPlot2(
  TB_ALL, reduction = "umap.cca.rna", group.by = "Sample_ID",
  label = TRUE, box = TRUE, repel = TRUE, label.size = 3, pt.size = 0.4,
  theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("RNA CCA — Sample")

p_rna_mnn_celltype <- DimPlot2(
  TB_ALL, reduction = "umap.mnn.rna", group.by = "predicted.celltype.l2",
  label = TRUE, box = TRUE, repel = TRUE, label.size = 3, pt.size = 0.4,
  theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("RNA MNN — Azimuth L2")

p_rna_mnn_cluster <- DimPlot2(
  TB_ALL, reduction = "umap.mnn.rna", group.by = "mnn_clusters_rna",
  label = TRUE, box = TRUE, repel = TRUE, label.size = 3, pt.size = 0.4,
  theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("RNA MNN — Clusters")

p_rna_mnn_sample <- DimPlot2(
  TB_ALL, reduction = "umap.mnn.rna", group.by = "Sample_ID",
  label = TRUE, box = TRUE, repel = TRUE, label.size = 3, pt.size = 0.4,
  theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("RNA MNN — Sample")

combined_rna <- wrap_plots(
  p_rna_cca_celltype, p_rna_cca_cluster, p_rna_cca_sample,
  p_rna_mnn_celltype, p_rna_mnn_cluster, p_rna_mnn_sample,
  ncol = 3
)

ggsave(
  filename = file.path(plot_dir, "TB_ALL_RNA_Integration.png"),
  plot = combined_rna, width = 26, height = 14, dpi = 300, bg = "white"
)

## --------------------------
## 6) ADT Integration
## --------------------------
DefaultAssay(TB_ALL) <- "ADT"

TB_ALL[["ADT"]] <- as(object = TB_ALL[["ADT"]], Class = "Assay5")
TB_ALL[["ADT"]] <- split(TB_ALL[["ADT"]], f = TB_ALL$Sample_ID)

# Variable features = all ADT proteins minus isotypes
prots <- rownames(TB_ALL[["ADT"]]$data)
isotype_genes <- c("Mouse-IgG1", "Mouse-IgG2a", "Mouse-IgG2b",
                   "Rat-IgG2b", "Rat-IgG1", "Rat-IgG2a", "Hamster-IgG")
prots <- setdiff(prots, isotype_genes)
VariableFeatures(TB_ALL) <- prots

TB_ALL <- ScaleData(TB_ALL)
TB_ALL <- RunPCA(TB_ALL, reduction.name = "apca")

# Unintegrated ADT
TB_ALL <- FindNeighbors(TB_ALL, dims = 1:30, reduction = "apca")
TB_ALL <- FindClusters(TB_ALL, resolution = 2, cluster.name = "unintegrated_clusters_adt")
TB_ALL <- RunUMAP(
  TB_ALL,
  dims           = 1:30,
  reduction      = "apca",
  reduction.name = "umap.unintegrated.adt"
)

# FastMNN on ADT
TB_ALL <- IntegrateLayers(
  TB_ALL,
  orig.reduction = "apca",
  features       = prots,
  method         = FastMNNIntegration,
  assay          = "ADT",
  new.reduction  = "integrated.mnn.adt"
)

# CCA on ADT
TB_ALL <- IntegrateLayers(
  TB_ALL,
  method         = CCAIntegration,
  orig.reduction = "apca",
  features       = prots,
  assay          = "ADT",
  new.reduction  = "integrated.cca.adt"
)

# CCA ADT — UMAP + clusters
TB_ALL <- FindNeighbors(TB_ALL, reduction = "integrated.cca.adt", dims = 1:30)
TB_ALL <- FindClusters(TB_ALL, resolution = 2, cluster.name = "cca_clusters_adt")
TB_ALL <- RunUMAP(
  TB_ALL,
  reduction      = "integrated.cca.adt",
  dims           = 1:30,
  reduction.name = "umap.cca.adt"
)

# FastMNN ADT — UMAP + clusters
TB_ALL <- FindNeighbors(TB_ALL, reduction = "integrated.mnn.adt", dims = 1:30)
TB_ALL <- FindClusters(TB_ALL, resolution = 2, cluster.name = "mnn_clusters_adt")
TB_ALL <- RunUMAP(
  TB_ALL,
  reduction      = "integrated.mnn.adt",
  dims           = 1:30,
  reduction.name = "umap.mnn.adt"
)

## --------------------------
## 7) ADT Integration Plots
## --------------------------
p_adt_cca_celltype <- DimPlot2(
  TB_ALL, reduction = "umap.cca.adt", group.by = "predicted.celltype.l2",
  label = TRUE, box = TRUE, repel = TRUE, label.size = 3, pt.size = 0.4,
  theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("ADT CCA — Azimuth L2")

p_adt_cca_cluster <- DimPlot2(
  TB_ALL, reduction = "umap.cca.adt", group.by = "cca_clusters_adt",
  label = TRUE, box = TRUE, repel = TRUE, label.size = 3, pt.size = 0.4,
  theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("ADT CCA — Clusters")

p_adt_cca_sample <- DimPlot2(
  TB_ALL, reduction = "umap.cca.adt", group.by = "Sample_ID",
  label = TRUE, box = TRUE, repel = TRUE, label.size = 3, pt.size = 0.4,
  theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("ADT CCA — Sample")

p_adt_mnn_celltype <- DimPlot2(
  TB_ALL, reduction = "umap.mnn.adt", group.by = "predicted.celltype.l2",
  label = TRUE, box = TRUE, repel = TRUE, label.size = 3, pt.size = 0.4,
  theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("ADT MNN — Azimuth L2")

p_adt_mnn_cluster <- DimPlot2(
  TB_ALL, reduction = "umap.mnn.adt", group.by = "mnn_clusters_adt",
  label = TRUE, box = TRUE, repel = TRUE, label.size = 3, pt.size = 0.4,
  theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("ADT MNN — Clusters")

p_adt_mnn_sample <- DimPlot2(
  TB_ALL, reduction = "umap.mnn.adt", group.by = "Sample_ID",
  label = TRUE, box = TRUE, repel = TRUE, label.size = 3, pt.size = 0.4,
  theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("ADT MNN — Sample")

combined_adt <- wrap_plots(
  p_adt_cca_celltype, p_adt_cca_cluster, p_adt_cca_sample,
  p_adt_mnn_celltype, p_adt_mnn_cluster, p_adt_mnn_sample,
  ncol = 3
)

ggsave(
  filename = file.path(plot_dir, "TB_ALL_ADT_Integration.png"),
  plot = combined_adt, width = 26, height = 14, dpi = 300, bg = "white"
)

## --------------------------
## 8) Save Intermediate — RNA + ADT integrated
## --------------------------
qs_save(TB_ALL, file = file.path(load.path, "TB_ALL_RNA_ADT.qs2"))
message("Saved: TB_ALL_RNA_ADT.qs2")

## --------------------------
## 9) WNN (Weighted Nearest Neighbors) — using CCA reductions
## --------------------------
TB_ALL <- FindMultiModalNeighbors(
  TB_ALL,
  reduction.list       = list("integrated.cca.rna", "integrated.cca.adt"),
  dims.list            = list(1:30, 1:20),
  modality.weight.name = "wnn.weight"
)

TB_ALL <- RunUMAP(
  TB_ALL,
  nn.name        = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key  = "wnnUMAP_"
)

# Multiple clustering algorithms × resolutions
TB_ALL <- FindClusters(TB_ALL, graph.name = "wsnn", algorithm = 2, resolution = 0.5, cluster.name = "snn.louvianmlr_0.5")
TB_ALL <- FindClusters(TB_ALL, graph.name = "wsnn", algorithm = 2, resolution = 1,   cluster.name = "snn.louvianmlr_1")
TB_ALL <- FindClusters(TB_ALL, graph.name = "wsnn", algorithm = 2, resolution = 1.5, cluster.name = "snn.louvianmlr_1.5")

TB_ALL <- FindClusters(TB_ALL, graph.name = "wsnn", algorithm = 3, resolution = 0.5, cluster.name = "snn.slm_0.5")
TB_ALL <- FindClusters(TB_ALL, graph.name = "wsnn", algorithm = 3, resolution = 1,   cluster.name = "snn.slm_1")
TB_ALL <- FindClusters(TB_ALL, graph.name = "wsnn", algorithm = 3, resolution = 1.5, cluster.name = "snn.slm_1.5")

# Order cluster factor levels numerically
cluster_cols <- c(
  "snn.louvianmlr_0.5", "snn.louvianmlr_1", "snn.louvianmlr_1.5",
  "snn.slm_0.5", "snn.slm_1", "snn.slm_1.5"
)

for (cluster_col in cluster_cols) {
  if (cluster_col %in% colnames(TB_ALL@meta.data)) {
    vals     <- TB_ALL[[cluster_col]][, 1]
    levs_num <- suppressWarnings(as.numeric(levels(vals)))
    if (!any(is.na(levs_num))) {
      TB_ALL[[cluster_col]][, 1] <- factor(vals, levels = sort(levs_num))
    } else {
      TB_ALL[[cluster_col]][, 1] <- factor(vals)
    }
  }
}

## --------------------------
## 10) WNN Plots
## --------------------------
p1 <- DimPlot2(TB_ALL, reduction = "wnn.umap", group.by = "predicted.celltype.l1",
               label = TRUE, box = TRUE, repel = TRUE, label.size = 5, pt.size = 0.4,
               theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("Azimuth L1")

p2 <- DimPlot2(TB_ALL, reduction = "wnn.umap", group.by = "predicted.celltype.l2",
               label = TRUE, box = TRUE, repel = TRUE, label.size = 5, pt.size = 0.4,
               theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("Azimuth L2")

p3 <- DimPlot2(TB_ALL, reduction = "wnn.umap", group.by = "predicted.celltype.l3",
               label = TRUE, box = TRUE, repel = TRUE, label.size = 5, pt.size = 0.4,
               theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("Azimuth L3")

p4 <- DimPlot2(TB_ALL, reduction = "wnn.umap", group.by = "snn.slm_0.5",
               label = TRUE, box = TRUE, repel = TRUE, label.size = 5, pt.size = 0.4,
               theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("SLM res 0.5")

p5 <- DimPlot2(TB_ALL, reduction = "wnn.umap", group.by = "snn.slm_1",
               label = TRUE, box = TRUE, repel = TRUE, label.size = 5, pt.size = 0.4,
               theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("SLM res 1.0")

p6 <- DimPlot2(TB_ALL, reduction = "wnn.umap", group.by = "snn.slm_1.5",
               label = TRUE, box = TRUE, repel = TRUE, label.size = 5, pt.size = 0.4,
               theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("SLM res 1.5")

p7 <- DimPlot2(TB_ALL, reduction = "wnn.umap", group.by = "snn.louvianmlr_0.5",
               label = TRUE, box = TRUE, repel = TRUE, label.size = 5, pt.size = 0.4,
               theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("Louvain MLR res 0.5")

p8 <- DimPlot2(TB_ALL, reduction = "wnn.umap", group.by = "snn.louvianmlr_1",
               label = TRUE, box = TRUE, repel = TRUE, label.size = 5, pt.size = 0.4,
               theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("Louvain MLR res 1.0")

p9 <- DimPlot2(TB_ALL, reduction = "wnn.umap", group.by = "snn.louvianmlr_1.5",
               label = TRUE, box = TRUE, repel = TRUE, label.size = 5, pt.size = 0.4,
               theme = list(NoLegend(), NoAxes(), theme_umap_arrows())
) + ggtitle("Louvain MLR res 1.5")

# Sample + IGRA splits on WNN
p_sample <- DimPlot2(TB_ALL, reduction = "wnn.umap", group.by = "Sample_ID",
                     label = TRUE, box = TRUE, repel = TRUE, label.size = 4, pt.size = 0.4,
                     theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("Sample ID")

p_igra_infant <- DimPlot2(TB_ALL, reduction = "wnn.umap", group.by = "IGRA_infant",
                          label = TRUE, box = TRUE, repel = TRUE, label.size = 5, pt.size = 0.4,
                          theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("Infant IGRA")

p_igra_mother <- DimPlot2(TB_ALL, reduction = "wnn.umap", group.by = "IGRA_mother",
                          label = TRUE, box = TRUE, repel = TRUE, label.size = 5, pt.size = 0.4,
                          theme = list(NoAxes(), theme_umap_arrows())
) + ggtitle("Maternal IGRA")

combined_wnn <- wrap_plots(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  p_sample, p_igra_infant, p_igra_mother,
  ncol = 3
)

ggsave(
  filename = file.path(plot_dir, "TB_ALL_WNN_Clusters.png"),
  plot = combined_wnn, width = 26, height = 28, dpi = 300, bg = "white"
)

## --------------------------
## 11) Clustree — Resolution Selection
## --------------------------
setwd(clustree_dir)

# Run clustering across a range of resolutions (SLM / algorithm 3)
resolutions <- seq(0.2, 2.0, by = 0.2)

for (res in resolutions) {
  TB_ALL <- FindClusters(
    TB_ALL,
    graph.name = "wsnn",
    resolution = res,
    algorithm  = 3   # SLM (Leiden-like)
  )
}

# Order factor levels numerically for clustree
for (col in grep("^wsnn_res\\.", colnames(TB_ALL@meta.data), value = TRUE)) {
  lvls <- as.character(sort(as.numeric(levels(TB_ALL[[col]][, 1]))))
  TB_ALL[[col]] <- factor(TB_ALL[[col]][, 1], levels = lvls)
}

# Generate clustree
p_clustree <- clustree(TB_ALL, prefix = "wsnn_res.")

ggsave(
  filename = file.path(clustree_dir, "TB_clustree.png"),
  plot = p_clustree, width = 15, height = 9, dpi = 300, bg = "white"
)

message("Clustree saved. Review to choose optimal resolution.")

## --------------------------
## 12) Final Resolution Plot (update after reviewing clustree)
## --------------------------

# ** PLACEHOLDER: Update "wsnn_res.0.6" to your chosen resolution **
chosen_res <- "wsnn_res.0.6"

if (chosen_res %in% colnames(TB_ALL@meta.data)) {
  p_final <- DimPlot2(
    TB_ALL,
    reduction  = "wnn.umap",
    group.by   = chosen_res,
    cols       = "default",
    label      = TRUE,
    box        = TRUE,
    repel      = TRUE,
    label.color = "black",
    label.size = 6,
    pt.size    = 1,
    raster     = FALSE,
    theme      = list(NoLegend(), NoAxes(), theme_umap_arrows())
  ) + ggtitle(paste0("TB WNN — ", chosen_res))
  
  ggsave(
    filename = file.path(plot_dir, paste0("TB_WNN_", chosen_res, ".png")),
    plot = p_final, dpi = 300, bg = "white"
  )
} else {
  message("Resolution column '", chosen_res, "' not found. ",
          "Run FindClusters with that resolution first.")
}

## --------------------------
## 13) Save Final WNN Object
## --------------------------
qs_save(TB_ALL, file = file.path(load.path, "TB_ALL_WNN.qs2"))
message("Done! Saved: TB_ALL_WNN.qs2")
