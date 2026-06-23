# ============================
# Integration for BCG Parse Dataset
# ============================

library(Seurat)
library(harmony)         # Batch correction
library(patchwork)
library(scCustomize)
library(SeuratWrappers)
library(SeuratExtend)
library(qs2)

# ----------------------------
# Base directory (BCG project)
# ----------------------------
base_dir <- "/home/akshay-iyer/Documents/BCG_Stim_Infants_12_44_wks_PARSE_GEX"

# Integration plots folder
integration_dir <- file.path(base_dir, "Integration_Plots")
dir.create(integration_dir, showWarnings = FALSE)

# ----------------------------
# Load doublet-clean, CC-scored object
# ----------------------------
seu_parsef_clean <- qs_read(
  file.path(base_dir, "saved_R_data", "BCG_PARSE_postQC_CellCycle_DoubletClean.qs2")
)

# ----------------------------
# Prepare RNA assay for integration
# ----------------------------
# Join layers first (safety)
seu_parsef_clean[["RNA"]] <- JoinLayers(seu_parsef_clean[["RNA"]])

# Split RNA assay by library (orig.ident = sample_id_timepoint_stim)
# e.g., "6044931_week12_BCG_treated"
seu_parsef_clean[["RNA"]] <- split(
  seu_parsef_clean[["RNA"]],
  f = seu_parsef_clean$orig.ident
)

# ----------------------------
# Standard preprocessing per sample
# ----------------------------
DefaultAssay(seu_parsef_clean) <- "RNA"

seu_parsef_clean <- NormalizeData(seu_parsef_clean)
seu_parsef_clean <- FindVariableFeatures(seu_parsef_clean)
seu_parsef_clean <- ScaleData(seu_parsef_clean)
seu_parsef_clean <- RunPCA(seu_parsef_clean)

seu_parsef_clean <- FindNeighbors(seu_parsef_clean, dims = 1:30, reduction = "pca")
seu_parsef_clean <- FindClusters(
  seu_parsef_clean,
  resolution   = 2,
  cluster.name = "unintegrated_clusters"
)
seu_parsef_clean <- RunUMAP(
  seu_parsef_clean,
  dims           = 1:30,
  reduction      = "pca",
  reduction.name = "umap.unintegrated"
)

# ----------------------------
# Integrate data
# ----------------------------

## 1) CCA Integration
seu_parsef_clean <- IntegrateLayers(
  object        = seu_parsef_clean,
  method        = CCAIntegration,
  orig.reduction = "pca",
  assay         = "RNA",
  new.reduction = "integrated.cca.rna"
)

## 2) Harmony Integration (on PCA)
seu_parsef_clean <- IntegrateLayers(
  object        = seu_parsef_clean,
  method        = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony"
)

## 3) FastMNN Integration
seu_parsef_clean <- IntegrateLayers(
  object        = seu_parsef_clean,
  method        = FastMNNIntegration,
  assay         = "RNA",
  new.reduction = "integrated.mnn.rna"
)

# ----------------------------
# CCA-based Neighbors, UMAP, Clustering
# ----------------------------

seu_parsef_clean <- FindNeighbors(
  seu_parsef_clean,
  reduction = "integrated.cca.rna",
  dims      = 1:30
)
seu_parsef_clean <- RunUMAP(
  seu_parsef_clean,
  reduction      = "integrated.cca.rna",
  dims           = 1:30,
  reduction.name = "umap.cca.rna"
)

# LouvainMLR clustering (algorithm = 2)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 2,
                                 resolution  = 0.5,
                                 cluster.name = "cca.snn.louvianmlr_0.5"
)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 2,
                                 resolution  = 1,
                                 cluster.name = "cca.snn.louvianmlr_1"
)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 2,
                                 resolution  = 1.5,
                                 cluster.name = "cca.snn.louvianmlr_1.5"
)

# SLM clustering (algorithm = 3)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 3,
                                 resolution  = 0.5,
                                 cluster.name = "cca.snn.slm_0.5"
)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 3,
                                 resolution  = 1,
                                 cluster.name = "cca.snn.slm_1"
)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 3,
                                 resolution  = 1.5,
                                 cluster.name = "cca.snn.slm_1.5"
)

# Quick check plot (can comment out if not needed)
DimPlot2(
  seu_parsef_clean,
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  reduction   = "umap.cca.rna",
  group.by    = "cca.snn.louvianmlr_0.5",
  cols        = "default",
  pt.size     = 1
)

# ---------------------------- #
#   CCA Integration Plots      #
# ---------------------------- #

# Row 1: Azimuth Predictions
p1 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.cca.rna",
  group.by    = "predicted.celltype.l1",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("Azimuth L1")

p2 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.cca.rna",
  group.by    = "predicted.celltype.l2",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("Azimuth L2")

p3 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.cca.rna",
  group.by    = "predicted.celltype.l3",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("Azimuth L3")

# Row 2: LouvainMLR
p4 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.cca.rna",
  group.by    = "cca.snn.louvianmlr_0.5",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("cca.snn.louvianmlr_0.5")

p5 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.cca.rna",
  group.by    = "cca.snn.louvianmlr_1",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("cca.snn.louvianmlr_1")

p6 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.cca.rna",
  group.by    = "cca.snn.louvianmlr_1.5",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("cca.snn.louvianmlr_1.5")

# Row 3: SLM
p7 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.cca.rna",
  group.by    = "cca.snn.slm_0.5",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("cca.snn.slm_0.5")

p8 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.cca.rna",
  group.by    = "cca.snn.slm_1",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("cca.snn.slm_1")

p9 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.cca.rna",
  group.by    = "cca.snn.slm_1.5",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("cca.snn.slm_1.5")

combined_plot_cca <- wrap_plots(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  ncol = 3, nrow = 3
)

ggsave(
  filename = file.path(integration_dir, "CCA_Clustering_Grid.png"),
  plot     = combined_plot_cca,
  width    = 24,
  height   = 17,
  dpi      = 300
)

# ============================ #
#       FastMNN Clustering     #
# ============================ #

seu_parsef_clean <- FindNeighbors(
  seu_parsef_clean,
  reduction = "integrated.mnn.rna",
  dims      = 1:30
)
seu_parsef_clean <- RunUMAP(
  seu_parsef_clean,
  reduction      = "integrated.mnn.rna",
  dims           = 1:30,
  reduction.name = "umap.mnn.rna"
)

# LouvainMLR
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 2,
                                 resolution  = 0.5,
                                 cluster.name = "mnn.snn.louvianmlr_0.5"
)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 2,
                                 resolution  = 1,
                                 cluster.name = "mnn.snn.louvianmlr_1"
)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 2,
                                 resolution  = 1.5,
                                 cluster.name = "mnn.snn.louvianmlr_1.5"
)

# SLM
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 3,
                                 resolution  = 0.5,
                                 cluster.name = "mnn.snn.slm_0.5"
)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 3,
                                 resolution  = 1,
                                 cluster.name = "mnn.snn.slm_1"
)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 3,
                                 resolution  = 1.5,
                                 cluster.name = "mnn.snn.slm_1.5"
)

# ---------------------------- #
#   FastMNN Integration Plots  #
# ---------------------------- #

p1 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.mnn.rna",
  group.by    = "predicted.celltype.l1",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("Azimuth L1")

p2 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.mnn.rna",
  group.by    = "predicted.celltype.l2",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("Azimuth L2")

p3 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.mnn.rna",
  group.by    = "predicted.celltype.l3",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("Azimuth L3")

p4 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.mnn.rna",
  group.by    = "mnn.snn.louvianmlr_0.5",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("mnn.snn.louvianmlr_0.5")

p5 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.mnn.rna",
  group.by    = "mnn.snn.louvianmlr_1",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("mnn.snn.louvianmlr_1")

p6 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.mnn.rna",
  group.by    = "mnn.snn.louvianmlr_1.5",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("mnn.snn.louvianmlr_1.5")

p7 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.mnn.rna",
  group.by    = "mnn.snn.slm_0.5",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("mnn.snn.slm_0.5")

p8 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.mnn.rna",
  group.by    = "mnn.snn.slm_1",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("mnn.snn.slm_1")

p9 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.mnn.rna",
  group.by    = "mnn.snn.slm_1.5",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("mnn.snn.slm_1.5")

combined_plot_mnn <- wrap_plots(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  ncol = 3, nrow = 3
)

ggsave(
  filename = file.path(integration_dir, "MNN_Clustering_Grid.png"),
  plot     = combined_plot_mnn,
  width    = 24,
  height   = 17,
  dpi      = 300
)

# ============================ #
#       Harmony Integration    #
# ============================ #

seu_parsef_clean <- FindNeighbors(
  seu_parsef_clean,
  reduction = "harmony",
  dims      = 1:30
)
seu_parsef_clean <- RunUMAP(
  seu_parsef_clean,
  reduction      = "harmony",
  dims           = 1:30,
  reduction.name = "umap.harmony"
)

# LouvainMLR
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 2,
                                 resolution  = 0.5,
                                 cluster.name = "harmony.snn.louvianmlr_0.5"
)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 2,
                                 resolution  = 1,
                                 cluster.name = "harmony.snn.louvianmlr_1"
)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 2,
                                 resolution  = 1.5,
                                 cluster.name = "harmony.snn.louvianmlr_1.5"
)

# SLM
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 3,
                                 resolution  = 0.5,
                                 cluster.name = "harmony.snn.slm_0.5"
)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 3,
                                 resolution  = 1,
                                 cluster.name = "harmony.snn.slm_1"
)
seu_parsef_clean <- FindClusters(seu_parsef_clean,
                                 algorithm   = 3,
                                 resolution  = 1.5,
                                 cluster.name = "harmony.snn.slm_1.5"
)

# ---------------------------- #
#   Harmony Integration Plots  #
# ---------------------------- #

p1 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.harmony",
  group.by    = "predicted.celltype.l1",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("Azimuth L1")

p2 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.harmony",
  group.by    = "predicted.celltype.l2",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("Azimuth L2")

p3 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.harmony",
  group.by    = "predicted.celltype.l3",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("Azimuth L3")

p4 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.harmony",
  group.by    = "harmony.snn.louvianmlr_0.5",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("harmony.snn.louvianmlr_0.5")

p5 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.harmony",
  group.by    = "harmony.snn.louvianmlr_1",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("harmony.snn.louvianmlr_1")

p6 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.harmony",
  group.by    = "harmony.snn.louvianmlr_1.5",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("harmony.snn.louvianmlr_1.5")

p7 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.harmony",
  group.by    = "harmony.snn.slm_0.5",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("harmony.snn.slm_0.5")

p8 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.harmony",
  group.by    = "harmony.snn.slm_1",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("harmony.snn.slm_1")

p9 <- DimPlot2(
  seu_parsef_clean,
  reduction   = "umap.harmony",
  group.by    = "harmony.snn.slm_1.5",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 1
) + ggtitle("harmony.snn.slm_1.5")

combined_plot_harmony <- wrap_plots(
  p1, p2, p3,
  p4, p5, p6,
  p7, p8, p9,
  ncol = 3, nrow = 3
)

ggsave(
  filename = file.path(integration_dir, "Harmony_Clustering_Grid.png"),
  plot     = combined_plot_harmony,
  width    = 24,
  height   = 17,
  dpi      = 300
)

# ----------------------------
# Save integrated object (qs2)
# ----------------------------
qs_save(
  seu_parsef_clean,
  file.path(base_dir, "saved_R_data", "BCG_PARSE_integrated_allMethods.qs2")
)
