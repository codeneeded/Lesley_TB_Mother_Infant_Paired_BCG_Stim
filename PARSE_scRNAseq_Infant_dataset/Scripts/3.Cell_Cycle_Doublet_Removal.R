# ========================================= #
#   Cell Cycle Effect & Doublet Detection   #
#   BCG_Stim_Infants_12_44_wks_PARSE_GEX    #
# ========================================= #

# ------------ #
# Libraries    #
# ------------ #
library(Seurat)
library(ggplot2)
library(scCustomize)
library(scDblFinder)
library(SingleCellExperiment)
library(Azimuth)
library(qs2)

# ---------------------------- #
# Paths & Output Dirs          #
# ---------------------------- #
base_dir       <- "C:/Users/ammas/Documents/BCG_Stim_Infants_12_44_wks_PARSE_GEX"
cell_cycle_dir <- file.path(base_dir, "QC_Plots", "Cell_Cycle")
doublet_dir    <- file.path(base_dir, "QC_Plots", "Doublets")

dir.create(cell_cycle_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(doublet_dir,    recursive = TRUE, showWarnings = FALSE)

# ---------------------------- #
# Load Filtered Seurat Object  #
# (post-QC, before cell cycle) #
# ---------------------------- #

# This is the object you saved with qs2:
# qs_save(filtered_final, "BCG_PARSE_filtered_postQC.qs2")
seurat_clean <- qs_read(
  file.path(base_dir, "saved_R_data", "BCG_PARSE_filtered_postQC.qs2")
)


# ---------------------------- #
# Cell Cycle Scoring           #
# ---------------------------- #
DefaultAssay(seurat_clean) <- "RNA"

# Optional: Azimuth PBMC reference mapping (gives predicted.celltype.l2)
# This step needs the Azimuth PBMC reference installed / accessible.
seurat_phase <- RunAzimuth(seurat_clean, reference = "pbmcref")

# Basic Seurat pipeline on the unintegrated object
seurat_phase <- NormalizeData(seurat_phase)
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst")
seurat_phase <- ScaleData(seurat_phase)

seurat_phase <- RunPCA(
  seurat_phase,
  features       = VariableFeatures(seurat_phase),
  ndims.print    = 6:10,
  nfeatures.print = 10
)

seurat_phase <- FindNeighbors(seurat_phase, dims = 1:30, reduction = "pca")
seurat_phase <- FindClusters(
  seurat_phase,
  resolution   = 2,
  cluster.name = "unintegrated_clusters"
)

seurat_phase <- RunUMAP(
  seurat_phase,
  dims           = 1:30,
  reduction      = "pca",
  reduction.name = "umap.unintegrated"
)

# Seurat’s built-in cell cycle gene lists
s.genes   <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Score cells for cell cycle phase
seurat_phase <- CellCycleScoring(
  seurat_phase,
  s.features   = s.genes,
  g2m.features = g2m.genes,
  set.ident    = TRUE
)

# ---------------------------- #
# Visualize Cell Cycle Effect  #
# ---------------------------- #

# Ridge plots for key cell-cycle markers
png(file.path(cell_cycle_dir, "RidgePlot_CellCycleMarkers.png"),
    width = 1800, height = 1200)
RidgePlot(
  seurat_phase,
  features = c("PCNA", "TOP2A", "MCM6", "MKI67"),
  ncol     = 2
)
dev.off()

# Violin plot of S.Score & G2M.Score by sample_id
png(file.path(cell_cycle_dir, "CellCycle_Scores_bySampleID.png"),
    width = 2000, height = 1400)
VlnPlot(
  seurat_phase,
  features = c("S.Score", "G2M.Score"),
  group.by = "sample_id",
  pt.size  = 0.1,
  ncol     = 2
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# PCA colored by Phase
p1 <- DimPlot_scCustom(
  seurat_phase,
  reduction  = "pca",
  group.by   = "Phase",
  label      = TRUE,
  repel      = TRUE,
  label.box  = TRUE,
  label.size = 3.5,
  pt.size    = 0.8
)
ggsave(
  filename = file.path(cell_cycle_dir, "PCA_by_CellCyclePhase.png"),
  plot     = p1,
  width    = 13,
  height   = 8,
  dpi      = 300
)

# Azimuth annotation over UMAP (unintegrated)
p2 <- DimPlot_scCustom(
  seurat_phase,
  reduction  = "umap.unintegrated",
  group.by   = "predicted.celltype.l2",
  label      = TRUE,
  repel      = TRUE,
  label.box  = TRUE,
  label.size = 3.5,
  pt.size    = 0.8
)
ggsave(
  filename = file.path(cell_cycle_dir, "Azimuth_Umap_Unintegrated.png"),
  plot     = p2,
  width    = 13,
  height   = 8,
  dpi      = 300
)

# Phase UMAP (unintegrated)
p3 <- DimPlot_scCustom(
  seurat_phase,
  reduction  = "umap.unintegrated",
  group.by   = "Phase",
  label      = TRUE,
  repel      = TRUE,
  label.box  = TRUE,
  label.size = 3.5,
  pt.size    = 0.8
)
ggsave(
  filename = file.path(cell_cycle_dir, "Cell_Cycle_Umap_Unintegrated.png"),
  plot     = p3,
  width    = 13,
  height   = 8,
  dpi      = 300
)

# Phase UMAP split by stimulation (BCG_treated vs media_treated)
p4 <- DimPlot_scCustom(
  seurat_phase,
  reduction  = "umap.unintegrated",
  group.by   = "Phase",
  split.by   = "stim",
  label      = TRUE,
  repel      = TRUE,
  label.box  = TRUE,
  label.size = 3.0,
  pt.size    = 0.7
)
ggsave(
  filename = file.path(cell_cycle_dir, "Cell_Cycle_Umap_Unintegrated_byStim.png"),
  plot     = p4,
  width    = 13,
  height   = 8,
  dpi      = 300
)

# ---------------------------- #
# Doublet Detection (scDblFinder)
# ---------------------------- #

# We do doublets per library (orig.ident), since each is a separate well:
# 6044931_week12_BCG_treated, 6044931_week12_media_treated, etc.
split_seurat <- SplitObject(seurat_phase, split.by = "orig.ident")
samples      <- names(split_seurat)

for (i in samples) {
  message("Running scDblFinder on: ", i)
  
  # Convert each Seurat object to SingleCellExperiment
  sce <- as.SingleCellExperiment(split_seurat[[i]])
  
  # Run scDblFinder
  sce <- scDblFinder(sce)
  
  # Transfer doublet metadata back into Seurat object
  # (scDblFinder preserves column order)
  split_seurat[[i]]$scDblFinder.score <- sce$scDblFinder.score
  split_seurat[[i]]$scDblFinder.class <- sce$scDblFinder.class
  
  # Plot doublet vs singlet on UMAP
  p_doublet <- DimPlot_scCustom(
    split_seurat[[i]],
    reduction  = "umap.unintegrated",
    group.by   = "scDblFinder.class",
    colors_use = c("singlet" = "#1f78b4", "doublet" = "#e31a1c"),
    label      = FALSE,
    pt.size    = 0.9
  ) + ggtitle(paste0(i, " — scDblFinder Doublets"))
  
  ggsave(
    filename = file.path(doublet_dir, paste0(i, "_Doublets_scDblFinder.png")),
    plot     = p_doublet,
    width    = 10,
    height   = 7,
    dpi      = 300
  )
}

# ---------------------------- #
# Remove Doublets              #
# ---------------------------- #
for (i in samples) {
  split_seurat[[i]] <- subset(
    split_seurat[[i]],
    subset = scDblFinder.class == "singlet"
  )
}

# Merge back post-doublet removal
seurat_phase_clean <- merge(split_seurat[[1]], y = split_seurat[-1])

# ---------------------------- #
# Save Output (qs2)            #
# ---------------------------- #
qs_save(
  seurat_phase_clean,
  file.path(base_dir, "saved_R_data", "BCG_PARSE_postQC_CellCycle_DoubletClean.qs2")
)
