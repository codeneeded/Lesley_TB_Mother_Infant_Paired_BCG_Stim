###############################################################################
#  TB CITE-seq Pipeline — Script 3: Cell Cycle Effect & Doublet Detection
#  Adapted from OPIS CITE-seq pipeline
###############################################################################

# ------------ #
# Libraries    #
# ------------ #
library(Seurat)
library(ggplot2)
library(scDblFinder)
library(Azimuth)
library(scCustomize)
library(patchwork)
library(qs2)

# ---------------------------- #
# Paths & Output Directories   #
# ---------------------------- #

base_dir       <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/TB_CITEseq"
cell_cycle_dir <- file.path(base_dir, "QC", "Cell_Cycle")
doublet_dir    <- file.path(base_dir, "QC", "Doublets")
saved_dir      <- file.path(base_dir, "saved_R_data")

dir.create(cell_cycle_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(doublet_dir,    recursive = TRUE, showWarnings = FALSE)

# ---------------------------- #
# Load Filtered Seurat Object  #
# ---------------------------- #

filtered_seurat <- qs_read(file.path(saved_dir, "TB_Seuratv5_filtered_seurat.qs2"))
seurat_clean    <- filtered_seurat

# ---------------------------- #
# Cell Cycle Scoring           #
# ---------------------------- #
DefaultAssay(seurat_clean) <- "RNA"

# Azimuth annotation using PBMC reference
seurat_phase <- RunAzimuth(seurat_clean, reference = "pbmcref")

# Basic Transformations
seurat_phase <- NormalizeData(seurat_phase)
seurat_phase <- FindVariableFeatures(seurat_phase, selection.method = "vst")
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(
  seurat_phase,
  features        = VariableFeatures(seurat_phase),
  ndims.print     = 6:10,
  nfeatures.print = 10
)
seurat_phase <- FindNeighbors(seurat_phase, dims = 1:30, reduction = "pca")
seurat_phase <- FindClusters(seurat_phase, resolution = 2, cluster.name = "unintegrated_clusters")
seurat_phase <- RunUMAP(
  seurat_phase,
  dims           = 1:30,
  reduction      = "pca",
  reduction.name = "umap.unintegrated"
)

# Get Seurat's built-in cell cycle gene lists
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

# Ridge plots for key markers
png(file.path(cell_cycle_dir, "RidgePlot_CellCycleMarkers.png"),
    width = 1800, height = 1200, res = 150)
RidgePlot(seurat_phase, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
dev.off()

# Violin plot of S.Score & G2M.Score by Sample_ID
png(file.path(cell_cycle_dir, "CellCycle_Scores_bySample.png"),
    width = 1800, height = 1200, res = 150)
VlnPlot(
  seurat_phase,
  features = c("S.Score", "G2M.Score"),
  group.by = "Sample_ID",
  pt.size  = 0.1
)
dev.off()

# Violin plot of S.Score & G2M.Score by Infant IGRA status
png(file.path(cell_cycle_dir, "CellCycle_Scores_byIGRA_infant.png"),
    width = 1400, height = 1000, res = 150)
VlnPlot(
  seurat_phase,
  features = c("S.Score", "G2M.Score"),
  group.by = "IGRA_infant",
  pt.size  = 0.1
)
dev.off()

# Violin plot of S.Score & G2M.Score by Maternal IGRA status
png(file.path(cell_cycle_dir, "CellCycle_Scores_byIGRA_mother.png"),
    width = 1400, height = 1000, res = 150)
VlnPlot(
  seurat_phase,
  features = c("S.Score", "G2M.Score"),
  group.by = "IGRA_mother",
  pt.size  = 0.1
)
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
  pt.size    = 1
)
ggsave(
  filename = file.path(cell_cycle_dir, "PCA_by_CellCyclePhase.png"),
  plot     = p1,
  width    = 13, height = 8, dpi = 300
)

# Azimuth annotated over UMAP (unintegrated)
p2 <- DimPlot_scCustom(
  seurat_phase,
  reduction  = "umap.unintegrated",
  group.by   = "predicted.celltype.l2",
  label      = TRUE,
  repel      = TRUE,
  label.box  = TRUE,
  label.size = 3.5,
  pt.size    = 1
)
ggsave(
  filename = file.path(cell_cycle_dir, "Azimuth_Umap_Unintegrated.png"),
  plot     = p2,
  width    = 13, height = 8, dpi = 300
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
  pt.size    = 1
)
ggsave(
  filename = file.path(cell_cycle_dir, "Cell_Cycle_Umap_Unintegrated.png"),
  plot     = p3,
  width    = 13, height = 8, dpi = 300
)

# Phase UMAP split by Infant IGRA status
p4 <- DimPlot_scCustom(
  seurat_phase,
  reduction  = "umap.unintegrated",
  group.by   = "Phase",
  split.by   = "IGRA_infant",
  label      = TRUE,
  repel      = TRUE,
  label.box  = TRUE,
  label.size = 3.5,
  pt.size    = 1
)
ggsave(
  filename = file.path(cell_cycle_dir, "Cell_Cycle_IGRA_infant_Umap_Unintegrated.png"),
  plot     = p4,
  width    = 16, height = 8, dpi = 300
)

# Phase UMAP split by Maternal IGRA status
p5 <- DimPlot_scCustom(
  seurat_phase,
  reduction  = "umap.unintegrated",
  group.by   = "Phase",
  split.by   = "IGRA_mother",
  label      = TRUE,
  repel      = TRUE,
  label.box  = TRUE,
  label.size = 3.5,
  pt.size    = 1
)
ggsave(
  filename = file.path(cell_cycle_dir, "Cell_Cycle_IGRA_mother_Umap_Unintegrated.png"),
  plot     = p5,
  width    = 16, height = 8, dpi = 300
)

# ---------------------------- #
# Doublet Detection (scDblFinder)
# ---------------------------- #

# Split by Sample_ID for per-sample doublet calling
split_seurat <- SplitObject(seurat_phase, split.by = "Sample_ID")
samples      <- names(split_seurat)

for (i in samples) {
  message("Running scDblFinder for sample: ", i)
  
  counts_mat <- GetAssayData(split_seurat[[i]], slot = "counts")
  sce        <- scDblFinder(counts_mat)
  
  # Store doublet info
  split_seurat[[i]]$scDblFinder.score <- sce$scDblFinder.score
  split_seurat[[i]]$scDblFinder.class <- sce$scDblFinder.class
  
  # UMAP by doublet classification
  p_doublet <- DimPlot_scCustom(
    split_seurat[[i]],
    reduction  = "umap.unintegrated",
    group.by   = "scDblFinder.class",
    colors_use = c("#1f78b4", "#e31a1c"),
    label      = TRUE,
    repel      = TRUE,
    label.box  = TRUE,
    label.size = 3.5,
    pt.size    = 1
  ) + ggtitle(paste0("Sample ", i, " — Doublets"))
  
  ggsave(
    filename = file.path(doublet_dir, paste0(i, "_Doublets.png")),
    plot     = p_doublet,
    width    = 13, height = 8, dpi = 300
  )
}

# --- Doublet summary table ---
doublet_summary <- do.call(rbind, lapply(samples, function(i) {
  tbl <- table(split_seurat[[i]]$scDblFinder.class)
  data.frame(
    Sample_ID = i,
    singlet   = as.integer(tbl["singlet"]),
    doublet   = as.integer(tbl["doublet"]),
    pct_doublet = round(100 * tbl["doublet"] / sum(tbl), 2),
    stringsAsFactors = FALSE
  )
}))
write.csv(doublet_summary,
          file = file.path(doublet_dir, "Doublet_Summary.csv"),
          row.names = FALSE)
message("Doublet summary:\n")
print(doublet_summary)

# ---------------------------- #
# Remove Doublets              #
# ---------------------------- #
for (i in samples) {
  split_seurat[[i]] <- subset(split_seurat[[i]], subset = scDblFinder.class == "singlet")
}

# Merge back post-doublet removal
seurat_phase_clean <- merge(split_seurat[[1]], y = split_seurat[-1])

message("Cells before doublet removal: ", ncol(seurat_phase))
message("Cells after doublet removal:  ", ncol(seurat_phase_clean))

# ---------------------------- #
# Save Output                  #
# ---------------------------- #
qs_save(
  seurat_phase_clean,
  file = file.path(saved_dir, "TB_Seuratv5_filtered_CellCycle_DoubletClean.qs2")
)

message("Done! Saved: TB_Seuratv5_filtered_CellCycle_DoubletClean.qs2")

