# ================================================================
# P1078 Infant (ex vivo, 44 wks)
# Annotated UMAPs — IGRA+ vs IGRA- INFANT comparison
#   - Fine annotation, split by infant IGRA status
#   - Main (collapsed) annotation, split by infant IGRA status
#   - UMAP coloured by infant IGRA status (overlay)
# Plotting style matched to existing Annotated_UMAPs (DimPlot2)
# ================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(qs2)
library(SeuratExtend)

# ----------------------------------------------------------------
# 0) Paths & parameters
# ----------------------------------------------------------------
base_dir  <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/10x Infant (P1078) ex vivo 44wks"
load.path <- file.path(base_dir, "saved_R_data")

# Annotated, cleaned object. Adjust filename here if this dataset's
# object is named differently.
obj_file  <- file.path(load.path, "TB_ALL_Annotated_Clean.qs2")

# Reduction to plot on:
#   WNN / CITE-seq        -> "wnn.umap"   (current default, matches prior plots)
#   RNA-only 10x run      -> "umap"  (or "umap.unintegrated")
reduction <- "wnn.umap"

# Output directory (as requested)
annot_umap_dir <- file.path(base_dir, "Pre_Annotation_Plots", "Annotated_UMAPs")
dir.create(annot_umap_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------------
# 1) Load annotated object
# ----------------------------------------------------------------
TB_ALL <- qs_read(obj_file)
DefaultAssay(TB_ALL) <- "RNA"

# ----------------------------------------------------------------
# 2) IGRA factor + subset to infants with a known IGRA call
#    (IGRA+ kept first so it is the "positive" reference)
# ----------------------------------------------------------------
TB_ALL$IGRA_infant <- factor(TB_ALL$IGRA_infant, levels = c("IGRA+", "IGRA-"))

cat("\n--- IGRA_infant (all cells) ---\n")
print(table(TB_ALL$IGRA_infant, useNA = "ifany"))

TB_inf <- subset(TB_ALL, subset = !is.na(IGRA_infant))
TB_inf$IGRA_infant   <- droplevels(TB_inf$IGRA_infant)
TB_inf$celltype_fine <- droplevels(TB_inf$celltype_fine)
TB_inf$celltype_main <- droplevels(TB_inf$celltype_main)

cat("\n--- IGRA_infant (cells used for plotting) ---\n")
print(table(TB_inf$IGRA_infant))

igra_cols <- c("IGRA+" = "#F54927", "IGRA-" = "#118080")

# ----------------------------------------------------------------
# 3) Fine annotation, SPLIT by infant IGRA status
#    DimPlot2's split.by keeps the cell-type colour mapping
#    consistent across both panels, so they are directly comparable.
# ----------------------------------------------------------------
p_fine_split <- DimPlot2(
  TB_inf,
  reduction  = reduction,
  group.by   = "celltype_fine",
  split.by   = "IGRA_infant",
  cols       = "light",
  label      = TRUE,
  box        = TRUE,
  repel      = TRUE,
  label.size = 3,
  pt.size    = 0.4,
  raster     = FALSE,
  theme      = list(NoLegend(), NoAxes(), theme_umap_arrows())
)
ggsave(
  file.path(annot_umap_dir, "UMAP_celltype_fine_by_IGRA_infant.png"),
  p_fine_split, dpi = 300, width = 18, height = 9, bg = "white"
)

# ----------------------------------------------------------------
# 4) Main (collapsed) annotation, SPLIT by infant IGRA status
# ----------------------------------------------------------------
p_main_split <- DimPlot2(
  TB_inf,
  reduction  = reduction,
  group.by   = "celltype_main",
  split.by   = "IGRA_infant",
  cols       = "light",
  label      = TRUE,
  box        = TRUE,
  repel      = TRUE,
  label.size = 4,
  pt.size    = 0.4,
  raster     = FALSE,
  theme      = list(NoAxes(), theme_umap_arrows())
)
ggsave(
  file.path(annot_umap_dir, "UMAP_celltype_main_by_IGRA_infant.png"),
  p_main_split, dpi = 300, width = 20, height = 9, bg = "white"
)

# ----------------------------------------------------------------
# 5) UMAP coloured by infant IGRA status (overlay, all infant cells)
#    Shows where IGRA+ vs IGRA- cells sit on the shared embedding.
# ----------------------------------------------------------------
p_igra <- DimPlot2(
  TB_inf,
  reduction  = reduction,
  group.by   = "IGRA_infant",
  cols       = igra_cols,
  label      = FALSE,
  pt.size    = 0.4,
  raster     = FALSE,
  theme      = list(NoAxes(), theme_umap_arrows())
) + ggtitle("P1078 infant ex vivo — infant IGRA status")
ggsave(
  file.path(annot_umap_dir, "UMAP_IGRA_infant_overlay.png"),
  p_igra, dpi = 300, width = 11, height = 8, bg = "white"
)

message("\n=== Annotated IGRA_infant UMAPs complete ===")
message("Saved under: ", annot_umap_dir)

