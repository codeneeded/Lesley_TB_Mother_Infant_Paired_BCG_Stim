###############################################
# BCG Parse — Cluster Plots, Vln Plots, DGE
# Using FastMNN (umap.mnn.rna) and mnn.snn.louvianmlr_1
###############################################

library(Seurat)
library(SeuratExtend)   # DimPlot2, VlnPlot2
library(ggplot2)
library(dplyr)
library(EnhancedVolcano) # optional if you want volcano later
library(qs2)             # qs2 package (no qs:: prefix, use qsread/qssave if needed)

# ---------------------------- #
# Paths
# ---------------------------- #
base_dir   <- "/home/akshay-iyer/Documents/BCG_Stim_Infants_12_44_wks_PARSE_GEX"
saved_dir  <- file.path(base_dir, "saved_R_data")

plots_dir        <- file.path(base_dir, "Cluster_Vln_Plots")
cluster_umap_dir <- file.path(plots_dir, "UMAPs")
vln_dir_qc       <- file.path(plots_dir, "Vln_QC")
vln_dir_genes    <- file.path(plots_dir, "Vln_Genes")

dir.create(plots_dir,        showWarnings = FALSE, recursive = TRUE)
dir.create(cluster_umap_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(vln_dir_qc,       showWarnings = FALSE, recursive = TRUE)
dir.create(vln_dir_genes,    showWarnings = FALSE, recursive = TRUE)

dge_dir <- file.path(base_dir, "DGE")
dir.create(dge_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------- #
# Load Integrated Object
# ---------------------------- #
# From previous step: BCG_PARSE_integrated_allMethods.qs2
seu <- qs_read(file.path(saved_dir, "BCG_PARSE_integrated_allMethods.qs2"))

seu$mnn.snn.louvianmlr_1 <- factor(
  seu$mnn.snn.louvianmlr_1,
  levels = sort(unique(seu$mnn.snn.louvianmlr_1))
)

# We’re using FastMNN + mnn.snn.louvianmlr_1
Idents(seu) <- "mnn.snn.louvianmlr_1"

# Optional: drop very tiny clusters (<100 cells)
cluster_sizes <- table(Idents(seu))
large_clusters <- names(cluster_sizes[cluster_sizes >= 100])
seu <- subset(seu, idents = large_clusters)

# Sanity: make sure metadata is factorized cleanly
seu$stim      <- factor(seu$stim)       # "BCG_treated", "media_treated"
seu$timepoint <- factor(seu$timepoint)  # "week12", "week44"
########## Fix cluster ordering

old <- levels(seu$mnn.snn.louvianmlr_1)
new <- seq_along(old) - 1  # new labels 0..N

# rename identities
seu <- RenameIdents(seu, setNames(new, old))

# store back into metadata
seu$mnn.snn.louvianmlr_1 <- Idents(seu)

# reorder again to be safe
seu$mnn.snn.louvianmlr_1 <- factor(
  seu$mnn.snn.louvianmlr_1,
  levels = sort(unique(seu$mnn.snn.louvianmlr_1))
)


########
###############################################
# 1) Cluster UMAP Plots (FastMNN)
###############################################

# UMAP colored by cluster (FastMNN)
p_umap_cluster <- DimPlot2(
  seu,
  reduction   = "umap.mnn.rna",
  group.by    = "mnn.snn.louvianmlr_1",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 0.6
) + ggtitle("FastMNN — Clusters (mnn.snn.louvianmlr_1)")

ggsave(
  filename = file.path(cluster_umap_dir, "UMAP_mnn_clusters.png"),
  plot     = p_umap_cluster,
  width    = 10,
  height   = 8,
  dpi      = 300,
  bg='white'
)

# UMAP split by stim
p_umap_stim <- DimPlot2(
  seu,
  reduction   = "umap.mnn.rna",
  group.by    = "mnn.snn.louvianmlr_1",
  split.by    = "stim",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 0.6
) + ggtitle("FastMNN — Clusters split by stim")

ggsave(
  filename = file.path(cluster_umap_dir, "UMAP_mnn_clusters_byStim.png"),
  plot     = p_umap_stim,
  width    = 12,
  height   = 8,
  dpi      = 300,
  bg='white'
)

# UMAP split by timepoint
p_umap_tp <- DimPlot2(
  seu,
  reduction   = "umap.mnn.rna",
  group.by    = "mnn.snn.louvianmlr_1",
  split.by    = "timepoint",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 0.6
) + ggtitle("FastMNN — Clusters split by timepoint")

ggsave(
  filename = file.path(cluster_umap_dir, "UMAP_mnn_clusters_byTimepoint.png"),
  plot     = p_umap_tp,
  width    = 12,
  height   = 8,
  dpi      = 300,
  bg='white'
)

seu$stim_time <- paste0(seu$timepoint, "_", seu$stim)
seu$stim_time <- factor(seu$stim_time)

p_umap_st_tp <- DimPlot2(
  seu,
  reduction   = "umap.mnn.rna",
  group.by    = "mnn.snn.louvianmlr_1",
  split.by    = "stim_time",
  label       = TRUE,
  box         = TRUE,
  repel       = TRUE,
  label.color = "black",
  pt.size     = 0.6
) + ggtitle("FastMNN — Clusters split by stim × timepoint")

ggsave(
  filename = file.path(cluster_umap_dir, "UMAP_mnn_clusters_byStim_Timepoint.png"),
  plot     = p_umap_st_tp,
  width    = 16,
  height   = 10,
  dpi      = 300,
  bg='white'
)

###############################################
# 2) Violin Plots
###############################################

# ---------- QC metrics by cluster ----------
qc_feats <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_hb", "percent_plat")

p_vln_qc <- VlnPlot2(
  seu,
  features  = qc_feats,
  group.by  = "mnn.snn.louvianmlr_1",
  pt.size   = 0.1,
  ncol      = 3,
  show.mean = TRUE
) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1)
)

ggsave(
  filename = file.path(vln_dir_qc, "Vln_QC_byCluster.png"),
  plot     = p_vln_qc,
  width    = 16,
  height   = 10,
  dpi      = 300,
  bg='white'
)

# ---------- Gene expression by cluster ----------
# reuse your RNA marker panel if you like, or trim this list
rna.features <- c(
  'ASCL2','BATF','BATF3','BCL6','C1QBP','CCL2','CCL3','CCL4L2','CCL5','CCND3','CD14','CD19','CD1C',
  'CD200','CD27','CD3D','CD3E','CD36','CD4','CD40','CD40LG','CD70','CD7','CD79A','CD8A','CD8B',
  'CLEC9A','CR2','CTLA4','CTSW','CXCL8','CXCR3','CXCR5','EBI3','ENTPD1','FABP5','FCGR2B','FCGR3A',
  'FCRL5','FOXP3','FASLG','GNLY','GP1BA','GP9','GATA3','GZMK','HAVCR2','HIF1A','HIST1H4C','HLA-DPA1',
  'HLA-DRA','HLA-DRB1','ICOS','IFI30','IFNG','IGFBP2','IGFBP4','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHM','IKZF2','IL10',
  'IL17A','IL18BP','IL18RAP','IL1B','IL21','IL2RA','IL2RB','IRF4','IRF8','ITGAX','JCHAIN','KLRB1',
  'KLRD1','KLRC1','LAG3','LDHA','LGALS1','LTA','LTB','MAF','MAL','MALAT1','MIR155HG','MKI67',
  'MT-ND1','MT-ND5','MS4A1','NELL2','NCAM1','NKG7','NR4A1','PDCD1','PF4','PPBP','PRDM1','PRF1',
  'RORC','SELL','SERPINA1','SERPING1','SH2D1A','TCF4','TCF7','TIGIT','TNF','TNFAIP2','TNFRSF18',
  'TNFRSF4','TNFRSF9','TOX','TBX21','TRBC1','TRDC','TRDV1','TRDV2','TRGC1','TRGC2','TRGV9','XBP1',
  'XCL1','XCL2','ZBTB16','ZEB2'
)

for (gene in rna.features) {
  if (!gene %in% rownames(seu)) next
  
  vln_file <- file.path(vln_dir_genes, paste0("Vln_", gene, "_byCluster.png"))
  if (file.exists(vln_file)) next
  
  p_vln_gene <- VlnPlot2(
    seu,
    features   = gene,
    group.by   = "mnn.snn.louvianmlr_1",
    show.mean  = TRUE,
    mean_colors = c("red", "blue")
  ) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + ggtitle(paste0("Violin: ", gene, " by cluster"))
  
  ggsave(
    filename = vln_file,
    plot     = p_vln_gene,
    width    = 10,
    height   = 8,
    dpi      = 300,
    bg='white'
  )
}

###############################################
# 3) DGE — stim vs unstim & timepoint effects
###############################################
DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)

# We’ll do DGE per cluster using MAST,
# controlling for nCount_RNA as a latent variable.

cluster_levels <- levels(Idents(seu))
timepoints     <- levels(seu$timepoint)     # "week12", "week44"
stims          <- levels(seu$stim)         # "BCG_treated", "media_treated"

# Subfolders for DGE outputs
dge_stim_vs_unstim_dir <- file.path(dge_dir, "Stim_vs_Unstim_byTimepoint")
dge_tp_dir             <- file.path(dge_dir, "Week44_vs_Week12_byStim")

dir.create(dge_stim_vs_unstim_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dge_tp_dir,             showWarnings = FALSE, recursive = TRUE)

min_cells_per_group <- 50

# ---------- 3A) Stim vs Unstim at each timepoint (within cluster) ----------
for (cl in cluster_levels) {
  for (tp in timepoints) {
    
    message("DGE: Cluster ", cl, " | Timepoint ", tp, " | BCG_treated vs media_treated")
    
    # subset: this cluster + this timepoint
    seu_sub <- subset(
      seu,
      subset = (mnn.snn.louvianmlr_1 == cl & timepoint == tp)
    )
    
    # total cells check
    if (ncol(seu_sub) == 0) {
      message("  -> Skipping: no cells in this subset")
      next
    }
    
    # need both stim groups present
    if (length(unique(seu_sub$stim)) < 2) {
      message("  -> Skipping: Only one stim level present")
      next
    }
    
    # check per-group counts
    tab_stim <- table(seu_sub$stim)
    if (any(tab_stim < min_cells_per_group)) {
      message("  -> Skipping: too few cells per stim group: ",
              paste(names(tab_stim), tab_stim, collapse = ", "))
      next
    }
    
    seu_sub$stim <- droplevels(factor(seu_sub$stim))
    
    # DGE: BCG_treated vs media_treated
    markers <- FindMarkers(
      seu_sub,
      ident.1         = "BCG_treated",
      ident.2         = "media_treated",
      group.by        = "stim",
      test.use        = "MAST",
      latent.vars     = c("nCount_RNA", "sample_id"),  # also paired across stim
      logfc.threshold = 0.1,
      min.pct         = 0.1
    )
    
    out_file <- file.path(
      dge_stim_vs_unstim_dir,
      paste0("DGE_cluster", cl, "_", tp, "_BCG_vs_media.csv")
    )
    write.csv(markers, out_file)
  }
}

# ---------- 3B) Week44 vs Week12 within each stim (per cluster) ----------
for (cl in cluster_levels) {
  for (st in stims) {
    
    message("DGE: Cluster ", cl, " | Stim ", st, " | week44 vs week12")
    
    seu_sub <- subset(
      seu,
      subset = (mnn.snn.louvianmlr_1 == cl & stim == st)
    )
    
    # total cells check
    if (ncol(seu_sub) == 0) {
      message("  -> Skipping: no cells in this subset")
      next
    }
    
    # need both timepoints present
    if (length(unique(seu_sub$timepoint)) < 2) {
      message("  -> Skipping: Only one timepoint present")
      next
    }
    
    tab_tp <- table(seu_sub$timepoint)
    if (any(tab_tp < min_cells_per_group)) {
      message("  -> Skipping: too few cells per timepoint: ",
              paste(names(tab_tp), tab_tp, collapse = ", "))
      next
    }
    
    seu_sub$timepoint <- droplevels(factor(seu_sub$timepoint))
    
    # DGE: week44 vs week12 (paired via sample_id)
    markers <- FindMarkers(
      seu_sub,
      ident.1         = "week44",
      ident.2         = "week12",
      group.by        = "timepoint",
      test.use        = "MAST",
      latent.vars     = c("nCount_RNA", "sample_id"),
      logfc.threshold = 0.1,
      min.pct         = 0.1
    )
    
    out_file <- file.path(
      dge_tp_dir,
      paste0("DGE_cluster", cl, "_", st, "_week44_vs_week12.csv")
    )
    write.csv(markers, out_file)
  }
}

# Save joined-layer object for all downstream analyses
qs_save(
  seu,
  file = file.path(base_dir, "saved_R_data", "seu_parsef_integrated_joined.qs2")
)

