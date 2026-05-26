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
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(stringr)
library(purrr)
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

# ---------- 3C) IGRA+vs- ----------
igra_map <- data.frame(
  sample_id = c("6044931","6063411","6063421","6063481","6061161","6062151","6065761","6067021","6071261","6044981"),
  Maternal_IGRA = c("Negative","Negative","Positive","Positive","Negative","Positive","Positive","Negative","Positive","Negative"),
  IGRA_status   = c("Negative","Negative","Negative","Negative","Positive","Positive","Positive","Positive","Positive","Positive"),
  stringsAsFactors = FALSE
)

idx <- match(as.character(seu$sample_id), igra_map$sample_id)

seu$Maternal_IGRA <- igra_map$Maternal_IGRA[idx]
seu$IGRA_status   <- igra_map$IGRA_status[idx]

# Edit these if your object uses different names
cluster_col <- "mnn.snn.louvianmlr_1"
igra_col    <- "IGRA_status"   # should contain IGRA+ / IGRA-
stim_col    <- "stim"          # "BCG_treated" / "media_treated"
tp_col      <- "timepoint"     # "week12" / "week44"
sample_col  <- "sample_id"     # paired / donor ID
latent_vars <- c("nCount_RNA")

stopifnot(all(c(cluster_col, igra_col, stim_col, tp_col, sample_col) %in% colnames(seu@meta.data)))
# --------- Output dirs ----------
dge_igra_dir <- file.path(dge_dir, "IGRApos_vs_IGRAneg_byStimTimepoint")
dir.create(dge_igra_dir, recursive = TRUE, showWarnings = FALSE)

heat_dir <- file.path(base_dir, "Heatmaps", "IGRA_Stim_Timepoint_AvgExpr")
dir.create(heat_dir, recursive = TRUE, showWarnings = FALSE)

min_cells_per_group <- 50
padj_cutoff <- 0.05

# --------- Comparisons requested ----------
# (IGRA+ vs IGRA-) within each (stim, timepoint), per cluster
requested_comparisons <- tibble::tribble(
  ~stim_level,      ~tp_level, ~comp_label,
  "media_treated",  "week12",  "IGRApos_Med_12wk_vs_IGRAneg_Med_12wk",
  "BCG_treated",    "week12",  "IGRApos_BCG_12wk_vs_IGRAneg_BCG_12wk",
  "media_treated",  "week44",  "IGRApos_Med_44wk_vs_IGRAneg_Med_44wk",
  "BCG_treated",    "week44",  "IGRApos_BCG_44wk_vs_IGRAneg_BCG_44wk"
)

# --------- Helper: safe FindMarkers ----------
run_igra_dge_one <- function(seu, cl, stim_level, tp_level, comp_label) {
  
  seu_sub <- subset(
    seu,
    subset = (!!as.name(cluster_col) == cl &
                !!as.name(stim_col)    == stim_level &
                !!as.name(tp_col)      == tp_level)
  )
  
  if (ncol(seu_sub) == 0) {
    message("  -> Skipping (no cells): cluster=", cl, " ", comp_label)
    return(NULL)
  }
  
  # ensure both IGRA levels exist
  lv <- unique(as.character(seu_sub@meta.data[[igra_col]]))
  if (length(lv) < 2) {
    message("  -> Skipping (one IGRA level): cluster=", cl, " ", comp_label)
    return(NULL)
  }
  
  tab <- table(seu_sub@meta.data[[igra_col]])
  if (any(tab < min_cells_per_group)) {
    message("  -> Skipping (too few cells): cluster=", cl, " ", comp_label,
            " | ", paste(names(tab), tab, collapse = ", "))
    return(NULL)
  }
  
  # enforce factor levels
  seu_sub@meta.data[[igra_col]] <- droplevels(factor(seu_sub@meta.data[[igra_col]]))
  
  # pick which string corresponds to pos/neg
  # (robust to "IGRA+" vs "pos" naming)
  levels_igra <- levels(seu_sub@meta.data[[igra_col]])
  ig_pos <- levels_igra[grepl("\\+|pos", levels_igra, ignore.case = TRUE)][1]
  ig_neg <- levels_igra[grepl("\\-|neg", levels_igra, ignore.case = TRUE)][1]
  
  if (is.na(ig_pos) || is.na(ig_neg)) {
    message("  -> Skipping (could not detect IGRA+/IGRA- labels): cluster=", cl, " ", comp_label,
            " | levels: ", paste(levels_igra, collapse = ", "))
    return(NULL)
  }
  
  markers <- FindMarkers(
    seu_sub,
    ident.1         = ig_pos,
    ident.2         = ig_neg,
    group.by        = igra_col,
    test.use        = "MAST",
    latent.vars     = latent_vars,
    logfc.threshold = 0.1,
    min.pct         = 0.1
  )
  
  if (is.null(markers) || nrow(markers) == 0) return(NULL)
  
  markers <- markers %>%
    tibble::rownames_to_column("gene") %>%
    mutate(
      cluster = as.character(cl),
      stim    = stim_level,
      timepoint = tp_level,
      comparison = comp_label
    )
  
  # save full table
  out_file <- file.path(dge_igra_dir, paste0("DGE_cluster", cl, "__", comp_label, ".csv"))
  write.csv(markers, out_file, row.names = FALSE)
  
  # return for aggregation
  markers
}

# --------- Run all clusters x 4 comparisons ----------
cluster_levels <- levels(factor(seu@meta.data[[cluster_col]]))
all_markers_list <- list()

for (cl in cluster_levels) {
  message("=== Cluster ", cl, " ===")
  for (i in seq_len(nrow(requested_comparisons))) {
    rr <- requested_comparisons[i, ]
    res <- run_igra_dge_one(
      seu = seu,
      cl  = cl,
      stim_level = rr$stim_level,
      tp_level   = rr$tp_level,
      comp_label = rr$comp_label
    )
    if (!is.null(res)) all_markers_list[[length(all_markers_list) + 1]] <- res
  }
}

all_markers <- dplyr::bind_rows(all_markers_list)


########## Combined CSV

# 1) Find all existing DGE cluster files
dge_files <- list.files(
  path = dge_igra_dir,
  pattern = "^DGE_cluster.*__IGRApos_.*\\.csv$",
  full.names = TRUE
)

if (length(dge_files) == 0) {
  stop("No files found in dge_igra_dir matching: DGE_cluster*__IGRApos_*.csv")
}

# 2) Helper: parse cluster + comparison from filename (works even if columns missing)
parse_from_filename <- function(f) {
  bn <- basename(f)
  cluster <- str_match(bn, "^DGE_cluster([^_]+)__")[,2]
  comparison <- str_match(bn, "__(.+)\\.csv$")[,2]
  list(cluster = cluster, comparison = comparison)
}

# 3) Read + standardize each file, add cluster/comparison if not already present
all_markers <- purrr::map_dfr(dge_files, function(f) {
  x <- readr::read_csv(f, show_col_types = FALSE)
  
  meta <- parse_from_filename(f)
  
  # Standardize gene column if needed
  if (!"gene" %in% names(x)) {
    # Some Seurat outputs are rownames saved into X column sometimes
    # If you see "X" with gene symbols, use that:
    if ("X" %in% names(x)) {
      x <- x %>% dplyr::rename(gene = X)
    }
  }
  
  # Ensure cluster/comparison exist
  if (!"cluster" %in% names(x))     x$cluster <- meta$cluster
  if (!"comparison" %in% names(x))  x$comparison <- meta$comparison
  
  # If stim/timepoint not present, infer from comparison label (optional)
  if (!"stim" %in% names(x)) {
    x$stim <- dplyr::case_when(
      str_detect(x$comparison, "_BCG_") ~ "BCG_treated",
      str_detect(x$comparison, "_Med_") ~ "media_treated",
      TRUE ~ NA_character_
    )
  }
  if (!"timepoint" %in% names(x)) {
    x$timepoint <- dplyr::case_when(
      str_detect(x$comparison, "_12wk_") ~ "week12",
      str_detect(x$comparison, "_44wk_") ~ "week44",
      TRUE ~ NA_character_
    )
  }
  
  # track file source (handy for debugging)
  x$source_file <- basename(f)
  
  x
})

# 4) Save combined LONG CSV (ALL results, no filtering)
combined_out <- file.path(dge_igra_dir, "IGRApos_vs_IGRAneg__ALL_clusters__ALL_results_LONG.csv")
write.csv(all_markers, combined_out, row.names = FALSE)

############################################################
# 5) HEATMAPS (avg expression) for selected cluster groups
#    Columns = IGRA (pos/neg) x Stim (Med/BCG) x Timepoint
#    For each cluster set: one heatmap per cluster (saved)
############################################################
############################################################
# STEP A (HEAVY): Select genes for 3 group heatmaps/dotplots
# - Finds markers across 8 group8 identities within each cluster
# - Allows BOTH directions (only.pos = FALSE) and ranks by abs(avg_log2FC)
# - Produces exactly 100 genes per group (CD8_GD, CD4, Monocytes)
# - Saves gene list + gene-to-source_cluster mapping for dotplot Option B
############################################################


# ---- REQUIRED existing objects/vars ----
# seu, heat_dir, cluster_col, igra_col, stim_col, tp_col, latent_vars
stopifnot(all(c(cluster_col, igra_col, stim_col, tp_col) %in% colnames(seu@meta.data)))

# ---- cluster sets (3 heatmaps total) ----
cluster_sets <- list(
  CD8_GD     = c("5","10","14","28"),
  CD4        = c("3","7","11","24"),
  Monocytes  = c("4","8","15")
)

# ---- output dirs ----
heat_marker_dir <- file.path(heat_dir, "MarkerSelection")
dir.create(heat_marker_dir, recursive = TRUE, showWarnings = FALSE)

# ---- detect IGRA labels robustly ----
ig_levels_all <- levels(factor(seu@meta.data[[igra_col]]))
ig_pos_all <- ig_levels_all[grepl("\\+|pos", ig_levels_all, ignore.case = TRUE)][1]
ig_neg_all <- ig_levels_all[grepl("\\-|neg", ig_levels_all, ignore.case = TRUE)][1]
if (is.na(ig_pos_all) || is.na(ig_neg_all)) {
  stop("Could not auto-detect IGRA+/IGRA- labels in ", igra_col,
       ". Levels: ", paste(ig_levels_all, collapse = ", "))
}

# ---- fixed 8-column order ----
group_levels <- c(
  paste(ig_pos_all, "media_treated", "week12", sep="__"),
  paste(ig_neg_all, "media_treated", "week12", sep="__"),
  paste(ig_pos_all, "BCG_treated",   "week12", sep="__"),
  paste(ig_neg_all, "BCG_treated",   "week12", sep="__"),
  paste(ig_pos_all, "media_treated", "week44", sep="__"),
  paste(ig_neg_all, "media_treated", "week44", sep="__"),
  paste(ig_pos_all, "BCG_treated",   "week44", sep="__"),
  paste(ig_neg_all, "BCG_treated",   "week44", sep="__")
)

# ---- knobs ----
genes_per_heatmap <- 100
top_n_per_group8_per_cluster <- 6   # per cluster, per group8 identity; union is capped to 100 anyway

# ---- helper: markers across group8 within ONE cluster ----
markers_group8_one_cluster <- function(seu, cl) {
  
  # Robust subsetting by cells (avoids Seurat subset NSE issues)
  cells_keep <- rownames(seu@meta.data)[seu@meta.data[[cluster_col]] == cl]
  if (length(cells_keep) == 0) return(NULL)
  
  seu_sub <- subset(seu, cells = cells_keep)
  if (ncol(seu_sub) == 0) return(NULL)
  
  # Write group8 directly into metadata
  seu_sub$group8 <- paste(
    seu_sub@meta.data[[igra_col]],
    seu_sub@meta.data[[stim_col]],
    seu_sub@meta.data[[tp_col]],
    sep = "__"
  )
  
  # Enforce your 8-column ordering (drops anything unexpected)
  seu_sub$group8 <- factor(seu_sub$group8, levels = group_levels)
  
  # Drop cells not in those 8 groups (should be none, but safe)
  keep <- !is.na(seu_sub$group8)
  seu_sub <- seu_sub[, keep]
  if (ncol(seu_sub) == 0) return(NULL)
  
  # Must have at least 2 groups present
  if (nlevels(droplevels(seu_sub$group8)) < 2) return(NULL)
  
  Idents(seu_sub) <- droplevels(seu_sub$group8)
  
  mk <- FindAllMarkers(
    seu_sub,
    test.use = "MAST",          # or "wilcox"
    latent.vars = latent_vars,  # if MAST; remove if wilcox
    only.pos = FALSE,
    logfc.threshold = 0.1,
    min.pct = 0.1
  )
  
  if (is.null(mk) || nrow(mk) == 0) return(NULL)
  
  mk %>%
    dplyr::mutate(
      source_cluster = as.character(cl),
      abs_log2FC = abs(avg_log2FC)
    )
}


##########


# ---- MAIN: for each of 3 groups, build and save gene list + mapping ----
for (set_name in names(cluster_sets)) {
  message("[HEAVY] Gene selection for: ", set_name)
  
  cl_vec <- cluster_sets[[set_name]]
  
  mk_all <- purrr::map(cl_vec, ~ markers_group8_one_cluster(seu, cl = .x)) %>%
    dplyr::bind_rows()
  
  if (nrow(mk_all) == 0) {
    message("  -> No markers found for set ", set_name, ". Skipping.")
    next
  }
  
  # pick top genes PER group8 identity PER cluster by abs effect size
  mk_top <- mk_all %>%
    group_by(source_cluster, cluster) %>%  # cluster here = group8 identity
    slice_max(order_by = abs_log2FC, n = top_n_per_group8_per_cluster, with_ties = FALSE) %>%
    ungroup()
  
  # rank genes globally within this set by best abs effect size
  gene_rank <- mk_top %>%
    group_by(gene) %>%
    summarise(best_abs_log2FC = max(abs_log2FC, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(best_abs_log2FC))
  
  genes_final <- head(gene_rank$gene, genes_per_heatmap)
  
  # map each gene to the cluster that contributed its strongest signal (for Option B dotplot grouping)
  gene_to_source <- mk_top %>%
    filter(gene %in% genes_final) %>%
    group_by(gene) %>%
    slice_max(order_by = abs_log2FC, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      gene,
      source_cluster,
      best_abs_log2FC = abs_log2FC
    ) %>%
    arrange(match(gene, genes_final))
  
  # save outputs
  genes_out <- file.path(heat_marker_dir, paste0("HeatGenes__", set_name, "__n", genes_per_heatmap, ".csv"))
  map_out   <- file.path(heat_marker_dir, paste0("HeatGenes__", set_name, "__n", genes_per_heatmap, "__gene_to_source_cluster.csv"))
  
  write.csv(data.frame(gene = genes_final), genes_out, row.names = FALSE)
  write.csv(gene_to_source, map_out, row.names = FALSE)
  
  message("  Saved: ", basename(genes_out))
  message("  Saved: ", basename(map_out))
}

############################################################
# STEP B (LIGHT): Render 3 heatmaps + 3 dotplots from saved genes
# - NO rerun of FindAllMarkers
# - Heatmaps: SeuratExtend::CalcStats(method="zscore") + Heatmap
# - Dotplots: SeuratExtend::DotPlot2 with Option B grouping (source_cluster)
############################################################


# ---- cluster sets (3 plots total each type) ----
cluster_sets <- list(
  CD8_GD     = c("5","10","14","28"),
  CD4        = c("3","7","11","24"),
  Monocytes  = c("4","8","15")
)

# ---- dirs ----
heat_marker_dir <- file.path(heat_dir, "MarkerSelection")
plot_dir <- file.path(heat_dir, "Plots_SE")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- detect IGRA labels (same logic as Step A) ----
ig_levels_all <- levels(factor(seu@meta.data[[igra_col]]))
ig_pos_all <- ig_levels_all[grepl("\\+|pos", ig_levels_all, ignore.case = TRUE)][1]
ig_neg_all <- ig_levels_all[grepl("\\-|neg", ig_levels_all, ignore.case = TRUE)][1]

group_levels <- c(
  paste(ig_pos_all, "media_treated", "week12", sep="__"),
  paste(ig_neg_all, "media_treated", "week12", sep="__"),
  paste(ig_pos_all, "BCG_treated",   "week12", sep="__"),
  paste(ig_neg_all, "BCG_treated",   "week12", sep="__"),
  paste(ig_pos_all, "media_treated", "week44", sep="__"),
  paste(ig_neg_all, "media_treated", "week44", sep="__"),
  paste(ig_pos_all, "BCG_treated",   "week44", sep="__"),
  paste(ig_neg_all, "BCG_treated",   "week44", sep="__")
)

# ---- helper: make group8 in meta and subset to a cluster set ----
subset_with_group8 <- function(seu, cl_vec) {
  
  # robust subset by cells (avoid NSE issues)
  cells_keep <- rownames(seu@meta.data)[seu@meta.data[[cluster_col]] %in% cl_vec]
  if (length(cells_keep) == 0) return(NULL)
  
  seu_sub <- subset(seu, cells = cells_keep)
  if (ncol(seu_sub) == 0) return(NULL)
  
  # IMPORTANT: write group8 into Seurat object's metadata
  seu_sub$group8 <- paste(
    seu_sub@meta.data[[igra_col]],
    seu_sub@meta.data[[stim_col]],
    seu_sub@meta.data[[tp_col]],
    sep = "__"
  )
  
  seu_sub$group8 <- factor(seu_sub$group8, levels = group_levels)
  
  keep <- !is.na(seu_sub$group8)
  seu_sub <- seu_sub[, keep]
  if (ncol(seu_sub) == 0) return(NULL)
  
  return(seu_sub)
}


# ---- helper: render heatmap (zscore) from gene list ----
render_heatmap <- function(seu_sub, genes, set_name) {
  
  zmat <- SeuratExtend::CalcStats(
    seu_sub,
    features = genes,
    group.by = "group8",
    slot = "data",
    assay = DefaultAssay(seu_sub),
    method = "zscore",
    n = Inf
  )
  
  # reorder columns manually
  cols_keep <- intersect(group_levels, colnames(zmat))
  zmat <- zmat[, cols_keep, drop = FALSE]
  
  p <- SeuratExtend::Heatmap(
    zmat,
    lab_fill = "zscore",
    color_scheme = "A"
  ) +
    ggtitle(paste0(set_name, " | zscore by group8")) +
    theme(
      plot.title  = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 8)  # helps fit 100 genes
    )
  
  p <- p + coord_cartesian(clip = "off")
  
  out_png <- file.path(plot_dir, paste0("Heatmap_SE__", set_name, ".png"))
  
  ggsave(
    filename = out_png,
    plot = p,
    width = 12, height = 15, units = "in", dpi = 300,
  )
}

# ---- helper: render dotplot with Option B grouping by source_cluster ----
render_dotplot_optionB <- function(seu_sub, gene_to_source, set_name) {
  
  Idents(seu_sub) <- droplevels(seu_sub@meta.data$group8)
  
  grouped_features <- split(gene_to_source$gene, gene_to_source$source_cluster)
  
  p <- SeuratExtend::DotPlot2(seu_sub, features = grouped_features) +
    ggtitle(paste0(set_name, " | DotPlot2 (group8) | grouped by source cluster")) +
    theme(
      plot.title  = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 8)
    )
  
  out_png <- file.path(plot_dir, paste0("DotPlot2_SE__", set_name, ".png"))
  
  ggsave(
    filename = out_png,
    plot = p,
    width = 12, height = 14, units = "in", dpi = 300, limitsize = FALSE
  )
}


# ---- MAIN: loop 3 groups, draw heatmap + dotplot each ----
# =========================
# RUN HEATMAPS ONLY
# =========================
for (set_name in names(cluster_sets)) {
  
  genes_file <- list.files(
    heat_marker_dir,
    pattern = paste0("^HeatGenes__", set_name, "__n[0-9]+\\.csv$"),
    full.names = TRUE
  )
  if (length(genes_file) != 1) {
    message("[HEATMAP] Skipping ", set_name, ": gene list not found (or multiple) in ", heat_marker_dir)
    next
  }
  
  genes <- readr::read_csv(genes_file, show_col_types = FALSE)$gene
  
  message("[HEATMAP] Rendering: ", set_name, " (", length(genes), " genes)")
  
  seu_sub <- subset_with_group8(seu, cluster_sets[[set_name]])
  if (is.null(seu_sub)) {
    message("  -> No cells for set ", set_name, "; skipping.")
    next
  }
  
  render_heatmap(seu_sub, genes, set_name)
}
# =========================
# RUN DOTPLOTS ONLY
# =========================
for (set_name in names(cluster_sets)) {
  
  genes_file <- list.files(
    heat_marker_dir,
    pattern = paste0("^HeatGenes__", set_name, "__n[0-9]+\\.csv$"),
    full.names = TRUE
  )
  map_file <- list.files(
    heat_marker_dir,
    pattern = paste0("^HeatGenes__", set_name, "__n[0-9]+__gene_to_source_cluster\\.csv$"),
    full.names = TRUE
  )
  
  if (length(genes_file) != 1 || length(map_file) != 1) {
    message("[DOTPLOT] Skipping ", set_name, ": gene list or mapping file not found (or multiple) in ", heat_marker_dir)
    next
  }
  
  genes <- readr::read_csv(genes_file, show_col_types = FALSE)$gene
  gene_to_source <- readr::read_csv(map_file, show_col_types = FALSE)
  
  # ensure mapping only contains genes in the list (safety)
  gene_to_source <- gene_to_source %>% filter(gene %in% genes)
  
  message("[DOTPLOT] Rendering: ", set_name, " (", length(genes), " genes)")
  
  seu_sub <- subset_with_group8(seu, cluster_sets[[set_name]])
  if (is.null(seu_sub)) {
    message("  -> No cells for set ", set_name, "; skipping.")
    next
  }
  
  render_dotplot_optionB(seu_sub, gene_to_source, set_name)
}


