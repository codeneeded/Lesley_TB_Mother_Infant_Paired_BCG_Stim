###############################################################################
## BCG Parse — DGE + condition heatmaps   (runs on the ANNOTATED object)
##
## Input : BCG_PARSE_integrated_annotated.qs2  (from BCG_PARSE_annotation.R)
## Group : cell_type (default) — pools same-named clusters and keeps rare types
##         testable in the 8-way IGRA split. Set group_col <- "cluster_celltype"
##         for a per-cluster pass, or "mnn.snn.louvianmlr_1" for raw clusters.
##
## Contrasts (per group, MAST):
##   3A  BCG_treated vs media_treated, within each timepoint
##   3B  week44 vs week12, within each stim
##   3C  IGRA+ vs IGRA-, within each (stim x timepoint)
##   4   condition heatmaps: genes varying across the 8 IGRA x stim x timepoint
##       groups, per lineage set (sets derived from the annotation)
##
## Latent vars: nCount_RNA always; sample_id (donor) added for 3A/3B where both
## arms share donors, but NOT for 3C — donor is collinear with IGRA status there.
###############################################################################

library(Seurat)
library(SeuratExtend)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(readr)
library(qs2)
# test.use = "MAST" needs the Bioconductor 'MAST' package installed.

# ============================ Config =========================================
base_dir  <- "/home/akshay-iyer/Documents/BCG_Stim_Infants_12_44_wks_PARSE_GEX"
saved_dir <- file.path(base_dir, "saved_R_data")
dge_dir   <- file.path(base_dir, "DGE")
heat_dir  <- file.path(base_dir, "Heatmaps", "IGRA_Stim_Timepoint")
for (d in c(dge_dir, heat_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

group_col   <- "cell_type"       # "cell_type" | "cluster_celltype" | "mnn.snn.louvianmlr_1"
igra_col    <- "IGRA_status"
stim_col    <- "stim"
tp_col      <- "timepoint"
sample_col  <- "sample_id"
min_cells   <- 50
run_heatmaps <- TRUE

safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)   # filename-safe labels

# ============================ Load annotated =================================
annot_path <- file.path(saved_dir, "BCG_PARSE_integrated_annotated.qs2")
if (!file.exists(annot_path))
  stop("Annotated object not found: ", annot_path, "  — run BCG_PARSE_annotation.R (PHASE 2) first.")
seu <- qs_read(annot_path)

DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu)
seu$stim      <- factor(seu$stim)
seu$timepoint <- factor(seu$timepoint)
stopifnot(group_col %in% colnames(seu@meta.data))

# ---- IGRA mapping (per donor) ----
igra_map <- data.frame(
  sample_id = c("6044931","6063411","6063421","6063481","6061161",
                "6062151","6065761","6067021","6071261","6044981"),
  IGRA_status = c("Negative","Negative","Negative","Negative","Positive",
                  "Positive","Positive","Positive","Positive","Positive"),
  stringsAsFactors = FALSE)
seu$IGRA_status <- igra_map$IGRA_status[match(as.character(seu$sample_id), igra_map$sample_id)]
stopifnot(all(c(group_col, igra_col, stim_col, tp_col, sample_col) %in% colnames(seu@meta.data)))

groups     <- levels(factor(seu@meta.data[[group_col]]))
timepoints <- levels(seu$timepoint)
stims      <- levels(seu$stim)

# ============================================================================
# Generic 2-group contrast within a stratum, for one cell group.
#   g        : value of group_col (one cell type)
#   strat    : named list of metadata col = value conditions defining the stratum
#   test_col : metadata column holding the two arms
#   id1,id2  : the two arms (ident.1 vs ident.2)
#   latent   : latent vars (kept only if present in metadata)
# Returns the markers data.frame (also written to outdir/fname), or NULL if the
# stratum is missing an arm or too small. Flat and individually testable.
# ============================================================================
run_contrast <- function(g, strat, test_col, id1, id2, latent, outdir, fname, extra = list()) {
  keep <- as.character(seu@meta.data[[group_col]]) == g
  for (nm in names(strat)) keep <- keep & as.character(seu@meta.data[[nm]]) == strat[[nm]]
  cells <- rownames(seu@meta.data)[keep]
  if (length(cells) < 2 * min_cells) return(NULL)
  
  sub <- subset(seu, cells = cells)
  arm <- droplevels(factor(sub@meta.data[[test_col]]))
  if (!all(c(id1, id2) %in% levels(arm))) return(NULL)
  tab <- table(arm)[c(id1, id2)]
  if (any(tab < min_cells)) return(NULL)
  sub@meta.data[[test_col]] <- arm
  
  lat <- intersect(latent, colnames(sub@meta.data))
  mk <- FindMarkers(sub, ident.1 = id1, ident.2 = id2, group.by = test_col,
                    test.use = "MAST", latent.vars = lat,
                    logfc.threshold = 0.1, min.pct = 0.1)
  if (is.null(mk) || nrow(mk) == 0) return(NULL)
  
  mk <- tibble::rownames_to_column(mk, "gene")
  mk[[group_col]] <- g
  for (nm in names(extra)) mk[[nm]] <- extra[[nm]]
  write.csv(mk, file.path(outdir, fname), row.names = FALSE)
  mk
}

# ---------------- 3A) BCG vs media, within each timepoint --------------------
dir3a <- file.path(dge_dir, "BCG_vs_media_byTimepoint"); dir.create(dir3a, showWarnings = FALSE)
res3a <- list()
for (g in groups) for (tp in timepoints) {
  message("3A | ", g, " | ", tp)
  res3a[[length(res3a) + 1]] <- run_contrast(
    g, strat = setNames(list(tp), tp_col), test_col = stim_col,
    id1 = "BCG_treated", id2 = "media_treated",
    latent = c("nCount_RNA", "sample_id"), outdir = dir3a,
    fname = paste0("DGE_", safe(g), "_", tp, "_BCG_vs_media.csv"),
    extra = list(timepoint = tp, comparison = "BCG_vs_media"))
}
res3a <- bind_rows(res3a)
if (nrow(res3a)) write.csv(res3a, file.path(dir3a, "ALL__BCG_vs_media_byTimepoint_LONG.csv"), row.names = FALSE)

# ---------------- 3B) week44 vs week12, within each stim ---------------------
dir3b <- file.path(dge_dir, "Week44_vs_Week12_byStim"); dir.create(dir3b, showWarnings = FALSE)
res3b <- list()
for (g in groups) for (st in stims) {
  message("3B | ", g, " | ", st)
  res3b[[length(res3b) + 1]] <- run_contrast(
    g, strat = setNames(list(st), stim_col), test_col = tp_col,
    id1 = "week44", id2 = "week12",
    latent = c("nCount_RNA", "sample_id"), outdir = dir3b,
    fname = paste0("DGE_", safe(g), "_", safe(st), "_week44_vs_week12.csv"),
    extra = list(stim = st, comparison = "week44_vs_week12"))
}
res3b <- bind_rows(res3b)
if (nrow(res3b)) write.csv(res3b, file.path(dir3b, "ALL__Week44_vs_Week12_byStim_LONG.csv"), row.names = FALSE)

# ---------------- 3C) IGRA+ vs IGRA-, within each (stim x timepoint) ---------
dir3c <- file.path(dge_dir, "IGRApos_vs_IGRAneg_byStimTimepoint"); dir.create(dir3c, showWarnings = FALSE)
ig_lv  <- levels(factor(seu@meta.data[[igra_col]]))
ig_pos <- ig_lv[grepl("\\+|pos", ig_lv, ignore.case = TRUE)][1]
ig_neg <- ig_lv[grepl("\\-|neg", ig_lv, ignore.case = TRUE)][1]
if (is.na(ig_pos) || is.na(ig_neg)) stop("Can't detect IGRA+/- labels: ", paste(ig_lv, collapse = ", "))

comparisons <- tibble::tribble(
  ~stim_level,     ~tp_level, ~label,
  "media_treated", "week12",  "IGRApos_vs_neg_Med_12wk",
  "BCG_treated",   "week12",  "IGRApos_vs_neg_BCG_12wk",
  "media_treated", "week44",  "IGRApos_vs_neg_Med_44wk",
  "BCG_treated",   "week44",  "IGRApos_vs_neg_BCG_44wk")

res3c <- list()
for (g in groups) for (i in seq_len(nrow(comparisons))) {
  cc <- comparisons[i, ]
  message("3C | ", g, " | ", cc$label)
  res3c[[length(res3c) + 1]] <- run_contrast(
    g, strat = setNames(list(cc$stim_level, cc$tp_level), c(stim_col, tp_col)),
    test_col = igra_col, id1 = ig_pos, id2 = ig_neg,
    latent = c("nCount_RNA"),                # NOT sample_id (collinear with IGRA)
    outdir = dir3c, fname = paste0("DGE_", safe(g), "__", cc$label, ".csv"),
    extra = list(stim = cc$stim_level, timepoint = cc$tp_level, comparison = cc$label))
}
res3c <- bind_rows(res3c)
if (nrow(res3c)) write.csv(res3c, file.path(dir3c, "ALL__IGRApos_vs_IGRAneg_LONG.csv"), row.names = FALSE)

message("DGE contrasts done.")

###############################################################################
# 4) Condition heatmaps: genes varying across the 8 IGRA x stim x timepoint
#    groups, per lineage set. Sets are DERIVED from the annotation (not hardcoded).
###############################################################################
if (run_heatmaps) {
  
  cts <- groups
  cluster_sets <- list(
    T_CD8_gdT = grep("CD8|gdT|MAIT", cts, value = TRUE, ignore.case = TRUE),
    T_CD4     = grep("CD4|Treg",     cts, value = TRUE, ignore.case = TRUE),
    NK        = grep("NK",           cts, value = TRUE, ignore.case = TRUE),
    B_Plasma  = grep("^B |Plasm",    cts, value = TRUE, ignore.case = TRUE),
    Myeloid   = grep("Mono|DC",      cts, value = TRUE, ignore.case = TRUE))
  cluster_sets <- cluster_sets[lengths(cluster_sets) > 0]
  message("Heatmap sets:\n", paste(names(cluster_sets),
                                   sapply(cluster_sets, paste, collapse = "; "), sep = ": ", collapse = "\n"))
  
  group_levels <- as.vector(t(outer(
    c(ig_pos, ig_neg),
    c(paste(c("media_treated","BCG_treated"), "week12", sep = "__"),
      paste(c("media_treated","BCG_treated"), "week44", sep = "__")),
    paste, sep = "__")))  # 8 ordered IGRA__stim__timepoint columns
  
  add_group8 <- function(obj) {
    obj$group8 <- factor(paste(obj@meta.data[[igra_col]], obj@meta.data[[stim_col]],
                               obj@meta.data[[tp_col]], sep = "__"), levels = group_levels)
    obj[, !is.na(obj$group8)]
  }
  
  genes_per_heatmap <- 100
  top_per_group8    <- 6
  
  for (set_name in names(cluster_sets)) {
    members <- cluster_sets[[set_name]]
    cells <- rownames(seu@meta.data)[as.character(seu@meta.data[[group_col]]) %in% members]
    if (length(cells) < 2 * min_cells) { message("[heatmap] skip ", set_name, " (too few cells)"); next }
    
    sub <- add_group8(subset(seu, cells = cells))
    if (ncol(sub) == 0 || nlevels(droplevels(sub$group8)) < 2) { message("[heatmap] skip ", set_name); next }
    Idents(sub) <- droplevels(sub$group8)
    
    mk <- FindAllMarkers(sub, test.use = "MAST", latent.vars = c("nCount_RNA"),
                         only.pos = FALSE, logfc.threshold = 0.1, min.pct = 0.1)
    if (is.null(mk) || nrow(mk) == 0) { message("[heatmap] no markers ", set_name); next }
    
    genes <- mk %>% mutate(a = abs(avg_log2FC)) %>% group_by(cluster) %>%
      slice_max(a, n = top_per_group8, with_ties = FALSE) %>% ungroup() %>%
      group_by(gene) %>% summarise(a = max(a), .groups = "drop") %>%
      arrange(desc(a)) %>% slice_head(n = genes_per_heatmap) %>% pull(gene)
    
    z <- SeuratExtend::CalcStats(sub, features = genes, group.by = "group8",
                                 slot = "data", method = "zscore", n = Inf)
    z <- z[, intersect(group_levels, colnames(z)), drop = FALSE]
    
    p <- SeuratExtend::Heatmap(z, lab_fill = "zscore", color_scheme = "A") +
      ggtitle(paste0(set_name, " | zscore by IGRA x stim x timepoint")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 7))
    ggsave(file.path(heat_dir, paste0("Heatmap__", set_name, ".png")),
           p, width = 11, height = 15, dpi = 300, bg = "white", limitsize = FALSE)
    write.csv(data.frame(gene = genes),
              file.path(heat_dir, paste0("HeatGenes__", set_name, ".csv")), row.names = FALSE)
    message("[heatmap] ", set_name, " (", length(genes), " genes)")
  }
}

message("All done.")
###############################################################################