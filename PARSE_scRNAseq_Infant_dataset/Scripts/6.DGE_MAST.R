###############################################################################
## BCG Parse — IGRA+ vs IGRA-  (collaborator request: 7 cell-type sets, 12 & 44 wk)
##
## Input : BCG_PARSE_integrated_annotated.qs2  (from BCG_PARSE_annotation.R PHASE 2)
##
## PRIMARY  : pseudobulk per donor (sum raw counts per donor x stim x timepoint
##            x set), DESeq2 ~ IGRA_status, run WITHIN BCG and WITHIN media
##            separately, at 12 wk and 44 wk separately. IGRA is a DONOR-level
##            variable (4 IGRA- vs 6 IGRA+ donors), so the donor is the unit of
##            replication — pseudobulk is the defensible test here.
## SPECIFIC : flag genes significant in BCG but NOT in media (or opposite sign)
##            => BCG-response-specific IGRA signal (collaborator's logic).
## SENSITIV.: MAST on cells (BCG only), restricted to the same sets, with a
##            per-donor detection check to expose pseudoreplicated 1-donor hits.
## PLOTS    : SeuratExtend CalcStats + Heatmap of BCG-specific genes across the
##            4 IGRA x stim groups within each timepoint, plus VlnPlot2 of tops.
##
## RUN ORDER: SECTION 0 first (every session), then 1 -> 2 -> 3 -> 4 -> 5.
## Each section is flat and re-runnable; SECTIONS 2/4 reload nothing destructive.
###############################################################################


## ============================================================================
## SECTION 0 — SETUP  (RUN FIRST)
## ============================================================================
library(Seurat)
library(SeuratExtend)
library(Matrix)
library(DESeq2)          # install.packages via BiocManager::install("DESeq2")
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(ggplot2)
library(qs2)
# MAST (SECTION 4) needs BiocManager::install("MAST")

base_dir   <- "/home/akshay-iyer/Documents/BCG_Stim_Infants_12_44_wks_PARSE_GEX"
saved_dir  <- file.path(base_dir, "saved_R_data")
out_dir    <- file.path(base_dir, "DGE", "IGRA_pos_vs_neg_7sets")
pb_dir     <- file.path(out_dir, "Pseudobulk_DESeq2")
spec_dir   <- file.path(out_dir, "BCG_response_specific")
mast_dir   <- file.path(out_dir, "MAST_sensitivity")
plot_dir   <- file.path(out_dir, "Plots")
for (d in c(pb_dir, spec_dir, mast_dir, plot_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

group_col  <- "cell_type"
igra_col   <- "IGRA_status"
stim_col   <- "stim"
tp_col     <- "timepoint"
sample_col <- "sample_id"

# thresholds
min_cells_per_donor   <- 10    # drop a donor's pseudobulk built from < this many cells
min_donors_per_group  <- 2     # need >= this many donors per IGRA arm to run DESeq2 (3+ is safer)
padj_thresh           <- 0.05
lfc_thresh            <- 0.585 # |log2FC| >= 0.585 ~ 1.5x
mast_min_cells        <- 30    # per IGRA arm, for the cell-level MAST sensitivity

safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

## ---- analysis units -> exact cell_type levels they pool ----
## EDIT these strings to match levels(seu$cell_type) on your object if needed.
target_sets <- list(
  ## --- requested by collaborator (combined as specified) ---
  "CD4_Naive_T"    = c("CD4 Naive T"),
  "CD4_TCM"        = c("CD4 TCM"),
  "CD8_Naive_T"    = c("CD8 Naive T"),
  "CD8_TCM_TEM"    = c("CD8 TCM", "CD8 TEM"),
  "Treg"           = c("Treg"),
  "gdT"            = c("gdT (Vd1)", "gdT (cytotoxic)", "gdT (Vg9Vd2)"),
  "CD14_Mono_Infl" = c("CD14 Monocyte", "CD14 Monocyte (inflammatory)"),
  ## --- remaining populations (lineage-level combining for power at N=4 vs 6) ---
  "MAIT"           = c("MAIT"),
  "NK"             = c("NK", "NK CD56bright", "NK (activated)"),
  "B"              = c("B naive", "B intermediate", "B memory"),
  "Plasmablast"    = c("Plasmablast"),
  "DC"             = c("pDC", "Activated DC")
)
## To split any combined set (e.g. NK -> CD56dim/bright/activated), replace its
## one line with separate entries — only that set re-runs, not the whole script.
## Intentionally EXCLUDED: cluster 15 "CD4 T (activated/stress)" (QC-flagged:
## IEG/MALAT1-hi, high-mito). Add it as a set only if you chose to keep it.

## ---- IGRA mapping (per donor) ----
igra_map <- data.frame(
  sample_id = c("6044931","6063411","6063421","6063481","6061161",
                "6062151","6065761","6067021","6071261","6044981"),
  IGRA_status = c("Negative","Negative","Negative","Negative","Positive",
                  "Positive","Positive","Positive","Positive","Positive"),
  stringsAsFactors = FALSE)

## ---- load annotated object ----
annot_path <- file.path(saved_dir, "BCG_PARSE_integrated_annotated.qs2")
if (!file.exists(annot_path))
  stop("Annotated object not found: ", annot_path, " — run BCG_PARSE_annotation.R PHASE 2 first.")
seu <- qs_read(annot_path)
DefaultAssay(seu) <- "RNA"
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])

seu$IGRA_status <- igra_map$IGRA_status[match(as.character(seu@meta.data[[sample_col]]), igra_map$sample_id)]
seu$IGRA_status <- factor(seu$IGRA_status, levels = c("Negative","Positive"))
seu$stim        <- factor(seu$stim)
seu$timepoint   <- factor(seu$timepoint)

## ---- validate requested cell_type levels exist, then build target_set column ----
have_levels <- levels(factor(seu@meta.data[[group_col]]))
want_levels <- unique(unlist(target_sets, use.names = FALSE))
missing_lv  <- setdiff(want_levels, have_levels)
if (length(missing_lv)) {
  message("Available cell_type levels:\n  ", paste(have_levels, collapse = "\n  "))
  stop("These requested cell_type levels were not found (edit target_sets to match): ",
       paste(missing_lv, collapse = " | "))
}

ct2set <- unlist(lapply(names(target_sets),
                        function(s) setNames(rep(s, length(target_sets[[s]])), target_sets[[s]])))
ct2set <- ct2set[!duplicated(names(ct2set))]
seu$target_set <- unname(ct2set[as.character(seu@meta.data[[group_col]])])  # NA if not requested

sets       <- names(target_sets)
timepoints <- levels(seu$timepoint)
stims      <- levels(seu$stim)
message("SETUP done. Sets: ", paste(sets, collapse = ", "),
        " | timepoints: ", paste(timepoints, collapse = ", "),
        " | stims: ", paste(stims, collapse = ", "))


## ============================================================================
## SECTION 1 — BUILD PSEUDOBULK MATRIX  (donor x stim x timepoint x set)
## ============================================================================
cnt <- LayerData(seu, assay = "RNA", layer = "counts")     # genes x cells (sparse)
md  <- seu@meta.data
md$pb_id <- paste(md$target_set, md[[sample_col]], as.character(md$stim),
                  as.character(md$timepoint), sep = "~~")

keep_cell <- !is.na(md$target_set)
md_k <- md[keep_cell, , drop = FALSE]
cnt_k <- cnt[, keep_cell, drop = FALSE]

pb_ids <- unique(md_k$pb_id)
pb_mat <- vapply(pb_ids,
                 function(id) Matrix::rowSums(cnt_k[, md_k$pb_id == id, drop = FALSE]),
                 numeric(nrow(cnt_k)))
rownames(pb_mat) <- rownames(cnt_k)
colnames(pb_mat) <- pb_ids
pb_mat <- round(pb_mat)

parts   <- do.call(rbind, strsplit(pb_ids, "~~", fixed = TRUE))
pb_meta <- data.frame(
  pb_id      = pb_ids,
  target_set = parts[, 1],
  sample_id  = parts[, 2],
  stim       = parts[, 3],
  timepoint  = parts[, 4],
  stringsAsFactors = FALSE)
pb_meta$n_cells     <- as.integer(table(md_k$pb_id)[pb_meta$pb_id])
pb_meta$IGRA_status <- igra_map$IGRA_status[match(pb_meta$sample_id, igra_map$sample_id)]

write.csv(pb_meta, file.path(pb_dir, "pseudobulk_sample_table.csv"), row.names = FALSE)
qs_save(list(pb_mat = pb_mat, pb_meta = pb_meta), file.path(pb_dir, "pseudobulk_checkpoint.qs2"))
message("SECTION 1 done. ", ncol(pb_mat), " pseudobulk samples x ", nrow(pb_mat), " genes.")
print(pb_meta %>% count(target_set, stim, timepoint, IGRA_status) %>%
        pivot_wider(names_from = IGRA_status, values_from = n, values_fill = 0))


## ============================================================================
## SECTION 2 — PSEUDOBULK DESeq2  (IGRA+ vs IGRA-, within BCG and within media)
## ============================================================================
run_pb_deseq <- function(set_name, tp, stim_level) {
  sel <- pb_meta$target_set == set_name & pb_meta$timepoint == tp &
    pb_meta$stim == stim_level & pb_meta$n_cells >= min_cells_per_donor
  cd <- pb_meta[sel, , drop = FALSE]
  if (nrow(cd) < 2) return(NULL)
  cd$IGRA_status <- factor(cd$IGRA_status, levels = c("Negative","Positive"))
  ng <- table(cd$IGRA_status)
  if (length(ng) < 2 || any(ng < min_donors_per_group)) {
    message("  skip ", set_name, " | ", tp, " | ", stim_level,
            "  (donors: -", ng["Negative"], " / +", ng["Positive"], ")")
    return(NULL)
  }
  mat <- pb_mat[, cd$pb_id, drop = FALSE]
  keep_g <- rowSums(mat >= 1) >= 2 & rowSums(mat) >= 10
  mat <- mat[keep_g, , drop = FALSE]
  if (nrow(mat) < 10) return(NULL)
  
  res <- tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat, colData = cd, design = ~ IGRA_status)
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    DESeq2::results(dds, contrast = c("IGRA_status", "Positive", "Negative"))
  }, error = function(e) { message("  DESeq2 error (", set_name, "|", tp, "|", stim_level, "): ", conditionMessage(e)); NULL })
  if (is.null(res)) return(NULL)
  
  out <- as.data.frame(res)
  out <- tibble::rownames_to_column(out, "gene")
  out$set <- set_name; out$timepoint <- tp; out$stim <- stim_level
  out$n_pos <- as.integer(ng["Positive"]); out$n_neg <- as.integer(ng["Negative"])
  out <- out[order(out$padj), ]
  write.csv(out, file.path(pb_dir, paste0("PB_", set_name, "_", tp, "_", safe(stim_level), ".csv")),
            row.names = FALSE)
  out
}

res_pb <- list()
for (s in sets) for (tp in timepoints) for (st in stims) {
  message("2 | ", s, " | ", tp, " | ", st)
  res_pb[[paste(s, tp, st, sep = "|")]] <- run_pb_deseq(s, tp, st)
}
res_pb_long <- bind_rows(Filter(Negate(is.null), res_pb))
if (nrow(res_pb_long))
  write.csv(res_pb_long, file.path(pb_dir, "ALL__pseudobulk_DESeq2_LONG.csv"), row.names = FALSE)
message("SECTION 2 done — pseudobulk DESeq2 results in ", pb_dir)


## ============================================================================
## SECTION 3 — BCG-RESPONSE-SPECIFIC FLAG  (sig in BCG, not media / opp. sign)
## ============================================================================
spec_all <- list()
for (s in sets) for (tp in timepoints) {
  bcg <- res_pb[[paste(s, tp, "BCG_treated",   sep = "|")]]
  med <- res_pb[[paste(s, tp, "media_treated", sep = "|")]]
  if (is.null(bcg)) { message("3 | no BCG result: ", s, " | ", tp); next }
  
  bcg2 <- bcg %>% transmute(gene, lfc_bcg = log2FoldChange, padj_bcg = padj)
  if (is.null(med)) {
    med2 <- tibble(gene = bcg2$gene, lfc_med = NA_real_, padj_med = NA_real_)
    message("3 | ", s, " | ", tp, "  (no media result — specificity vs media unverified)")
  } else {
    med2 <- med %>% transmute(gene, lfc_med = log2FoldChange, padj_med = padj)
  }
  
  m <- bcg2 %>% left_join(med2, by = "gene")
  m$sig_bcg   <- !is.na(m$padj_bcg) & m$padj_bcg < padj_thresh & abs(m$lfc_bcg) >= lfc_thresh
  m$sig_media <- !is.na(m$padj_med) & m$padj_med < padj_thresh & abs(m$lfc_med) >= lfc_thresh
  m$opp_sign  <- !is.na(m$lfc_med) & sign(m$lfc_bcg) != sign(m$lfc_med)
  m$bcg_specific <- m$sig_bcg & (!m$sig_media | m$opp_sign)
  m$set <- s; m$timepoint <- tp
  m <- m %>% arrange(desc(bcg_specific), padj_bcg)
  
  write.csv(m, file.path(spec_dir, paste0("SPEC_", s, "_", tp, ".csv")), row.names = FALSE)
  spec_all[[paste(s, tp)]] <- m %>% filter(bcg_specific)
}
spec_long <- bind_rows(spec_all)
if (nrow(spec_long))
  write.csv(spec_long, file.path(spec_dir, "ALL__BCG_response_specific_LONG.csv"), row.names = FALSE)
message("SECTION 3 done. BCG-response-specific genes per set x timepoint:")
print(spec_long %>% count(set, timepoint, name = "n_bcg_specific"))


## ============================================================================
## SECTION 4 — MAST SENSITIVITY  (cells, BCG only; flags 1-donor pseudoreplication)
## ============================================================================
run_mast_igra <- function(set_name, tp) {
  keep <- !is.na(seu$target_set) & seu$target_set == set_name &
    as.character(seu$timepoint) == tp & as.character(seu$stim) == "BCG_treated"
  cells <- colnames(seu)[keep]
  if (length(cells) < 2 * mast_min_cells) return(NULL)
  sub <- subset(seu, cells = cells)
  ig  <- droplevels(factor(sub$IGRA_status))
  if (!all(c("Positive","Negative") %in% levels(ig))) return(NULL)
  if (any(table(ig) < mast_min_cells)) return(NULL)
  sub$IGRA_status <- ig
  
  mk <- FindMarkers(sub, ident.1 = "Positive", ident.2 = "Negative", group.by = "IGRA_status",
                    test.use = "MAST", latent.vars = "nCount_RNA",
                    logfc.threshold = 0.1, min.pct = 0.1)
  if (is.null(mk) || nrow(mk) == 0) return(NULL)
  mk <- tibble::rownames_to_column(mk, "gene")
  
  # per-donor detection for the significant hits (>1 donor in each arm = trustworthy)
  sig <- mk$gene[!is.na(mk$p_val_adj) & mk$p_val_adj < padj_thresh]
  if (length(sig)) {
    dat   <- LayerData(sub, assay = "RNA", layer = "data")[sig, , drop = FALSE]
    donor <- as.character(sub@meta.data[[sample_col]])
    armv  <- as.character(sub$IGRA_status)
    det   <- (dat > 0)
    ndon  <- function(arm) vapply(sig, function(g) {
      d <- unique(donor[armv == arm & det[g, ]]); length(d)
    }, integer(1))
    dpos <- ndon("Positive"); dneg <- ndon("Negative")
    mk$n_donors_pos_detect <- NA_integer_; mk$n_donors_neg_detect <- NA_integer_
    mk$n_donors_pos_detect[match(sig, mk$gene)] <- dpos
    mk$n_donors_neg_detect[match(sig, mk$gene)] <- dneg
    mk$one_donor_flag <- NA
    mk$one_donor_flag[match(sig, mk$gene)] <- (pmin(dpos, dneg) <= 1)
  }
  mk$set <- set_name; mk$timepoint <- tp
  write.csv(mk, file.path(mast_dir, paste0("MAST_", set_name, "_", tp, "_BCG.csv")), row.names = FALSE)
  mk
}

res_mast <- list()
for (s in sets) for (tp in timepoints) {
  message("4 | ", s, " | ", tp)
  res_mast[[paste(s, tp)]] <- run_mast_igra(s, tp)
}

# overlap: pseudobulk(BCG) vs MAST, per set x timepoint
overlap <- list()
for (s in sets) for (tp in timepoints) {
  pb <- res_pb[[paste(s, tp, "BCG_treated", sep = "|")]]
  ms <- res_mast[[paste(s, tp)]]
  pb_sig <- if (!is.null(pb)) pb$gene[!is.na(pb$padj) & pb$padj < padj_thresh & abs(pb$log2FoldChange) >= lfc_thresh] else character(0)
  ms_sig <- if (!is.null(ms)) ms$gene[!is.na(ms$p_val_adj) & ms$p_val_adj < padj_thresh] else character(0)
  overlap[[paste(s, tp)]] <- tibble(
    set = s, timepoint = tp,
    n_pseudobulk_sig = length(pb_sig),
    n_MAST_sig       = length(ms_sig),
    n_shared         = length(intersect(pb_sig, ms_sig)))
}
overlap_tbl <- bind_rows(overlap)
write.csv(overlap_tbl, file.path(mast_dir, "OVERLAP_pseudobulk_vs_MAST.csv"), row.names = FALSE)
message("SECTION 4 done. Pseudobulk vs MAST overlap (large MAST/pb gap = pseudoreplication):")
print(overlap_tbl)


## ============================================================================
## SECTION 5 — PLOTS (SeuratExtend): BCG-specific genes across IGRA x stim
## ============================================================================
plot_set_tp <- function(set_name, tp, top_n = 40, n_vln = 6) {
  genes <- spec_all[[paste(set_name, tp)]]
  if (is.null(genes) || nrow(genes) == 0) { message("5 | no BCG-specific genes: ", set_name, " | ", tp); return(invisible(NULL)) }
  genes <- genes %>% arrange(desc(abs(lfc_bcg))) %>% slice_head(n = top_n) %>% pull(gene)
  
  keep  <- !is.na(seu$target_set) & seu$target_set == set_name & as.character(seu$timepoint) == tp
  cells <- colnames(seu)[keep]
  if (length(cells) < 2 * mast_min_cells) return(invisible(NULL))
  sub <- subset(seu, cells = cells)
  sub$grp4 <- factor(paste(sub$IGRA_status, sub$stim, sep = " | "))
  
  # heatmap of z-scored mean expression across the 4 IGRA x stim groups
  z <- SeuratExtend::CalcStats(sub, features = genes, group.by = "grp4",
                               slot = "data", method = "zscore", n = Inf)
  p_heat <- SeuratExtend::Heatmap(z, lab_fill = "zscore", color_scheme = "A") +
    ggtitle(paste0(set_name, " | ", tp, " | BCG-specific IGRA genes")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 7))
  ggsave(file.path(plot_dir, paste0("Heatmap_", set_name, "_", tp, ".png")),
         p_heat, width = 9, height = 12, dpi = 300, bg = "white", limitsize = FALSE)
  
  # violins for the top genes, BCG cells only, split by IGRA
  bcg_cells <- colnames(sub)[as.character(sub$stim) == "BCG_treated"]
  if (length(bcg_cells) >= 2 * mast_min_cells) {
    sub_bcg <- subset(sub, cells = bcg_cells)
    vg <- head(genes, n_vln)
    p_vln <- SeuratExtend::VlnPlot2(sub_bcg, features = vg, group.by = "IGRA_status",
                                    cols = "default", show.mean = TRUE, stat.method = "wilcox.test") +
      ggtitle(paste0(set_name, " | ", tp, " | BCG cells, IGRA+ vs IGRA-"))
    ggsave(file.path(plot_dir, paste0("Vln_", set_name, "_", tp, ".png")),
           p_vln, width = 14, height = 8, dpi = 300, bg = "white")
  }
  message("5 | plotted ", set_name, " | ", tp, " (", length(genes), " genes)")
}

for (s in sets) for (tp in timepoints) plot_set_tp(s, tp)
message("SECTION 5 done — plots under ", plot_dir)
message("ALL DONE.")
###############################################################################