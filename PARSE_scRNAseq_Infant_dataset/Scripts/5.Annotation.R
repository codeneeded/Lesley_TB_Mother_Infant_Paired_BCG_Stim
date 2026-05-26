###############################################################################
## BCG Parse — ANNOTATION   (tables -> annotate -> labeled plots)
##
## RUN ORDER (each block is independent and flat — no nested wrappers):
##   1. SETUP    : libs + config + gene panels.  RUN THIS FIRST, every session.
##   2. PHASE 1  : write annotation tables + template, checkpoint object.
##                 -> then fill Tables/ANNOTATION_TEMPLATE_fill_me.csv,
##                    save as Tables/ANNOTATION_filled.csv
##   3. PHASE 2  : import the filled CSV, drop flagged clusters, save annotated.
##   4. PHASE 3  : labeled plots.
##
## PHASE 2 reloads the PHASE-1 checkpoint at its start, so it is deterministic
## and safe to re-run: a failed run can never leave you with a half-mutated
## object. Each block stops early with a clear message instead of burying the
## logic inside a giant if/else.
###############################################################################


## ============================================================================
## SETUP  — RUN FIRST (defines everything PHASE 1/2/3 rely on)
## ============================================================================
library(Seurat)
library(SeuratExtend)
library(scCustomize)    # optional — feature plots fall back to Seurat if absent
library(tidyverse)
library(patchwork)
library(viridis)
library(pheatmap)       # optional — publication heatmap
library(qs2)

base_dir   <- "/home/akshay-iyer/Documents/BCG_Stim_Infants_12_44_wks_PARSE_GEX"
saved_dir  <- file.path(base_dir, "saved_R_data")

cluster_col        <- "mnn.snn.louvianmlr_1"
umap_reduction     <- "umap.mnn.rna"
condition_col      <- "stim"
sample_col         <- "sample_id"
tp_col             <- "timepoint"
plot_group         <- "cell_type"   # "cell_type" (merge same-named clusters) or "cluster_celltype"
run_findallmarkers <- FALSE          # slow; OFF. If TRUE, install 'presto' for ~50x speedup.
run_pergene_plots  <- TRUE

integrated_qs   <- file.path(saved_dir, "BCG_PARSE_integrated_allMethods.qs2")
preannot_qs     <- file.path(saved_dir, "BCG_PARSE_PreAnnotation.qs2")
annotated_qs    <- file.path(saved_dir, "BCG_PARSE_integrated_annotated.qs2")

annot_dir    <- file.path(base_dir, "Annotation")
tables_dir   <- file.path(annot_dir, "Tables")
cluster_dir  <- file.path(annot_dir, "Cluster_Plots")
vln_dir      <- file.path(annot_dir, "Violin_byCluster")
vln_cond_dir <- file.path(annot_dir, "Violin_byStim")
feature_dir  <- file.path(annot_dir, "Feature_Plots")
for (d in c(tables_dir, cluster_dir, vln_dir, vln_cond_dir, feature_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- marker panel (grouped) + gene lists used by tables and plots ----
marker_panel <- list(
  Pan_T          = c("CD3D","CD3E","TRBC1","TRAC"),
  Lineage        = c("CD4","CD8A","CD8B"),
  Naive_Memory   = c("TCF7","SELL","CD27","NELL2","MAL"),
  Effector_Cyto  = c("GZMK","GNLY","NKG7","PRF1","KLRD1","CTSW","KLRB1"),
  Activation_Exh = c("PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","ENTPD1","TOX","ICOS","TNFRSF9"),
  Treg           = c("FOXP3","IL2RA","IKZF2"),
  gdT            = c("TRDC","TRGC1","TRGC2","TRDV1","TRDV2","TRGV9"),
  NK             = c("NCAM1","KLRC1","XCL1","XCL2","FCGR3A"),
  B              = c("CD19","MS4A1","CD79A","CR2","FCRL5","TCF4"),
  Plasma         = c("JCHAIN","XBP1","PRDM1","IGHG1","IGHA1","IGHM"),
  Mono_Mac       = c("CD14","SERPINA1","IFI30","LGALS1","ITGAX"),
  DC             = c("CD1C","CLEC9A","IRF8"),
  Platelet       = c("PPBP","PF4","GP9","GP1BA"),
  Lineage_TFs    = c("TBX21","GATA3","RORC","BCL6","MAF","IRF4","ZEB2","ZBTB16"),
  Cytokines      = c("IFNG","TNF","IL10","IL17A","IL21","IL1B","CXCL8","LTB"),
  Proliferation  = c("MKI67","HIST1H4C")
)

provided_genes <- c(
  "ASCL2","BATF","BATF3","BCL6","C1QBP","CCL2","CCL3","CCL4L2","CCL5","CCND3","CD14","CD19","CD1C",
  "CD200","CD27","CD3D","CD3E","CD36","CD4","CD40","CD40LG","CD70","CD7","CD79A","CD8A","CD8B",
  "CLEC9A","CR2","CTLA4","CTSW","CXCL8","CXCR3","CXCR5","EBI3","ENTPD1","FABP5","FCGR2B","FCGR3A",
  "FCRL5","FOXP3","FASLG","GNLY","GP1BA","GP9","GATA3","GZMK","HAVCR2","HIF1A","HIST1H4C","HLA-DPA1",
  "HLA-DRA","HLA-DRB1","ICOS","IFI30","IFNG","IGFBP2","IGFBP4","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3",
  "IGHG4","IGHM","IKZF2","IL10","IL17A","IL18BP","IL18RAP","IL1B","IL21","IL2RA","IL2RB","IRF4","IRF8",
  "ITGAX","JCHAIN","KLRB1","KLRD1","KLRC1","LAG3","LDHA","LGALS1","LTA","LTB","MAF","MAL","MALAT1",
  "MIR155HG","MKI67","MT-ND1","MT-ND5","MS4A1","NELL2","NCAM1","NKG7","NR4A1","PDCD1","PF4","PPBP",
  "PRDM1","PRF1","RORC","SELL","SERPINA1","SERPING1","SH2D1A","TCF4","TCF7","TIGIT","TNF","TNFAIP2",
  "TNFRSF18","TNFRSF4","TNFRSF9","TOX","TBX21","TRBC1","TRDC","TRDV1","TRDV2","TRGC1","TRGC2","TRGV9",
  "XBP1","XCL1","XCL2","ZBTB16","ZEB2"
)

lineage_sigs <- list(
  T        = c("CD3D","CD3E","CD3G","TRBC1","TRAC"),
  gdT      = c("TRDC","TRGC1","TRGC2","TRDV1","TRDV2"),
  B        = c("CD19","MS4A1","CD79A","CD79B","CR2"),
  Plasma   = c("JCHAIN","XBP1","PRDM1","MZB1","IGHG1"),
  NK       = c("GNLY","NKG7","KLRD1","KLRC1","NCAM1","PRF1"),
  Monocyte = c("CD14","FCGR3A","LYZ","SERPINA1","IFI30","S100A8","S100A9"),
  DC       = c("CD1C","CLEC9A","IRF8","FCER1A"),
  Platelet = c("PPBP","PF4","GP9","GP1BA")
)

# helper: commit cluster identity + numeric level ordering on an object
prep_clusters <- function(obj) {
  obj@meta.data[[cluster_col]] <- droplevels(factor(obj@meta.data[[cluster_col]]))
  lv <- levels(obj@meta.data[[cluster_col]])
  nlv <- suppressWarnings(as.numeric(lv)); if (!any(is.na(nlv))) lv <- lv[order(nlv)]
  obj@meta.data[[cluster_col]] <- factor(obj@meta.data[[cluster_col]], levels = lv)
  Idents(obj) <- cluster_col
  obj
}
message("SETUP done.")


## ============================================================================
## PHASE 1 — TABLES   (run, then fill the template and save ANNOTATION_filled.csv)
## ============================================================================
seu <- qs_read(integrated_qs)
DefaultAssay(seu) <- "RNA"
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
seu <- prep_clusters(seu)

# panels restricted to genes present (depends on `seu`, so computed here)
marker_panel  <- lapply(marker_panel, intersect, y = rownames(seu[["RNA"]]))
marker_panel  <- marker_panel[vapply(marker_panel, length, integer(1)) > 0]
rna.features  <- unique(unlist(marker_panel, use.names = FALSE))
gene2group    <- unlist(lapply(names(marker_panel),
                               function(g) setNames(rep(g, length(marker_panel[[g]])), marker_panel[[g]])))
gene2group    <- gene2group[!duplicated(names(gene2group))]
heat_feats    <- names(gene2group)
pergene_genes <- intersect(unique(c(provided_genes, rna.features)), rownames(seu[["RNA"]]))
lineage_sigs  <- lapply(lineage_sigs, intersect, y = rownames(seu[["RNA"]]))
lineage_sigs  <- lineage_sigs[vapply(lineage_sigs, length, integer(1)) >= 2]

# lineage module scores -> dominant lineage
seu <- AddModuleScore(seu, features = lineage_sigs, name = "LINSCORE_", assay = "RNA")
score_cols <- paste0("score_", names(lineage_sigs))
colnames(seu@meta.data)[match(paste0("LINSCORE_", seq_along(lineage_sigs)), colnames(seu@meta.data))] <- score_cols
seu$lineage_score_call <- factor(
  names(lineage_sigs)[max.col(as.matrix(seu@meta.data[, score_cols, drop = FALSE]), ties.method = "first")],
  levels = names(lineage_sigs))

cl_score <- seu@meta.data %>% group_by(cluster = as.character(.data[[cluster_col]])) %>%
  summarise(across(all_of(score_cols), ~ round(mean(.x), 4)), .groups = "drop")
write.csv(cl_score, file.path(tables_dir, "lineage_module_scores_byCluster.csv"), row.names = FALSE)

lin_call <- seu@meta.data %>%
  count(cluster = as.character(.data[[cluster_col]]), lineage_score_call) %>%
  group_by(cluster) %>% mutate(frac = n / sum(n)) %>%
  slice_max(frac, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(cluster, dominant_lineage = as.character(lineage_score_call), dominant_frac = round(frac, 3))
write.csv(lin_call, file.path(tables_dir, "cluster_dominant_lineage.csv"), row.names = FALSE)

size_tbl <- seu@meta.data %>% group_by(cluster = as.character(.data[[cluster_col]])) %>%
  summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  mutate(pct_of_total = round(100 * n_cells / sum(n_cells), 2))
write.csv(size_tbl, file.path(tables_dir, "cluster_sizes.csv"), row.names = FALSE)

comp_tbl <- function(col) seu@meta.data %>%
  count(cluster = as.character(.data[[cluster_col]]), grp = as.character(.data[[col]])) %>%
  pivot_wider(names_from = grp, values_from = n, values_fill = 0)
write.csv(comp_tbl(condition_col), file.path(tables_dir, "cluster_counts_by_stim.csv"), row.names = FALSE)
write.csv(comp_tbl(tp_col),        file.path(tables_dir, "cluster_counts_by_timepoint.csv"), row.names = FALSE)
write.csv(comp_tbl(sample_col),    file.path(tables_dir, "cluster_counts_by_sample.csv"), row.names = FALSE)

dp <- Seurat::DotPlot(seu, features = pergene_genes, group.by = cluster_col)$data %>%
  rename(gene = features.plot, cluster = id)
write.csv(dp, file.path(tables_dir, "MarkerPanel_dotplot_long.csv"), row.names = FALSE)

avg_all <- AverageExpression(seu, assays = "RNA", features = pergene_genes, group.by = cluster_col, layer = "data")$RNA
colnames(avg_all) <- sub("^g", "", colnames(avg_all))
np <- suppressWarnings(as.numeric(colnames(avg_all))); if (!any(is.na(np))) avg_all <- avg_all[, order(np), drop = FALSE]
write.csv(round(as.matrix(avg_all), 4), file.path(tables_dir, "avg_expression_byCluster.csv"))

if (run_findallmarkers) {
  all_markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                                test.use = "wilcox", max.cells.per.ident = 500)
  write.csv(all_markers, file.path(tables_dir, "FindAllMarkers.csv"), row.names = FALSE)
}

top_str <- dp %>% filter(pct.exp >= 10) %>% group_by(cluster) %>%
  slice_max(order_by = avg.exp.scaled, n = 15, with_ties = FALSE) %>%
  summarise(top_markers = paste(gene, collapse = ", "), .groups = "drop") %>%
  mutate(cluster = as.character(cluster))
template <- size_tbl %>% left_join(lin_call, by = "cluster") %>% left_join(top_str, by = "cluster") %>%
  arrange(suppressWarnings(as.numeric(cluster)), cluster) %>%
  mutate(cell_type = "", lineage = "", state = "", confidence = "", flag = "", notes = "")
write.csv(template, file.path(tables_dir, "ANNOTATION_TEMPLATE_fill_me.csv"), row.names = FALSE)

qs_save(seu, preannot_qs)
message("PHASE 1 done. Tables in ", tables_dir,
        "\n  -> fill ANNOTATION_TEMPLATE_fill_me.csv, save as ANNOTATION_filled.csv, then run PHASE 2.")


## ============================================================================
## PHASE 2 — IMPORT ANNOTATION   (reloads PHASE-1 checkpoint; flat; re-runnable)
## ============================================================================
annot_file <- file.path(tables_dir, "ANNOTATION_filled.csv")
if (!file.exists(annot_file))
  stop("Not found: ", annot_file,
       "\n  Fill ANNOTATION_TEMPLATE_fill_me.csv and save it as ANNOTATION_filled.csv first.")

seu <- qs_read(preannot_qs)          # deterministic clean start
seu <- prep_clusters(seu)

annot <- readr::read_csv(annot_file, show_col_types = FALSE) %>% mutate(cluster = as.character(cluster))

missing <- setdiff(as.character(levels(seu@meta.data[[cluster_col]])), annot$cluster)
if (length(missing)) stop("Clusters missing from annotation CSV: ", paste(missing, collapse = ", "))

drop_flags <- c("doublet","lowqual","low_quality","remove","junk","drop")
to_drop    <- annot$cluster[tolower(trimws(annot$flag)) %in% drop_flags]
blank_lab  <- setdiff(annot$cluster[is.na(annot$cell_type) | trimws(annot$cell_type) == ""], to_drop)
if (length(blank_lab)) stop("Blank cell_type with no removal flag: ", paste(blank_lab, collapse = ", "))

cl_chr <- as.character(seu@meta.data[[cluster_col]])
m <- match(cl_chr, annot$cluster)
seu$cell_type        <- annot$cell_type[m]
seu$lineage          <- annot$lineage[m]
seu$state            <- annot$state[m]
seu$cluster_celltype <- paste0(cl_chr, ": ", seu$cell_type)

if (length(to_drop)) {
  message("Dropping flagged clusters: ", paste(to_drop, collapse = ", "))
  seu <- subset(seu, cells = rownames(seu@meta.data)[!cl_chr %in% to_drop])
  seu@meta.data[[cluster_col]] <- droplevels(factor(seu@meta.data[[cluster_col]]))
  annot  <- annot[!annot$cluster %in% to_drop, ]
  cl_chr <- as.character(seu@meta.data[[cluster_col]])
}

lin_rank <- c(T = 1, gdT = 2, NK = 3, B = 4, Plasma = 5, Myeloid = 6, DC = 7, Unresolved = 8)
ord <- annot %>% mutate(cn = suppressWarnings(as.numeric(cluster)),
                        lr = ifelse(lineage %in% names(lin_rank), lin_rank[lineage], 99)) %>%
  arrange(lr, cn)
seu$cell_type        <- factor(seu$cell_type,        levels = unique(ord$cell_type))
seu$cluster_celltype <- factor(seu$cluster_celltype, levels = paste0(ord$cluster, ": ", ord$cell_type))
seu$lineage          <- factor(seu$lineage,          levels = intersect(names(lin_rank), unique(ord$lineage)))
Idents(seu) <- plot_group

qs_save(seu, annotated_qs)
message("PHASE 2 done -> ", basename(annotated_qs), ". Cells per cell_type:")
print(table(seu$cell_type))


## ============================================================================
## PHASE 3 — PLOTS   (needs SETUP run; uses `seu` from PHASE 2)
## ============================================================================
if (!"cell_type" %in% colnames(seu@meta.data))
  stop("No 'cell_type' in object — run PHASE 2 first.")

# annotated UMAPs
ggsave(file.path(cluster_dir, "UMAP_annotation.png"),
       DimPlot2(seu, reduction = umap_reduction, group.by = plot_group,
                label = TRUE, repel = TRUE, label.size = 3.5) + ggtitle("Annotation"),
       width = 14, height = 9, dpi = 300, bg = "white")
ggsave(file.path(cluster_dir, "UMAP_annotation_vs_stim_timepoint.png"),
       wrap_plots(
         DimPlot2(seu, reduction = umap_reduction, group.by = plot_group,
                  label = TRUE, repel = TRUE, label.size = 3) + ggtitle("Cell type"),
         DimPlot2(seu, reduction = umap_reduction, group.by = condition_col) + ggtitle("Stim"),
         DimPlot2(seu, reduction = umap_reduction, group.by = tp_col) + ggtitle("Timepoint"),
         ncol = 3),
       width = 26, height = 8, dpi = 300, bg = "white")

# composition by cell type
distr <- function(group_vec, fname, w = 16)
  ggsave(file.path(cluster_dir, fname),
         ClusterDistrBar(group_vec, seu@meta.data[[plot_group]], flip = FALSE, border = "black") +
           theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)),
         width = w, height = 9, dpi = 300, bg = "white")
distr(seu@meta.data[[sample_col]],    "Composition_by_sample.png")
distr(seu@meta.data[[condition_col]], "Composition_by_stim.png")
distr(seu@meta.data[[tp_col]],        "Composition_by_timepoint.png")
if ("IGRA_status" %in% colnames(seu@meta.data))
  distr(seu@meta.data[["IGRA_status"]], "Composition_by_IGRA.png")

# grouped summary DotPlot
ggsave(file.path(cluster_dir, "DotPlot_marker_panel.png"),
       DotPlot2(seu, features = marker_panel, group.by = plot_group) +
         theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)),
       width = 16, height = 20, dpi = 300, bg = "white", limitsize = FALSE)

# publication heatmap (mean expr by cell type, row z-score)
avg_ct <- AverageExpression(seu, assays = "RNA", features = heat_feats, group.by = plot_group, layer = "data")$RNA
avg_ct <- log1p(as.matrix(avg_ct))[heat_feats, , drop = FALSE]
group_colors <- c(
  Pan_T="#5E60CE", Lineage="#6930C3", Naive_Memory="#74C2E1", Effector_Cyto="#EF5350",
  Activation_Exh="#A1887F", Treg="#F06292", gdT="#26C6DA", NK="#9C6FD6",
  B="#42A5F5", Plasma="#1E88E5", Mono_Mac="#FFB300", DC="#FB8C00",
  Platelet="#8D6E63", Lineage_TFs="#66BB6A", Cytokines="#26A69A", Proliferation="#FF8A65")
group_colors <- group_colors[names(group_colors) %in% gene2group]
row_annot <- data.frame(Group = factor(gene2group[heat_feats], levels = names(group_colors)), row.names = heat_feats)
gaps <- which(head(as.character(row_annot$Group), -1) != tail(as.character(row_annot$Group), -1))
png(file.path(cluster_dir, "Heatmap_marker_panel.png"), width = 13, height = 20, units = "in", res = 300, bg = "white")
pheatmap(avg_ct, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_row = row_annot, annotation_colors = list(Group = group_colors), gaps_row = gaps,
         color = colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(100),
         border_color = NA, fontsize_row = 9, fontsize_col = 11,
         main = "Mean expression (row z-score) by cell type")
dev.off()

# per-gene Vln + Feature plots (grouped by cell type)
if (run_pergene_plots) {
  pal <- viridis(n = 10, option = "A")
  have_sccustom <- requireNamespace("scCustomize", quietly = TRUE)
  if (!have_sccustom) message("scCustomize not installed -> using Seurat::FeaturePlot.")
  feature_plot <- function(gene) {
    if (have_sccustom)
      scCustomize::FeaturePlot_scCustom(seu, reduction = umap_reduction, features = gene, colors_use = pal, order = TRUE)
    else
      Seurat::FeaturePlot(seu, reduction = umap_reduction, features = gene, order = TRUE) +
      viridis::scale_color_viridis(option = "A")
  }
  message("Rendering per-gene plots for ", length(pergene_genes), " genes ...")
  for (g in pergene_genes) {
    f1 <- file.path(vln_dir,      paste0(g, "_Vln.png"))
    f2 <- file.path(vln_cond_dir, paste0(g, "_Vln_byStim.png"))
    f3 <- file.path(feature_dir,  paste0(g, "_Feature.png"))
    if (!file.exists(f1))
      ggsave(f1, VlnPlot2(seu, features = g, group.by = plot_group, cols = "default", show.mean = TRUE) +
               ggtitle(paste("RNA |", g)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
             dpi = 300, width = 14, height = 8, bg = "white")
    if (!file.exists(f2))
      ggsave(f2, VlnPlot2(seu, features = g, group.by = plot_group, cols = "default",
                          split.by = condition_col, stat.method = "wilcox.test") +
               ggtitle(paste("RNA |", g, "| split by stim")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
             dpi = 300, width = 14, height = 8, bg = "white")
    if (!file.exists(f3))
      ggsave(f3, feature_plot(g), dpi = 300, width = 8, height = 7, bg = "white")
  }
  message("Per-gene plots done.")
}
message("PHASE 3 done — labeled plots under ", annot_dir)
###############################################################################

