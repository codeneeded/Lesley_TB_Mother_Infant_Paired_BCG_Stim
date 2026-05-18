# ===============================
# TB_ALL: DGE/DPE by IGRA status
# Per-cluster MAST tests for:
#   - Infant IGRA status (IGRA_infant)
#   - Mother IGRA status (IGRA_mother)
# Both comparing IGRA+ vs IGRA-
# ===============================

library(Seurat)
library(dplyr)
library(ggplot2)
library(qs2)
library(SeuratExtend)

# -------------------------------------------------
# 0) Load annotated, cleaned object
#    (clusters 14, 23, 24 already removed)
# -------------------------------------------------

base_dir   <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/TB_CITEseq"
load.path  <- file.path(base_dir, "saved_R_data")

TB_ALL <- qs_read(file.path(load.path, "TB_ALL_Annotated_Clean.qs2"))

# ---- Join split layers (Seurat v5 requirement for FindMarkers) ----
message("Joining RNA layers...")
TB_ALL[["RNA"]] <- JoinLayers(TB_ALL[["RNA"]])

if ("ADT" %in% Assays(TB_ALL)) {
  message("Joining ADT layers...")
  TB_ALL[["ADT"]] <- JoinLayers(TB_ALL[["ADT"]])
}

DefaultAssay(TB_ALL) <- "RNA"
Idents(TB_ALL)      <- "celltype_fine"

DefaultAssay(TB_ALL) <- "RNA"
Idents(TB_ALL)      <- "celltype_fine"   # change to "celltype_main" if you prefer collapsed labels

# Sanity check — confirm IGRA columns exist and look right
cat("\n--- IGRA_infant ---\n")
print(table(TB_ALL$IGRA_infant, useNA = "ifany"))
cat("\n--- IGRA_mother ---\n")
print(table(TB_ALL$IGRA_mother, useNA = "ifany"))
cat("\n--- IGRA_infant x IGRA_mother (cell-level) ---\n")
print(table(TB_ALL$IGRA_infant, TB_ALL$IGRA_mother, useNA = "ifany"))

# ======================================
# 1) Optional: cluster distribution plots
#    by IGRA status (sample-level boxplots)
# ======================================

post_annot_dir <- file.path(base_dir, "Annotation", "Post-Annotation")
dir.create(post_annot_dir, recursive = TRUE, showWarnings = FALSE)

# Make IGRA factors with a consistent order so OUD+ / IGRA+ is the "positive" reference
TB_ALL$IGRA_infant <- factor(TB_ALL$IGRA_infant, levels = c("IGRA+", "IGRA-"))
TB_ALL$IGRA_mother <- factor(TB_ALL$IGRA_mother, levels = c("IGRA+", "IGRA-"))

igra_cols <- c("IGRA+" = "#F54927", "IGRA-" = "#1180808B")

# -- Infant IGRA
ClusterDistrPlot(
  origin      = TB_ALL$Sample_ID,
  cluster     = TB_ALL$celltype_fine,
  condition   = TB_ALL$IGRA_infant,
  stat.method = "wilcox.test",
  cols        = igra_cols
)
ggsave(
  filename = file.path(post_annot_dir, "TB_ClusterDistr_IGRA_infant_BOXPLOT.png"),
  width = 14, height = 9, dpi = 300, bg = "white"
)

p_inf <- ClusterDistrPlot(
  origin = TB_ALL$Sample_ID, cluster = TB_ALL$celltype_fine,
  condition = TB_ALL$IGRA_infant, stat.method = "wilcox.test", plot = FALSE
)
write.csv(as.data.frame(p_inf),
          file.path(post_annot_dir, "TB_ClusterDistr_IGRA_infant.csv"),
          row.names = TRUE)

# -- Mother IGRA
ClusterDistrPlot(
  origin      = TB_ALL$Sample_ID,
  cluster     = TB_ALL$celltype_fine,
  condition   = TB_ALL$IGRA_mother,
  stat.method = "wilcox.test",
  cols        = igra_cols
)
ggsave(
  filename = file.path(post_annot_dir, "TB_ClusterDistr_IGRA_mother_BOXPLOT.png"),
  width = 14, height = 9, dpi = 300, bg = "white"
)

p_mom <- ClusterDistrPlot(
  origin = TB_ALL$Sample_ID, cluster = TB_ALL$celltype_fine,
  condition = TB_ALL$IGRA_mother, stat.method = "wilcox.test", plot = FALSE
)
write.csv(as.data.frame(p_mom),
          file.path(post_annot_dir, "TB_ClusterDistr_IGRA_mother.csv"),
          row.names = TRUE)

# ======================================
# 2) Helper: cluster-wise DE function
#    (same shape as the OPIS pipeline)
# ======================================

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
count_sig <- function(df, pcol = "p_val_adj", thr = 0.05) sum(df[[pcol]] < thr, na.rm = TRUE)

run_clusterwise_de <- function(obj,
                               assay,
                               group_col,
                               ident1,
                               ident2,
                               out_dir,
                               latent  = c("nCount_RNA"),
                               min_cells_per_grp = 10,
                               title_prefix = "",
                               file_stub = NULL,
                               annotation_col = "celltype_fine") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  DefaultAssay(obj) <- assay
  Idents(obj) <- annotation_col
  
  de_list   <- list()
  count_tbl <- data.frame()
  
  cl_levels <- levels(obj)
  
  for (cl in cl_levels) {
    message("[", assay, "] Processing cluster: ", cl)
    
    cl_cells <- WhichCells(obj, idents = cl)
    sc <- subset(obj, cells = cl_cells)
    sc$.__grp <- sc[[group_col]]
    
    tab <- table(sc$.__grp)
    
    # Require both groups present and enough cells in each
    if (all(c(ident1, ident2) %in% names(tab)) &&
        all(tab[c(ident1, ident2)] >= min_cells_per_grp)) {
      
      de <- tryCatch(
        FindMarkers(
          sc,
          ident.1     = ident1,
          ident.2     = ident2,
          group.by    = ".__grp",
          test.use    = "MAST",
          latent.vars = latent
        ),
        error = function(e) {
          message("  FindMarkers failed for ", cl, ": ", conditionMessage(e))
          NULL
        }
      )
      
      if (!is.null(de)) {
        comp_name <- if (is.null(file_stub)) paste0(ident1, "_vs_", ident2) else file_stub
        out_csv   <- file.path(out_dir, paste0(safe_name(cl), "_", comp_name, ".csv"))
        
        write.csv(de, out_csv)
        de_list[[cl]] <- de
        
        count_tbl <- rbind(
          count_tbl,
          data.frame(Cluster = cl, DE_Genes = count_sig(de))
        )
      }
    } else {
      message("  Skipping ", cl,
              " | counts: ", paste(names(tab), tab, sep = "=", collapse = ", "))
    }
  }
  
  # Barplot of DE counts per cluster
  if (nrow(count_tbl) > 0) {
    p <- ggplot(count_tbl, aes(x = reorder(Cluster, -DE_Genes), y = DE_Genes, fill = Cluster)) +
      geom_bar(stat = "identity") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      xlab("Cluster") +
      ylab("# of DE Features (adj. p < 0.05)") +
      ggtitle(paste0(title_prefix, ": Clusters Ranked by DE Features (",
                     ident1, " vs ", ident2, " | ", assay, ")"))
    
    stub <- if (is.null(file_stub)) paste0(ident1, "vs", ident2) else file_stub
    
    ggsave(
      filename = file.path(out_dir, paste0("Cluster_DE_Counts_", stub, "_", assay, ".png")),
      plot   = p,
      width  = 11,
      height = 6,
      dpi    = 400,
      bg     = "white"
    )
    
    write.csv(
      count_tbl,
      file.path(out_dir, paste0("DE_Counts_Summary_", stub, "_", assay, ".csv")),
      row.names = FALSE
    )
  }
  
  invisible(list(de = de_list, counts = count_tbl))
}

# ======================================
# 3) DGE (RNA) + DPE (ADT) by IGRA_infant
# ======================================

de_base <- file.path(base_dir, "Differential_Expression")

dge_infant_dir <- file.path(de_base, "DGE", "IGRA_infant_Pos_vs_Neg")
dpe_infant_dir <- file.path(de_base, "DPE", "IGRA_infant_Pos_vs_Neg")

TB_infant <- subset(TB_ALL, subset = !is.na(IGRA_infant))

cat("\n--- IGRA_infant counts (post-NA filter) ---\n")
print(table(TB_infant$IGRA_infant))

dge_infant <- run_clusterwise_de(
  obj          = TB_infant,
  assay        = "RNA",
  group_col    = "IGRA_infant",
  ident1       = "IGRA+",
  ident2       = "IGRA-",
  out_dir      = dge_infant_dir,
  latent       = c("nCount_RNA"),
  title_prefix = "TB IGRA_infant (RNA)",
  file_stub    = "IGRA_infant_Pos_vs_Neg_RNA"
)

dpe_infant <- run_clusterwise_de(
  obj          = TB_infant,
  assay        = "ADT",
  group_col    = "IGRA_infant",
  ident1       = "IGRA+",
  ident2       = "IGRA-",
  out_dir      = dpe_infant_dir,
  latent       = c("nCount_RNA", "nCount_ADT"),
  title_prefix = "TB IGRA_infant (ADT)",
  file_stub    = "IGRA_infant_Pos_vs_Neg_ADT"
)

# ======================================
# 4) DGE (RNA) + DPE (ADT) by IGRA_mother
# ======================================

dge_mother_dir <- file.path(de_base, "DGE", "IGRA_mother_Pos_vs_Neg")
dpe_mother_dir <- file.path(de_base, "DPE", "IGRA_mother_Pos_vs_Neg")

TB_mother <- subset(TB_ALL, subset = !is.na(IGRA_mother))

cat("\n--- IGRA_mother counts (post-NA filter) ---\n")
print(table(TB_mother$IGRA_mother))

dge_mother <- run_clusterwise_de(
  obj          = TB_mother,
  assay        = "RNA",
  group_col    = "IGRA_mother",
  ident1       = "IGRA+",
  ident2       = "IGRA-",
  out_dir      = dge_mother_dir,
  latent       = c("nCount_RNA"),
  title_prefix = "TB IGRA_mother (RNA)",
  file_stub    = "IGRA_mother_Pos_vs_Neg_RNA"
)

dpe_mother <- run_clusterwise_de(
  obj          = TB_mother,
  assay        = "ADT",
  group_col    = "IGRA_mother",
  ident1       = "IGRA+",
  ident2       = "IGRA-",
  out_dir      = dpe_mother_dir,
  latent       = c("nCount_RNA", "nCount_ADT"),
  title_prefix = "TB IGRA_mother (ADT)",
  file_stub    = "IGRA_mother_Pos_vs_Neg_ADT"
)

message("\n=== DGE/DPE pipeline complete ===")
message("Results saved under: ", de_base)

