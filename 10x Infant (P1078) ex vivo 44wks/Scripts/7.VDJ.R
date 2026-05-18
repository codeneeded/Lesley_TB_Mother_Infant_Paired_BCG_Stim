###############################################################################
#  TB CITE-seq Pipeline — Script 7: TCR Repertoire Analysis
#  Adapted from TARA/EARTH TCR pipeline using scRepertoire
###############################################################################

library(Seurat)
library(tidyverse)
library(ggplot2)
library(scRepertoire)
library(igraph)
library(Cairo)
library(RColorBrewer)
library(Polychrome)
library(patchwork)
library(Trex)
library(qs2)

###############################################################################
#                         PATHS & SETUP
###############################################################################

base_dir  <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/TB_CITEseq"
load.path <- file.path(base_dir, "saved_R_data")
in.path   <- "/media/akshay-iyer/Elements/10x_TB_organized/Cell_Ranger_Out/"

# VDJ output directories
vdj_dir       <- file.path(base_dir, "VDJ", "TCR")
clonal_viz    <- file.path(vdj_dir, "Clonal_Visualizations")
cd3_comp      <- file.path(vdj_dir, "CD3_Composition")
clonal_div    <- file.path(vdj_dir, "Clonal_Diversity")
network_dir   <- file.path(vdj_dir, "Network_Analysis")
seurat_plots  <- file.path(vdj_dir, "Seurat_Plots")
hyperexp_dir  <- file.path(seurat_plots, "Hyperexpansion")
occupy_dir    <- file.path(seurat_plots, "Clones_per_Cluster")
overlay_dir   <- file.path(seurat_plots, "Clonal_Overlay")
trex_dir      <- file.path(vdj_dir, "Trex")
trex_pie_dir  <- file.path(trex_dir, "PID_Pie_Chart")

for (d in c(clonal_viz, cd3_comp, clonal_div, network_dir,
            hyperexp_dir, occupy_dir, overlay_dir, trex_dir, trex_pie_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

###############################################################################
#  SAMPLE METADATA
###############################################################################

# Sample names (Cell Ranger multi folder names)
f_names <- c("Bone-19321-001", "Bone-19321-002", "Bone-19321-003", "Bone-19321-004",
             "Bone-19430-001", "Bone-19430-002", "Bone-19429-001", "Bone-19429-002",
             "Bone-19649-001", "Bone-19708-001")

# Sample ID mapping
sample_meta <- data.frame(
  orig.ident  = f_names,
  Sample_ID   = c("6062151", "6067021", "6063421", "6063411",
                  "6061161", "6063481", "6044931", "6044981",
                  "6071261", "6065761"),
  IGRA_infant = c("+", "+", "-", "-", "+", "-", "-", "-", "+", "+"),
  IGRA_mother = c("+", "-", "+", "-", "+", "+", "-", "+", "+", "+"),
  stringsAsFactors = FALSE
)

###############################################################################
#  STEP 1: LOAD FILTERED CONTIG ANNOTATIONS
###############################################################################

message("Loading TCR contig annotations...")

for (name in f_names) {
  t_file <- file.path(in.path, name, "outs", "per_sample_outs", name,
                      "vdj_t", "filtered_contig_annotations.csv")
  
  if (!file.exists(t_file)) {
    warning("TCR file not found for: ", name, " — skipping")
    next
  }
  
  t_data <- read.csv(t_file)
  assign(paste0(name, ".TCR"), t_data)
  message("  Loaded: ", name, " (", nrow(t_data), " contigs)")
}

###############################################################################
#  STEP 2: CREATE CONTIG LIST & COMBINE TCR
###############################################################################

TB.TCR.names   <- paste0(f_names, ".TCR")
TB.contig_list <- as.list(mget(TB.TCR.names))

# Combine using S_ prefixed Sample_ID to ensure character (not numeric)
sample_labels <- paste0("S_", sample_meta$Sample_ID)

combined.TCR.TB <- combineTCR(
  TB.contig_list,
  samples = sample_labels
)

# Add IGRA status variables
combined.TCR.TB <- addVariable(
  combined.TCR.TB,
  variable.name = "IGRA_infant",
  variables     = ifelse(sample_meta$IGRA_infant == "+", "IGRA+", "IGRA-")
)

combined.TCR.TB <- addVariable(
  combined.TCR.TB,
  variable.name = "IGRA_mother",
  variables     = ifelse(sample_meta$IGRA_mother == "+", "IGRA+", "IGRA-")
)

combined.TCR.TB <- addVariable(
  combined.TCR.TB,
  variable.name = "IGRA_group",
  variables     = paste0("Inf:", ifelse(sample_meta$IGRA_infant == "+", "+", "-"),
                         "_Mat:", ifelse(sample_meta$IGRA_mother == "+", "+", "-"))
)

message("Combined TCR for ", length(combined.TCR.TB), " samples")

###############################################################################
#  STEP 3: CLONAL VISUALIZATIONS
###############################################################################

setwd(clonal_viz)

# --- Unique Clones ---
clonalQuant(combined.TCR.TB, cloneCall = "strict", chain = "both", scale = FALSE)
ggsave("TB_Unique_Clones_Strict_raw.png", width = 18, height = 10, bg = "white")

clonalQuant(combined.TCR.TB, cloneCall = "strict", chain = "both", scale = TRUE)
ggsave("TB_Unique_Clones_Strict_scaled.png", width = 18, height = 10, bg = "white")

# By IGRA infant
clonalQuant(combined.TCR.TB, cloneCall = "strict", chain = "both",
            group.by = "IGRA_infant", scale = FALSE)
ggsave("TB_Unique_Clones_byIGRA_infant.png", width = 18, height = 10, bg = "white")

# By IGRA mother
clonalQuant(combined.TCR.TB, cloneCall = "strict", chain = "both",
            group.by = "IGRA_mother", scale = FALSE)
ggsave("TB_Unique_Clones_byIGRA_mother.png", width = 18, height = 10, bg = "white")

# --- Clonal Abundance ---
clonalAbundance(combined.TCR.TB, cloneCall = "strict", scale = FALSE)
ggsave("TB_Clonal_Abundance_raw.png", width = 16, height = 12, bg = "white")

clonalAbundance(combined.TCR.TB, cloneCall = "strict", scale = TRUE)
ggsave("TB_Clonal_Abundance_scaled.png", width = 16, height = 12, bg = "white")

# --- Clonal Length ---
clonalLength(combined.TCR.TB, cloneCall = "aa", chain = "both")
ggsave("TB_Clonal_Length_raw.png", width = 16, height = 12, bg = "white")

clonalLength(combined.TCR.TB, cloneCall = "aa", chain = "both", scale = TRUE)
ggsave("TB_Clonal_Length_scaled.png", width = 16, height = 12, bg = "white")

# --- Clonal Homeostasis ---
clonalHomeostasis(combined.TCR.TB, cloneCall = "strict",
                  cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1,
                                Large = 0.3, Hyperexpanded = 1))
ggsave("TB_Clonal_Homeostasis.png", width = 20, height = 12, bg = "white")

###############################################################################
#  STEP 4: CDR3 COMPOSITION
###############################################################################

setwd(cd3_comp)

percentAA(combined.TCR.TB, chain = "TRA", aa.length = 20)
ggsave("TB_Percent_AA_TRA.png", width = 26, height = 24, bg = "white")

percentAA(combined.TCR.TB, chain = "TRB", aa.length = 20)
ggsave("TB_Percent_AA_TRB.png", width = 26, height = 24, bg = "white")

positionalEntropy(combined.TCR.TB, chain = "both", aa.length = 20)
ggsave("TB_Positional_Entropy.png", width = 20, height = 15, bg = "white")

# V/J Gene Usage Heatmaps
for (gene in c("TRAV", "TRBV", "TRAJ", "TRBJ")) {
  vizGenes(combined.TCR.TB, x.axis = gene, y.axis = NULL,
           plot = "heatmap", scale = TRUE)
  ggsave(paste0("TB_Heatmap_", gene, ".png"), width = 12, height = 9, bg = "white")
}

# K-mer Analysis
for (chain in c("TRA", "TRB")) {
  percentKmer(combined.TCR.TB, cloneCall = "aa", chain = chain,
              motif.length = 3, top.motifs = 25)
  ggsave(paste0("TB_Kmer_", chain, ".png"), width = 12, height = 9, bg = "white")
}

###############################################################################
#  STEP 5: CLONAL DIVERSITY & OVERLAP
###############################################################################

setwd(clonal_div)

clonalDiversity(combined.TCR.TB, cloneCall = "strict")
ggsave("TB_Clonal_Diversity.png", width = 14, height = 9, bg = "white")

clonalOverlap(combined.TCR.TB, cloneCall = "strict", method = "morisita")
ggsave("TB_Clonal_Overlap_morisita.png", width = 18, height = 13, bg = "white")

clonalOverlap(combined.TCR.TB, cloneCall = "strict", method = "raw")
ggsave("TB_Clonal_Overlap_raw.png", width = 18, height = 13, bg = "white")

###############################################################################
#  STEP 6: CLONAL NETWORK CLUSTERING
###############################################################################

setwd(network_dir)

custom_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#A6D854", "#FFD92F", "#E5C494"
)

set_colors_and_legend <- function(igraph.object) {
  col_legend   <- factor(igraph::V(igraph.object)$group)
  color.legend <- unique(igraph::V(igraph.object)$group)
  ordered_idx  <- order(nchar(color.legend), decreasing = TRUE)
  color.legend <- color.legend[ordered_idx]
  num_colors   <- length(color.legend)
  col_samples  <- custom_colors[seq_len(num_colors)][ordered_idx]
  list(col_legend = col_legend, color.legend = color.legend, col_samples = col_samples)
}

plot_igraph <- function(igraph.object, col_samples, color.legend, title_text) {
  plot(igraph.object,
       vertex.size     = sqrt(igraph::V(igraph.object)$size),
       vertex.label    = NA,
       edge.arrow.size = 0,
       edge.curved     = 0.3,
       vertex.color    = col_samples)
  legend("topleft", legend = color.legend, pch = 16,
         col = unique(col_samples), bty = "n")
  title(title_text, cex.main = 1.5, font.main = 2)
}

# All samples network
for (chain_info in list(
  list(chain = "TRA", threshold = 0.90),
  list(chain = "TRB", threshold = 0.85)
)) {
  ig <- clonalCluster(
    combined.TCR.TB,
    chain        = chain_info$chain,
    sequence     = "aa",
    group.by     = "sample",
    threshold    = chain_info$threshold,
    exportGraph  = TRUE
  )
  
  colors <- set_colors_and_legend(ig)
  
  CairoPNG(paste0("TB_ALL_", chain_info$chain, ".png"),
           width = 7200, height = 5500, res = 500)
  plot_igraph(ig, colors$col_samples, colors$color.legend,
              paste0("TB All Samples ", chain_info$chain,
                     " (Threshold ", chain_info$threshold, ")"))
  dev.off()
}

###############################################################################
#  STEP 7: MERGE TCR WITH SEURAT OBJECT
###############################################################################

message("Loading Seurat WNN object...")
TB_ALL <- qs_read(file.path(load.path, "TB_ALL_WNN.qs2"))

# Fix barcodes to match scRepertoire format: S_SampleID_BARCODE-1
# scRepertoire barcodes look like: S_6062151_BARCODE-1 (S_ prefix from combineTCR)
# Seurat barcodes currently: Bone-19321-001_BARCODE-1 (from merge with add.cell.id)
barcodes <- colnames(TB_ALL)

# Map orig.ident to S_ prefixed Sample_ID (matching combineTCR labels)
sid_map <- setNames(paste0("S_", sample_meta$Sample_ID), sample_meta$orig.ident)

# Extract barcode suffix (the actual cell barcode after Bone-XXXXX-XXX_)
barcode_suffix <- gsub("^Bone-[0-9]+-[0-9]+_", "", barcodes)

# Build new barcodes matching scRepertoire format: S_6062151_ACGTACGT-1
new_barcodes <- paste0(sid_map[TB_ALL$orig.ident], "_", barcode_suffix)
TB_ALL <- RenameCells(TB_ALL, new.names = new_barcodes)

# Verify match
message("Sample Seurat barcodes: ", head(colnames(TB_ALL), 3))
message("Sample scRepertoire barcodes: ", head(combined.TCR.TB[[1]]$barcode, 3))

# Combine expression
TB_ALL <- combineExpression(
  combined.TCR.TB,
  TB_ALL,
  cloneCall  = "strict",
  chain      = "both",
  group.by   = "sample",
  cloneSize  = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500),
  proportion = FALSE
)

message("TCR data merged with Seurat. Cells with clonotype: ",
        sum(!is.na(TB_ALL$CTstrict)))

###############################################################################
#  STEP 8: SEURAT HYPEREXPANSION PLOTS
###############################################################################

setwd(hyperexp_dir)

colorblind_vector <- hcl.colors(n = 7, palette = "plasma", fixup = TRUE)

DimPlot(TB_ALL, group.by = "cloneSize", reduction = "wnn.umap") +
  scale_color_manual(values = rev(colorblind_vector[c(1, 3, 4, 5, 7)]))
ggsave("TB_Hyperexpansion.png", width = 10, height = 8, bg = "white")

DimPlot(TB_ALL, group.by = "cloneSize", reduction = "wnn.umap",
        split.by = "IGRA_infant") +
  scale_color_manual(values = rev(colorblind_vector[c(1, 3, 4, 5, 7)]))
ggsave("TB_Hyperexpansion_byIGRA_infant.png", width = 18, height = 8, bg = "white")

DimPlot(TB_ALL, group.by = "cloneSize", reduction = "wnn.umap",
        split.by = "IGRA_mother") +
  scale_color_manual(values = rev(colorblind_vector[c(1, 3, 4, 5, 7)]))
ggsave("TB_Hyperexpansion_byIGRA_mother.png", width = 18, height = 8, bg = "white")

###############################################################################
#  STEP 9: CLONAL OCCUPANCY PER CLUSTER
###############################################################################

setwd(occupy_dir)

# By WNN cluster
clonalOccupy(TB_ALL, x.axis = "seurat_clusters", label = FALSE)
ggsave("TB_Clonal_Occupancy_clusters.png", width = 17, height = 11, bg = "white")

clonalOccupy(TB_ALL, x.axis = "seurat_clusters", proportion = TRUE, label = FALSE)
ggsave("TB_Clonal_Occupancy_clusters_proportion.png", width = 17, height = 11, bg = "white")

table_cl <- clonalOccupy(TB_ALL, x.axis = "seurat_clusters", exportTable = TRUE)
write.csv(table_cl, "TB_Clones_per_Cluster.csv", row.names = FALSE)

# By Azimuth L2
clonalOccupy(TB_ALL, x.axis = "predicted.celltype.l2", label = FALSE)
ggsave("TB_Clonal_Occupancy_Azimuth.png", width = 27, height = 11, bg = "white")

clonalOccupy(TB_ALL, x.axis = "predicted.celltype.l2", proportion = TRUE, label = FALSE)
ggsave("TB_Clonal_Occupancy_Azimuth_proportion.png", width = 27, height = 11, bg = "white")

table_az <- clonalOccupy(TB_ALL, x.axis = "predicted.celltype.l2", exportTable = TRUE)
write.csv(table_az, "TB_Clones_per_Cluster_Azimuth.csv", row.names = FALSE)

###############################################################################
#  STEP 10: CLONAL OVERLAY
###############################################################################

setwd(overlay_dir)

clonalOverlay(TB_ALL, reduction = "wnn.umap", cutpoint = 1,
              bins = 10, facet.by = "Sample_ID") +
  guides(color = "none")
ggsave("TB_Clonal_Overlay_bySample.png", width = 22, height = 17, bg = "white")

clonalOverlay(TB_ALL, reduction = "wnn.umap", cutpoint = 10,
              bins = 25, facet.by = "IGRA_infant") +
  guides(color = "none")
ggsave("TB_Clonal_Overlay_byIGRA_infant.png", width = 17, height = 9, bg = "white")

clonalOverlay(TB_ALL, reduction = "wnn.umap", cutpoint = 10,
              bins = 25, facet.by = "IGRA_mother") +
  guides(color = "none")
ggsave("TB_Clonal_Overlay_byIGRA_mother.png", width = 17, height = 9, bg = "white")

###############################################################################
#  STEP 11: TREX — EPITOPE DATABASE ANNOTATION
###############################################################################

setwd(trex_dir)

TB_TRB_0 <- annotateDB(TB_ALL, chains = "TRB")
TB_TRB_1 <- annotateDB(TB_ALL, chains = "TRB", edit.distance = 1)
TB_TRB_2 <- annotateDB(TB_ALL, chains = "TRB", edit.distance = 2)

# Extract metadata for each edit distance
extract_trex <- function(obj, cluster_col = "seurat_clusters") {
  obj@meta.data %>%
    tibble::rownames_to_column("cells") %>%
    dplyr::select(
      cells,
      Sample_ID,
      IGRA_infant,
      IGRA_mother,
      cluster = all_of(cluster_col),
      CTstrict,
      clonalFrequency,
      TRB_Epitope.target,
      TRB_Epitope.sequence,
      TRB_Epitope.species,
      TRB_Tissue,
      TRB_Cell.type,
      TRB_Database
    ) %>%
    dplyr::filter(!is.na(clonalFrequency) & clonalFrequency > 1)
}

TB_TRB_df_0 <- extract_trex(TB_TRB_0)
TB_TRB_df_1 <- extract_trex(TB_TRB_1)
TB_TRB_df_2 <- extract_trex(TB_TRB_2)

# Label by edit distance
TB_TRB_df_0_lab <- TB_TRB_df_0 %>% rename_with(~ paste0(., "_ED0"), starts_with("TRB_"))
TB_TRB_df_1_lab <- TB_TRB_df_1 %>% rename_with(~ paste0(., "_ED1"), starts_with("TRB_"))
TB_TRB_df_2_lab <- TB_TRB_df_2 %>% rename_with(~ paste0(., "_ED2"), starts_with("TRB_"))

key_cols <- c("cells", "Sample_ID", "IGRA_infant", "IGRA_mother",
              "cluster", "CTstrict", "clonalFrequency")

TB_TRB_combined <- TB_TRB_df_0_lab %>%
  full_join(TB_TRB_df_1_lab, by = key_cols) %>%
  full_join(TB_TRB_df_2_lab, by = key_cols)

write_csv(TB_TRB_combined, "TB_Trex_TRB_Epitope_Database.csv")

# Pivot to long format for plots
TB_TRB_long <- TB_TRB_combined %>%
  pivot_longer(
    cols = starts_with("TRB_Epitope.target"),
    names_to = "Edit_Distance",
    names_pattern = "TRB_Epitope.target_(.*)",
    values_to = "Epitope_Target"
  ) %>%
  mutate(Epitope_Target = ifelse(is.na(Epitope_Target), "Unknown", Epitope_Target)) %>%
  distinct(CTstrict, Edit_Distance, Epitope_Target, clonalFrequency)

# All targets
ggplot(TB_TRB_long, aes(x = Epitope_Target, y = clonalFrequency, fill = CTstrict)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Edit_Distance, ncol = 1) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(y = "Clonal Frequency", x = "Epitope Target",
       title = "TB: Clonal Frequency per Epitope Target by Edit Distance")
ggsave("TB_TRB_ClonalFreq_vs_Target.png", width = 18, height = 13, dpi = 300, bg = "white")

# Excluding Unknown
TB_TRB_long %>%
  filter(Epitope_Target != "Unknown") %>%
  ggplot(aes(x = Epitope_Target, y = clonalFrequency, fill = CTstrict)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Edit_Distance, ncol = 1) +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(y = "Clonal Frequency", x = "Epitope Target",
       title = "TB: Epitope Targets (Excluding Unknown)")
ggsave("TB_TRB_ClonalFreq_vs_Target_noUnknown.png", width = 20, height = 13, dpi = 300, bg = "white")

# Summary table by Sample and IGRA
TB_TRB_long_meta <- TB_TRB_combined %>%
  pivot_longer(
    cols = starts_with("TRB_Epitope.target"),
    names_to = "Edit_Distance",
    names_pattern = "TRB_Epitope.target_(.*)",
    values_to = "Epitope_Target"
  ) %>%
  mutate(Epitope_Target = ifelse(is.na(Epitope_Target), "Unknown", Epitope_Target)) %>%
  distinct(Sample_ID, IGRA_infant, IGRA_mother, CTstrict,
           Edit_Distance, Epitope_Target, clonalFrequency)

TB_TRB_summary <- TB_TRB_long_meta %>%
  group_by(Edit_Distance, Sample_ID, IGRA_infant, IGRA_mother, Epitope_Target) %>%
  summarise(Total_Clonal_Frequency = sum(clonalFrequency, na.rm = TRUE), .groups = "drop") %>%
  arrange(Edit_Distance, Sample_ID, desc(Total_Clonal_Frequency))

write_csv(TB_TRB_summary, "TB_Trex_Epitope_Specificity_BySample.csv")

###############################################################################
#  STEP 12: EPITOPE SPECIES WAFFLE CHARTS PER SAMPLE
#  1 square = 1 unique expanded clone. All grids same size.
#  All epitope species shown (no 'Other' grouping).
###############################################################################

waffle_dir <- file.path(trex_dir, "Waffle_Charts")
dir.create(waffle_dir, recursive = TRUE, showWarnings = FALSE)
setwd(waffle_dir)

# Pivot species from combined Trex data
epitope_target_long <- TB_TRB_combined %>%
  tibble::rownames_to_column("row_id") %>%
  pivot_longer(cols = starts_with("TRB_Epitope.target"),
               names_to = "Edit_Distance", names_pattern = "TRB_Epitope.target_(.*)",
               values_to = "Epitope_Target")

epitope_species_long <- TB_TRB_combined %>%
  tibble::rownames_to_column("row_id") %>%
  pivot_longer(cols = starts_with("TRB_Epitope.species"),
               names_to = "Edit_Distance", names_pattern = "TRB_Epitope.species_(.*)",
               values_to = "Epitope_Species")

TB_TRB_long_species <- epitope_target_long %>%
  select(row_id, cells, Sample_ID, IGRA_infant, CTstrict, clonalFrequency,
         Edit_Distance, Epitope_Target) %>%
  left_join(
    epitope_species_long %>%
      select(row_id, cells, CTstrict, Edit_Distance, Epitope_Species),
    by = c("row_id", "cells", "CTstrict", "Edit_Distance")
  ) %>%
  mutate(Epitope_Target  = ifelse(is.na(Epitope_Target), "Unknown", Epitope_Target),
         Epitope_Species = ifelse(is.na(Epitope_Species), "Unknown", Epitope_Species)) %>%
  select(-row_id) %>%
  distinct()

# Save long species table
write_csv(TB_TRB_long_species, file.path(trex_dir, "TB_Trex_TRB_Species_Long.csv"))

# --- Build waffle data (ED2 only, all samples) ---
# Use ED2 for broadest matching; change to "ED0" or "ED1" if preferred
waffle_ed <- "ED2"

waffle_data <- TB_TRB_long_species %>%
  filter(Edit_Distance == waffle_ed) %>%
  filter(clonalFrequency > 1) %>%
  distinct(Sample_ID, CTstrict, Epitope_Species) %>%
  mutate(
    # Simplify multi-species labels (take first species listed)
    Epitope_Simple = sapply(strsplit(Epitope_Species, ";"), `[`, 1),
    Epitope_Simple = trimws(Epitope_Simple),
    # Collapse into TB, CMV, Other, Unknown
    Epitope_Clean = case_when(
      Epitope_Simple == "Unknown"                                                    ~ "Unknown",
      grepl("tuberculosis|M\\.tuberculosis|Mtb", Epitope_Simple, ignore.case = TRUE) ~ "TB",
      grepl("CMV|Cytomegalovirus", Epitope_Simple, ignore.case = TRUE)               ~ "CMV",
      TRUE ~ "Other"
    )
  )

# Count unique clones per sample × category
waffle_counts <- waffle_data %>%
  group_by(Sample_ID, Epitope_Clean) %>%
  summarise(n_clones = n(), .groups = "drop")

# Fixed order and palette
species_order <- c("TB", "CMV", "Other", "Unknown")
waffle_counts$Epitope_Clean <- factor(waffle_counts$Epitope_Clean, levels = species_order)

species_palette <- c(
  "TB"      = "#984EA3",
  "CMV"     = "#377EB8",
  "Other"   = "#A65628",
  "Unknown" = "#D9D9D9"
)

# --- Uniform grid dimensions ---
ncol_waffle <- 10
all_samples <- unique(waffle_counts$Sample_ID)

max_clones <- waffle_counts %>%
  group_by(Sample_ID) %>%
  summarise(total = sum(n_clones)) %>%
  pull(total) %>%
  max()

max_squares <- ceiling(max_clones / ncol_waffle) * ncol_waffle
nrow_waffle <- max_squares / ncol_waffle

cat("Waffle grid:", ncol_waffle, "cols x", nrow_waffle, "rows =", max_squares, "squares\n")

# --- Helper: make tile dataframe padded to uniform grid ---
make_waffle_df <- function(counts_vec, ncol, total_squares, level_order) {
  categories <- rep(names(counts_vec), times = counts_vec)
  n_pad      <- total_squares - length(categories)
  if (n_pad > 0) categories <- c(categories, rep("_empty", n_pad))
  n <- length(categories)
  data.frame(
    category = factor(categories, levels = c(level_order, "_empty")),
    x = ((seq_len(n) - 1) %% ncol) + 1,
    y = ((seq_len(n) - 1) %/% ncol) + 1
  )
}

palette_ext <- c(species_palette, "_empty" = NA)

# --- Generate per-sample waffle plots ---
waffle_plots <- list()

for (sid in all_samples) {
  sid_data <- waffle_counts %>%
    filter(Sample_ID == sid) %>%
    arrange(Epitope_Clean)
  
  waf_vec      <- setNames(sid_data$n_clones, as.character(sid_data$Epitope_Clean))
  total_clones <- sum(waf_vec)
  tile_df      <- make_waffle_df(waf_vec, ncol = ncol_waffle,
                                 total_squares = max_squares, level_order = species_order)
  
  p <- ggplot(tile_df, aes(x = x, y = y, fill = category)) +
    geom_tile(color = "white", linewidth = 0.8) +
    scale_fill_manual(values = palette_ext, drop = FALSE,
                      na.value = "white", guide = "none") +
    scale_y_reverse() +
    coord_equal(xlim = c(0.5, ncol_waffle + 0.5),
                ylim = c(nrow_waffle + 0.5, 0.5)) +
    labs(title = paste0("S_", sid),
         subtitle = paste0("n = ", total_clones, " expanded clones")) +
    theme_void() +
    theme(
      plot.title      = element_text(size = 24, face = "bold", hjust = 0.5),
      plot.subtitle   = element_text(size = 16, hjust = 0.5, color = "grey40",
                                     margin = margin(t = 4)),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin     = margin(8, 8, 8, 8)
    )
  
  waffle_plots[[sid]] <- p
  
  # Individual PNG
  ggsave(file.path(waffle_dir, paste0("Waffle_", sid, ".png")),
         plot = p, width = 7, height = 7, dpi = 300, bg = "white")
}

# --- Shared legend panel ---
n_items   <- length(species_order)
legend_df <- data.frame(
  Epitope = factor(species_order, levels = species_order),
  x = seq_len(n_items),
  y = 1
)

legend_panel <- ggplot(legend_df, aes(x = x, y = y, fill = Epitope)) +
  geom_point(shape = 22, size = 10, color = "grey30", stroke = 0.5) +
  geom_text(aes(label = Epitope), vjust = 2.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = species_palette) +
  scale_x_continuous(limits = c(0, n_items + 1), expand = c(0, 0)) +
  coord_cartesian(ylim = c(-0.5, 2), clip = "off") +
  theme_void() +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(5, 8, 20, 8))

# --- Combined figure: all samples + legend ---
n_plot_cols <- min(5, length(all_samples))
waffle_grid <- wrap_plots(waffle_plots, ncol = n_plot_cols)

fig_combined <- waffle_grid / legend_panel +
  plot_layout(heights = c(1, 0.18))

ggsave(file.path(waffle_dir, "Waffle_EpitopeSpecificity_AllSamples.png"),
       plot = fig_combined,
       width = n_plot_cols * 8, height = ceiling(length(all_samples) / n_plot_cols) * 8 + 3,
       dpi = 300, bg = "white")

# --- Also split by IGRA infant status ---
for (igra_status in c("IGRA+", "IGRA-")) {
  igra_samples <- sample_meta$Sample_ID[
    ifelse(sample_meta$IGRA_infant == "+", "IGRA+", "IGRA-") == igra_status
  ]
  igra_plots <- waffle_plots[as.character(igra_samples)]
  igra_plots <- igra_plots[!sapply(igra_plots, is.null)]
  
  if (length(igra_plots) > 0) {
    n_cols <- min(5, length(igra_plots))
    igra_grid <- wrap_plots(igra_plots, ncol = n_cols) / legend_panel +
      plot_layout(heights = c(1, 0.18))
    
    safe_igra <- gsub("\\+", "pos", gsub("-", "neg", igra_status))
    ggsave(file.path(waffle_dir, paste0("Waffle_Infant_", safe_igra, ".png")),
           plot = igra_grid,
           width = n_cols * 8,
           height = ceiling(length(igra_plots) / n_cols) * 8 + 3,
           dpi = 300, bg = "white")
  }
}

message("Waffle charts saved to: ", waffle_dir)

###############################################################################
#  STEP 13: SAVE
###############################################################################

qs_save(TB_ALL, file = file.path(load.path, "TB_ALL_TCR.qs2"))
message("Done! Saved: TB_ALL_TCR.qs2")

# Also save Trex-annotated objects
qs_save(TB_TRB_0, file = file.path(load.path, "TB_TRB_ED0.qs2"))
qs_save(TB_TRB_1, file = file.path(load.path, "TB_TRB_ED1.qs2"))
qs_save(TB_TRB_2, file = file.path(load.path, "TB_TRB_ED2.qs2"))
message("Trex objects saved.")
