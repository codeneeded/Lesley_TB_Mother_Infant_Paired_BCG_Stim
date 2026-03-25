###############################################################################
#  TB CITE-seq Pipeline — Script 4: Isotype Correction
#  Adapted from OPIS CITE-seq pipeline
###############################################################################

library(Seurat)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(readxl)
library(qs2)

###############################################################################
#                         PATHS & GLOBAL VARIABLES
###############################################################################

base_dir  <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/TB_CITEseq"
load.path <- file.path(base_dir, "saved_R_data")
iso_dir   <- file.path(base_dir, "Integration", "Isotype")

dir.create(iso_dir, recursive = TRUE, showWarnings = FALSE)
setwd(iso_dir)

###############################################################################
#                         LOAD SEURAT OBJECT
###############################################################################

seurat_layer <- qs_read(file.path(load.path, "TB_Seuratv5_filtered_CellCycle_DoubletClean.qs2"))

# Drop all assays except RNA and ADT (remove Azimuth score assays etc.)
assays_to_keep <- c("RNA", "ADT")
current_assays <- names(seurat_layer@assays)
assays_to_drop <- setdiff(current_assays, assays_to_keep)

for (assay in assays_to_drop) {
  seurat_layer[[assay]] <- NULL
}

message("Assays retained: ", paste(names(seurat_layer@assays), collapse = ", "))

###############################################################################
#                  ISOTYPE CONTROL THRESHOLD VISUALIZATION
###############################################################################

isotype_genes <- c('Mouse-IgG1', 'Mouse-IgG2a', 'Mouse-IgG2b',
                   'Rat-IgG2b', 'Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')

# Extract isotype data and reshape for plotting
isotype_data <- seurat_layer[["ADT"]]$data[isotype_genes, ]

plot_data <- as.data.frame(t(isotype_data)) %>%
  rownames_to_column("CellBarcode") %>%
  pivot_longer(cols = -CellBarcode, names_to = "Isotype", values_to = "Expression")

# Add Sample_ID for faceting
plot_data$Sample_ID <- seurat_layer@meta.data[plot_data$CellBarcode, "Sample_ID", drop = TRUE]

# Calculate 99% threshold per isotype per sample
threshold_data_list <- lapply(isotype_genes, function(isotype) {
  thresholds <- tapply(
    seurat_layer[["ADT"]]$data[isotype, ],
    seurat_layer@meta.data$Sample_ID,
    function(x) quantile(x, 0.99)
  )
  data.frame(Isotype = isotype, Sample_ID = names(thresholds),
             Threshold = as.numeric(thresholds))
})
threshold_data <- do.call(rbind, threshold_data_list)

# Plot
p_iso <- ggplot(plot_data, aes(x = Isotype, y = Expression)) +
  geom_violin(scale = "width", fill = "#92C5DE", color = "grey40", linewidth = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.3, color = "grey50") +
  geom_line(data = threshold_data,
            aes(y = Threshold, group = Sample_ID),
            color = "#D6604D", linetype = "solid", linewidth = 0.6) +
  facet_wrap(~ Sample_ID) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y   = element_text(size = 11),
    strip.text     = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Isotype Control Expression with 99% Thresholds",
       y = "DSB Normalized Expression", x = NULL)

ggsave(file.path(iso_dir, "Isotype_Threshold.png"),
       plot = p_iso, dpi = 300, width = 17, height = 10, bg = "white")

###############################################################################
#                         ISOTYPE CORRECTION
###############################################################################

# Load protein-to-isotype mapping
# ** UPDATE this path to your TB isotype mapping file **
# The Excel file should have two columns: 'names' (protein) and 'isotype'
isotype_df <- read_excel(file.path(base_dir, "Isotype.xlsx"))

# Named vector: protein -> isotype
protein_to_isotype <- setNames(isotype_df$isotype, isotype_df$names)

# Check which proteins are in the ADT assay
adt_proteins <- rownames(seurat_layer[["ADT"]]@data)

missing_proteins <- setdiff(names(protein_to_isotype), adt_proteins)
if (length(missing_proteins) > 0) {
  message("Proteins in Isotype.xlsx but not in ADT (skipped):\n",
          paste(missing_proteins, collapse = ", "))
}
protein_to_isotype <- protein_to_isotype[names(protein_to_isotype) %in% adt_proteins]

# Calculate 99% threshold per isotype per sample
names(isotype_genes) <- isotype_genes
thresholds_list <- lapply(isotype_genes, function(isotype) {
  tapply(
    seurat_layer[["ADT"]]@data[isotype, ],
    seurat_layer@meta.data$Sample_ID,
    function(x) quantile(x, 0.99)
  )
})
names(thresholds_list) <- isotype_genes

# Ensure mapped isotypes exist in thresholds
protein_to_isotype <- protein_to_isotype[protein_to_isotype %in% names(thresholds_list)]

# Copy object and apply correction
seurat_isotype <- seurat_layer

for (sample in unique(seurat_layer@meta.data$Sample_ID)) {
  sample_cells <- which(seurat_layer@meta.data$Sample_ID == sample)
  
  for (protein in names(protein_to_isotype)) {
    iso <- protein_to_isotype[[protein]]
    
    if (!protein %in% rownames(seurat_isotype[["ADT"]]@data)) next
    if (is.null(thresholds_list[[iso]][sample])) next
    
    thr <- thresholds_list[[iso]][sample]
    
    # Subtract threshold
    seurat_isotype[["ADT"]]@data[protein, sample_cells] <-
      seurat_isotype[["ADT"]]@data[protein, sample_cells] - thr
    
    # Zero out negatives
    neg_idx <- which(seurat_isotype[["ADT"]]@data[protein, ] < 0)
    if (length(neg_idx) > 0) {
      seurat_isotype[["ADT"]]@data[protein, neg_idx] <- 0
    }
  }
}

###############################################################################
#                              SAVE
###############################################################################

qs_save(seurat_isotype, file = file.path(load.path, "TB_Seurat_isotype.qs2"))
message("Done! Saved: TB_Seurat_isotype.qs2")
