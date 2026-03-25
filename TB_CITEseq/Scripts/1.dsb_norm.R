###############################################################################
#  TB CITE-seq Pipeline — Script 1: DSB Normalization
#  Adapted from OPIS CITE-seq pipeline
###############################################################################

# Load Required Libraries
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(dsb)
library(data.table)
library(ggplot2)
library(qs2)

###############################################################################
#                         PATHS & GLOBAL VARIABLES
###############################################################################

# --- INPUT: Cell Ranger output ---
in.path <- "/media/akshay-iyer/Elements/10x_TB_organized/Cell_Ranger_Out/"

# --- OUTPUT: Where to save R objects ---
out.path <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/TB_CITEseq/saved_R_data/"

# Create output directory if needed
dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

# Cell Ranger MULTI output h5 file paths
# Raw matrix: shared across all samples in the gem well
# Filtered matrix: per-sample, lives under per_sample_outs/<sample_name>/count/
r.str <- "/outs/multi/count/raw_feature_bc_matrix.h5"
# NOTE: Filtered path is built per-sample inside the loop as:
# /outs/per_sample_outs/<SAMPLE>/count/sample_filtered_feature_bc_matrix.h5

# Get sample folder names
get_folder_names <- function(in.path) {
  folder_names <- list.dirs(in.path, full.names = FALSE, recursive = FALSE)
  folder_names <- folder_names[folder_names != ""]
  return(folder_names)
}

f_names <- get_folder_names(in.path)
print(f_names)

###############################################################################
#  STEP 0: DIAGNOSTIC — Check ADT feature names before committing to mapping
#  Run this section FIRST, inspect output, then proceed.
###############################################################################

# Read one sample to inspect ADT feature names
test_sample <- f_names[1]
message("Reading test sample: ", test_sample)

test_raw <- Read10X_h5(paste0(in.path, test_sample, r.str))

# Print all Antibody Capture feature names
cat("\n========== ADT Feature Names in", test_sample, "==========\n")
adt_names <- rownames(test_raw$`Antibody Capture`)
cat(paste(adt_names, collapse = "\n"), "\n")
cat("\nTotal ADT features:", length(adt_names), "\n")

# --- Compare against the OPIS panel_to_gene mapping ---
# This is the same mapping from OPIS; we check overlap
panel_to_gene <- c(
  "Hu.CD101"           = "CD101",
  "Hu.CD103"           = "ITGAE",
  "Hu.CD105_43A3"      = "ENG",
  "Hu.CD107a"          = "LAMP1",
  "Hu.CD112"           = "NECTIN2",
  "Hu.CD119"           = "IFNGR1",
  "Hu.CD11a"           = "ITGAL",
  "Hu.CD11b"           = "ITGAM",
  "Hu.CD11c"           = "ITGAX",
  "Hu.CD122"           = "IL2RB",
  "Hu.CD123"           = "IL3RA",
  "Hu.CD124"           = "IL4R",
  "Hu.CD127"           = "IL7R",
  "Hu.CD13"            = "ANPEP",
  "Hu.CD134"           = "TNFRSF4",
  "Hu.CD137"           = "TNFRSF9",
  "Hu.CD141"           = "THBD",
  "Hu.CD146"           = "MCAM",
  "Hu.CD14_M5E2"       = "CD14",
  "Hu.CD152"           = "CTLA4",
  "Hu.CD154"           = "CD40LG",
  "Hu.CD155"           = "PVR",
  "Hu.CD158"           = "KIR2DL1",
  "Hu.CD158b"          = "KIR2DL3",
  "Hu.CD158e1"         = "KIR3DL1",
  "Hu.CD16"            = "FCGR3A",
  "Hu.CD161"           = "KLRB1",
  "Hu.CD163"           = "CD163",
  "Hu.CD169"           = "SIGLEC1",
  "Hu.CD18"            = "ITGB2",
  "Hu.CD183"           = "CXCR3",
  "Hu.CD185"           = "CXCR5",
  "Hu.CD19"            = "CD19",
  "Hu.CD194"           = "CCR4",
  "Hu.CD195"           = "CCR5",
  "Hu.CD196"           = "CCR6",
  "Hu.CD1c"            = "CD1C",
  "Hu.CD1d"            = "CD1D",
  "Hu.CD2"             = "CD2",
  "Hu.CD20_2H7"        = "MS4A1",
  "Hu.CD21"            = "CR2",
  "Hu.CD22"            = "CD22",
  "Hu.CD223"           = "LAG3",
  "Hu.CD224"           = "GGT1",
  "Hu.CD226_11A8"      = "CD226",
  "Hu.CD23"            = "FCER2",
  "Hu.CD24"            = "CD24",
  "Hu.CD244"           = "CD244",
  "Hu.CD25"            = "IL2RA",
  "Hu.CD26"            = "DPP4",
  "Hu.CD267"           = "TNFRSF13B",
  "Hu.CD268"           = "TNFRSF13C",
  "Hu.CD27"            = "CD27",
  "Hu.CD270"           = "TNFRSF14",
  "Hu.CD272"           = "BTLA",
  "Hu.CD274"           = "CD274",
  "Hu.CD279"           = "PDCD1",
  "Hu.CD28"            = "CD28",
  "Hu.CD29"            = "ITGB1",
  "Hu.CD303"           = "CLEC4C",
  "Hu.CD31"            = "PECAM1",
  "Hu.CD314"           = "KLRK1",
  "Hu.CD319"           = "SLAMF7",
  "Hu.CD32"            = "FCGR2A",
  "Hu.CD328"           = "SIGLEC7",
  "Hu.CD33"            = "CD33",
  "Hu.CD335"           = "NCR1",
  "Hu.CD35"            = "CR1",
  "Hu.CD352"           = "SLAMF6",
  "Hu.CD36"            = "CD36",
  "Hu.CD38_HIT2"       = "CD38",
  "Hu.CD39"            = "ENTPD1",
  "Hu.CD3_UCHT1"       = "CD3D",
  "Hu.CD40"            = "CD40",
  "Hu.CD41"            = "ITGA2B",
  "Hu.CD42b"           = "GP1BB",
  "Hu.CD45RA"          = "CD45RA",
  "Hu.CD45RO"          = "CD45RO",
  "Hu.CD45_HI30"       = "CD45",
  "Hu.CD47"            = "CD47",
  "Hu.CD48"            = "CD48",
  "Hu.CD49a"           = "ITGA1",
  "Hu.CD49b"           = "ITGA2",
  "Hu.CD49d"           = "ITGA4",
  "Hu.CD4_RPA.T4"      = "CD4",
  "Hu.CD5"             = "CD5",
  "Hu.CD52"            = "CD52",
  "Hu.CD54"            = "ICAM1",
  "Hu.CD56"            = "NCAM1",
  "Hu.CD57"            = "B3GAT1",
  "Hu.CD58"            = "CD58",
  "Hu.CD62L"           = "SELL",
  "Hu.CD62P"           = "SELP",
  "Hu.CD64"            = "FCGR1A",
  "Hu.CD69"            = "CD69",
  "Hu.CD7"             = "CD7",
  "Hu.CD71"            = "TFRC",
  "Hu.CD73"            = "NT5E",
  "Hu.CD79b"           = "CD79B",
  "Hu.CD8"             = "CD8A",
  "Hu.CD81"            = "CD81",
  "Hu.CD82"            = "CD82",
  "Hu.CD83"            = "CD83",
  "Hu.CD85j"           = "LILRB1",
  "Hu.CD86"            = "CD86",
  "Hu.CD88"            = "C5AR1",
  "Hu.CD94"            = "KLRD1",
  "Hu.CD95"            = "FAS",
  "Hu.CD99"            = "CD99",
  "Hu.CLEC12A"         = "CLEC12A",
  "Hu.CX3CR1"          = "CX3CR1",
  "Hu.FceRIa"          = "FCER1A",
  "Hu.GPR56"           = "ADGRG1",
  "Hu.HLA.ABC"         = "HLA-A",
  "Hu.HLA.DR"          = "HLA-DRA",
  "Hu.HLA.E"           = "HLA-E",
  "Hu.Ig.LightChain.k" = "IGKC",
  "Hu.Ig.LightChain.l" = "IGLC1",
  "Hu.IgD"             = "IGHD",
  "Hu.IgM"             = "IGHM",
  "Hu.KLRG1"           = "KLRG1",
  "Hu.LOX.1"           = "OLR1",
  "Hu.TCR.AB"          = "TCR-AB",
  "Hu.TCR.Va7.2"       = "TCR-vA7.2",
  "Hu.TCR.Vd2"         = "TCR-vD2",
  "Hu.TIGIT"           = "TIGIT",
  "HuMs.CD44"          = "CD44",
  "HuMs.CD49f"         = "ITGA6",
  "HuMs.integrin.b7"   = "ITGB7",
  "HuMsRt.CD278"       = "ICOS",
  "Isotype_HTK888"     = "Hamster-IgG",
  "Isotype_MOPC.173"   = "Mouse-IgG2a",
  "Isotype_MOPC.21"    = "Mouse-IgG1",
  "Isotype_MPC.11"     = "Mouse-IgG2b",
  "Isotype_RTK2071"    = "Rat-IgG1",
  "Isotype_RTK2758"    = "Rat-IgG2a",
  "Isotype_RTK4530"    = "Rat-IgG2b"
)

# Check which names match, which don't
matched   <- adt_names[adt_names %in% names(panel_to_gene)]
unmatched <- adt_names[!adt_names %in% names(panel_to_gene)]

cat("\n========== MAPPING DIAGNOSTIC ==========\n")
cat("Matched:  ", length(matched), "/", length(adt_names), "\n")
if (length(unmatched) > 0) {
  cat("\n*** UNMATCHED ADT names (need manual mapping or name fix): ***\n")
  cat(paste(unmatched, collapse = "\n"), "\n")
} else {
  cat("All ADT names match the panel_to_gene mapping.\n")
}

# Also check: are there panel_to_gene keys NOT present in the data?
missing_from_data <- setdiff(names(panel_to_gene), adt_names)
if (length(missing_from_data) > 0) {
  cat("\n*** Panel entries NOT found in TB data (may be absent from this panel): ***\n")
  cat(paste(missing_from_data, collapse = "\n"), "\n")
}

cat("\n=== STOP HERE. Review the output above. ===\n")
cat("If all names match, proceed with the rest of the script.\n")
cat("If there are unmatched names, update panel_to_gene before continuing.\n\n")

###############################################################################
#  STEP 1: Load Data + Basic Pre-processing
#  (Only run after confirming panel_to_gene mapping is correct)
###############################################################################

for (i in f_names) {
  message("Loading: ", i)
  
  # Read data from 10x MULTI outputs
  # Raw: shared across gem well
  raw   <- Read10X_h5(paste0(in.path, i, r.str))
  # Filtered: per-sample path
  f.str.i <- paste0("/outs/per_sample_outs/", i, "/count/sample_filtered_feature_bc_matrix.h5")
  cells <- Read10X_h5(paste0(in.path, i, f.str.i))
  
  # Replace antibody capture names with standardized names
  raw$`Antibody Capture`@Dimnames[[1]]   <- unname(panel_to_gene[raw$`Antibody Capture`@Dimnames[[1]]])
  cells$`Antibody Capture`@Dimnames[[1]] <- unname(panel_to_gene[cells$`Antibody Capture`@Dimnames[[1]]])
  
  # Define cell-containing barcodes and background
  stained_cells <- colnames(cells$`Gene Expression`)
  background    <- setdiff(colnames(raw$`Gene Expression`), stained_cells)
  
  # Assign final processed data
  prot <- raw$`Antibody Capture`
  rna  <- raw$`Gene Expression`
  
  # Create metadata
  rna.size     <- log10(Matrix::colSums(rna))
  prot.size    <- log10(Matrix::colSums(prot))
  nCount_RNA   <- Matrix::colSums(rna)
  nCount_ADT   <- Matrix::colSums(prot)
  nFeature_RNA <- Matrix::colSums(rna > 0)
  nFeature_ADT <- Matrix::colSums(prot > 0)
  mtgene       <- grep(pattern = "^MT-", rownames(rna), value = TRUE)
  mt.prop      <- Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
  
  md <- as.data.frame(cbind(nCount_RNA, nFeature_RNA, nCount_ADT, nFeature_ADT,
                            rna.size, prot.size, mt.prop))
  md$drop.class <- ifelse(rownames(md) %in% stained_cells, 'cell', 'background')
  md <- md[md$rna.size > 0 & md$prot.size > 0, ]
  
  assign(paste0(i, ".md"), md)
  assign(paste0(i, ".prot"), prot)
  assign(paste0(i, ".rna"), rna)
}

###############################################################################
#  STEP 2: Droplet QC Plots
###############################################################################

# QC directories
base_path           <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/TB_CITEseq"
qc_dir              <- file.path(base_path, "QC")
droplet_settings_dir <- file.path(qc_dir, "Droplet_Settings")
dir.create(droplet_settings_dir, recursive = TRUE, showWarnings = FALSE)

for (i in f_names) {
  md   <- get(paste0(i, ".md"))
  prot <- get(paste0(i, ".prot"))
  rna  <- get(paste0(i, ".rna"))
  
  png(file.path(droplet_settings_dir, paste0(i, "_genevsprotlibsize.png")),
      width = 800, height = 600)
  p <- ggplot(md, aes(x = log10(nFeature_RNA), y = prot.size)) +
    theme_bw() +
    geom_bin2d(bins = 300) +
    scale_fill_viridis_c(option = "C") +
    facet_wrap(~drop.class) +
    ggtitle(paste0(i, " — Gene vs Protein Library Size"))
  print(p)
  dev.off()
}

###############################################################################
#  STEP 3: Define Background Droplets for DSB
###############################################################################

# ** IMPORTANT: Review the droplet plots from Step 2 before setting thresholds.
# The thresholds below are carried from OPIS — adjust if your TB data
# shows different background/cell distributions.

process_entry <- function(entry_name) {
  md <- get(paste0(entry_name, ".md"))
  
  # Background thresholds — REVIEW AND ADJUST if needed
  filtered_md <- md[md$prot.size > 1.5 & md$prot.size < 4 & md$rna.size < 2.5, ]
  background_drops <- rownames(filtered_md)
  
  assign(paste0(entry_name, ".background_drops"), background_drops, envir = .GlobalEnv)
}

for (entry_name in f_names) {
  process_entry(entry_name)
}

###############################################################################
#  STEP 4: QC Filtering + DSB Normalization
###############################################################################

setwd(droplet_settings_dir)

for (i in f_names) {
  message("Processing QC for: ", i)
  
  md               <- get(paste0(i, ".md"))
  prot             <- get(paste0(i, ".prot"))
  rna              <- get(paste0(i, ".rna"))
  background_drops <- get(paste0(i, ".background_drops"))
  
  background.adt.mtx <- as.matrix(prot[, background_drops])
  cellmd             <- md[md$drop.class == 'cell', ]
  
  # MAD-based filtering (±3 MAD from median)
  rna.mult  <- 3 * mad(cellmd$rna.size)
  prot.mult <- 3 * mad(cellmd$prot.size)
  rna.lower  <- median(cellmd$rna.size)  - rna.mult
  rna.upper  <- median(cellmd$rna.size)  + rna.mult
  prot.lower <- median(cellmd$prot.size) - prot.mult
  prot.upper <- median(cellmd$prot.size) + prot.mult
  
  qc_cells <- rownames(
    cellmd[cellmd$prot.size > prot.lower &
             cellmd$prot.size < prot.upper &
             cellmd$rna.size  > rna.lower  &
             cellmd$rna.size  < rna.upper  &
             cellmd$mt.prop   < 0.25, ]
  )
  
  # QC threshold plots
  png(file.path(droplet_settings_dir, paste0(i, "_qc_thresholds.png")),
      width = 800, height = 600)
  plot_aes <- list(theme_bw(),
                   geom_point(shape = 21, stroke = 0, size = 0.7),
                   scale_fill_viridis_c(option = "C"))
  p1 <- ggplot(cellmd, aes(x = rna.size))  + geom_histogram(bins = 50) + theme_bw() + xlab("log10 RNA library size")
  p2 <- ggplot(cellmd, aes(x = mt.prop))   + geom_histogram(bins = 50) + theme_bw() + xlab("mitochondrial read proportion")
  p3 <- ggplot(cellmd, aes(x = log10(nFeature_RNA), y = rna.size, fill = mt.prop))  + plot_aes
  p4 <- ggplot(cellmd, aes(x = nFeature_RNA, y = prot.size, fill = mt.prop))        + plot_aes
  print(p1 + p2 + p3 + p4)
  dev.off()
  
  cell.adt.raw <- as.matrix(prot[, qc_cells])
  cell.rna.raw <- rna[, qc_cells]
  cellmd       <- cellmd[qc_cells, ]
  
  # Check protein staining
  pm  <- sort(apply(cell.adt.raw, 1, max))
  pm2 <- apply(background.adt.mtx, 1, max)
  
  assign(paste0(i, ".cell.adt.raw"),       cell.adt.raw)
  assign(paste0(i, ".cell.rna.raw"),       cell.rna.raw)
  assign(paste0(i, ".background.adt.mtx"), background.adt.mtx)
  assign(paste0(i, ".cellmd"),             cellmd)
  assign(paste0(i, ".pm"),                 pm)
}

###############################################################################
#  STEP 5: DSB Normalization
###############################################################################

# Isotype controls (same panel)
isotype.controls <- c('Mouse-IgG1', 'Mouse-IgG2a', 'Mouse-IgG2b',
                      'Rat-IgG2b', 'Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')

drop_log <- list()

for (i in f_names) {
  message("Running DSB for: ", i)
  
  cell.adt.raw       <- get(paste0(i, ".cell.adt.raw"))
  background.adt.mtx <- get(paste0(i, ".background.adt.mtx"))
  
  # Compute max and 99th percentile
  max_cell <- apply(cell.adt.raw, 1, max)
  max_bg   <- apply(background.adt.mtx, 1, max)
  q99_cell <- apply(cell.adt.raw, 1, function(x) quantile(x, 0.99, na.rm = TRUE))
  q99_bg   <- apply(background.adt.mtx, 1, function(x) quantile(x, 0.99, na.rm = TRUE))
  
  # Identify proteins to drop
  zero_drop    <- (max_cell == 0 & max_bg == 0)
  lowq99_drop  <- (q99_cell < 1  & q99_bg < 1)
  non_staining <- zero_drop | lowq99_drop
  
  # Never drop isotypes
  non_staining[rownames(cell.adt.raw) %in% isotype.controls] <- FALSE
  
  dropped <- rownames(cell.adt.raw)[non_staining]
  
  if (length(dropped) > 0) {
    message("  Dropping proteins in ", i, ": ", paste(dropped, collapse = ", "))
    
    drop_log[[i]] <- data.frame(
      sample   = i,
      protein  = dropped,
      reason   = ifelse(zero_drop[dropped], "zero", "low_q99"),
      max_cell = max_cell[dropped],
      max_bg   = max_bg[dropped],
      q99_cell = q99_cell[dropped],
      q99_bg   = q99_bg[dropped],
      stringsAsFactors = FALSE
    )
    
    cell.adt.raw       <- cell.adt.raw[!non_staining, , drop = FALSE]
    background.adt.mtx <- background.adt.mtx[!non_staining, , drop = FALSE]
  }
  
  # Run DSB
  cells.dsb.norm <- DSBNormalizeProtein(
    cell_protein_matrix      = cell.adt.raw,
    empty_drop_matrix        = background.adt.mtx,
    denoise.counts           = TRUE,
    use.isotype.control      = TRUE,
    isotype.control.name.vec = isotype.controls
  )
  
  cells.dsb.norm <- Matrix(as.matrix(cells.dsb.norm), sparse = TRUE)
  assign(paste0(i, ".cells.dsb.norm"), cells.dsb.norm)
}

# Save drop log
if (length(drop_log) > 0) {
  drop_log_df <- do.call(rbind, drop_log)
  write.csv(drop_log_df,
            file = file.path(droplet_settings_dir, "DSB_dropped_proteins_log.csv"),
            row.names = FALSE)
  message("Saved: DSB_dropped_proteins_log.csv")
} else {
  message("No proteins were dropped across samples.")
}

###############################################################################
#  STEP 6: Create + Merge Seurat Objects
###############################################################################

for (i in f_names) {
  cellmd         <- get(paste0(i, ".cellmd"))
  cell.adt.raw   <- get(paste0(i, ".cell.adt.raw"))
  cells.dsb.norm <- get(paste0(i, ".cells.dsb.norm"))
  cell.rna.raw   <- get(paste0(i, ".cell.rna.raw"))
  
  stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.adt.raw))))
  stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.rna.raw))))
  
  s <- CreateSeuratObject(counts   = cell.rna.raw,
                          meta.data = cellmd,
                          assay     = "RNA",
                          min.cells = 20)
  
  s[["ADT"]]    <- CreateAssayObject(data = cells.dsb.norm)
  s$orig.ident  <- i
  
  assign(paste0(i, ".seurat"), s)
}

# Merge all sample Seurat objects
merged_seurat <- merge(
  x = get(paste0(f_names[1], ".seurat")),
  y = lapply(f_names[-1], function(name) get(paste0(name, ".seurat"))),
  add.cell.id = f_names
)

###############################################################################
#  STEP 7: Add Basic Metadata
###############################################################################

# Genes per UMI
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

# Percent mitochondrial
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^MT-", col.name = "percent_mito")

# Percent ribosomal
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^RP[SL]", col.name = "percent_ribo")

# Percent hemoglobin
merged_seurat <- PercentageFeatureSet(merged_seurat, "^HB[^(P)]", col.name = "percent_hb")

# Percent platelet
merged_seurat <- PercentageFeatureSet(merged_seurat, "PECAM1|PF4", col.name = "percent_plat")

levels(as.factor(merged_seurat$orig.ident))

###############################################################################
#  STEP 7b: IGRA Status + Sample ID Metadata
###############################################################################

# Sample metadata lookup table
sample_meta <- data.frame(
  orig.ident    = c("Bone-19321-001", "Bone-19321-002", "Bone-19321-003", "Bone-19321-004",
                    "Bone-19430-001", "Bone-19430-002", "Bone-19429-001", "Bone-19429-002",
                    "Bone-19649-001", "Bone-19708-001"),
  Sample_ID     = c("6062151", "6067021", "6063421", "6063411",
                    "6061161", "6063481", "6044931", "6044981",
                    "6071261", "6065761"),
  IGRA_infant   = c("+", "+", "-", "-", "+", "-", "-", "-", "+", "+"),
  IGRA_mother   = c("+", "-", "+", "-", "+", "+", "-", "+", "+", "+"),
  stringsAsFactors = FALSE
)

# Map metadata onto Seurat object
cell_orig <- merged_seurat$orig.ident
idx       <- match(cell_orig, sample_meta$orig.ident)

# Sample ID
merged_seurat$Sample_ID <- sample_meta$Sample_ID[idx]

# Infant IGRA status
merged_seurat$IGRA_infant <- ifelse(sample_meta$IGRA_infant[idx] == "+", "IGRA+", "IGRA-")
merged_seurat$IGRA_infant <- factor(merged_seurat$IGRA_infant, levels = c("IGRA-", "IGRA+"))

# Maternal IGRA status
merged_seurat$IGRA_mother <- ifelse(sample_meta$IGRA_mother[idx] == "+", "IGRA+", "IGRA-")
merged_seurat$IGRA_mother <- factor(merged_seurat$IGRA_mother, levels = c("IGRA-", "IGRA+"))

# Combined IGRA group (useful for stratified analyses)
merged_seurat$IGRA_group <- paste0("Infant:", merged_seurat$IGRA_infant,
                                   "_Mother:", merged_seurat$IGRA_mother)
merged_seurat$IGRA_group <- factor(merged_seurat$IGRA_group)

# Quick check
message("Metadata added. Cell counts per IGRA group:")
print(table(merged_seurat$IGRA_infant, merged_seurat$IGRA_mother,
            dnn = c("Infant_IGRA", "Mother_IGRA")))

###############################################################################
#  STEP 8: Save Merged Object
###############################################################################

qs_save(merged_seurat,
        file = file.path(out.path, "TB_Seuratv5_CITEseq_dsbnorm_merged.qs2"))

message("Done! Merged Seurat object saved.")
