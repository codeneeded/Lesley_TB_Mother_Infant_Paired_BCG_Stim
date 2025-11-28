library(Seurat)
library(qs2)

############################ Paths ############################

# Project working dir
setwd("/home/akshay-iyer/Documents/BCG_Stim_Infants_12_44_wks_PARSE_GEX")

# Where to save the merged Seurat object
save_dir  <- "/home/akshay-iyer/Documents/BCG_Stim_Infants_12_44_wks_PARSE_GEX/saved_R_data"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

save_path <- file.path(save_dir, "BCG_PARSE_filtered_merged.qs2")

# Base directory that contains:
#   process/
#   all-sample/
#   6071261_week12_media_treated/
#   ...
base_dir <- "/media/akshay-iyer/Elements/Parse_Data_Bone_20169/GEX_Combined"

############################ List sample folders ############################

sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

# Drop non-sample folders
sample_dirs <- sample_dirs[!basename(sample_dirs) %in% c("process", "all-sample")]

cat("Found sample folders:\n")
print(basename(sample_dirs))

############################ Build Seurat objects per sample ############################

seurat_list <- list()

for (sample_path in sample_dirs) {
  sample_name <- basename(sample_path)  # e.g. "6071261_week12_media_treated"
  
  # Split into SAMPLEID_TIMEPOINT_STIM
  parts <- strsplit(sample_name, "_")[[1]]
  sample_id <- parts[1]                                            # "6071261"
  timepoint <- ifelse(length(parts) >= 2, parts[2], NA_character_) # "week12"
  stim      <- ifelse(length(parts) > 2,
                      paste(parts[3:length(parts)], collapse = "_"),
                      NA_character_)                               # "media_treated"
  
  dge_path <- file.path(sample_path, "DGE_filtered")
  
  # Skip if required files are missing
  if (!file.exists(file.path(dge_path, "count_matrix.mtx"))) {
    warning(paste("Skipping", sample_name, "- missing DGE_filtered data at", dge_path))
    next
  }
  
  cat("Loading:", sample_name, "\n")
  
  # Read matrix
  mat <- ReadParseBio(dge_path)
  
  # Fix empty gene names
  if (any(rownames(mat) == "")) {
    rownames(mat)[rownames(mat) == ""] <- "unknown"
  }
  
  # Read cell metadata
  cell_meta <- read.csv(file.path(dge_path, "cell_metadata.csv"), row.names = 1)
  
  # Create Seurat object
  seu <- CreateSeuratObject(
    counts    = mat,
    meta.data = cell_meta,
    project   = sample_name
  )
  
  # Assign sample information
  seu$orig.ident <- sample_name
  seu$sample_id  <- sample_id
  seu$timepoint  <- timepoint
  seu$stim       <- stim
  
  # Store
  seurat_list[[sample_name]] <- seu
  
  cat("✅ Loaded", sample_name, "\n")
}

############################ Merge all Seurat objects ############################

if (length(seurat_list) == 0) {
  stop("No Seurat objects were created. Check directory structure and DGE_filtered paths.")
}

seu_merged <- Reduce(function(x, y) merge(x, y), seurat_list)
seu_merged@project.name <- "BCG_Parse_Merged_Filtered"

############################ QC Metrics ############################

# Mitochondrial content (genes starting with "MT-")
seu_merged[["percent_mito"]] <- PercentageFeatureSet(seu_merged, pattern = "^MT-")

# Ribosomal content (RPS/RPL)
seu_merged[["percent_ribo"]] <- PercentageFeatureSet(seu_merged, pattern = "^RP[SL]")

# Hemoglobin content (HB*, excluding HBP)
hb_genes <- grep("^HB", rownames(seu_merged), value = TRUE)
hb_genes <- hb_genes[!grepl("^HBP", hb_genes)]
seu_merged[["percent_hb"]] <- PercentageFeatureSet(seu_merged, features = hb_genes)

# Platelet-ish content (PECAM1 or PF4)
seu_merged[["percent_plat"]] <- PercentageFeatureSet(seu_merged, pattern = "PECAM1|PF4")

############################ Save ############################

qs_save(seu_merged, file = save_path)
#data <- qs_read("myfile.qs2")
cat("💾 Saved merged filtered Seurat object to:\n", save_path, "\n")
