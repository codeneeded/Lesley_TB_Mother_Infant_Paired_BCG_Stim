## ---------------------------
## Project: Lesley TB — Maternal→Infant (12w BCG CD8 proliferation)
## Script: 02_build_master_Q1_prolif12w_MaternalJoin_SUFFIX.R
## ---------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(readr); library(purrr); library(digest)
})
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(glmnet)
library(pROC)
library(readr)
library(purrr)

## ---- Paths to saved inputs ----
saved_dir <- "C:/Users/ammas/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/saved_R_dat"

# Cytokine (preprocessed LONG with MED-normalization & parsed gates)
cyto_rds <- file.path(saved_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.rds")
cyto_csv <- file.path(saved_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.csv")
cytokine_long <- if (file.exists(cyto_rds)) readRDS(cyto_rds) else read_csv(cyto_csv, show_col_types = FALSE)

# Plasma (wide, z-scored within TP)
plasma_rds <- file.path(saved_dir, "Plasma_Biomarker_wide_scaled.rds")
plasma_csv <- file.path(saved_dir, "Plasma_Biomarker_wide_scaled.csv")
plasma_wide <- if (file.exists(plasma_rds)) readRDS(plasma_rds) else read_csv(plasma_csv, show_col_types = FALSE)

## ---- Helpers ----
norm_tp <- function(x) {
  z <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    z %in% c("entry") ~ "Entry",
    z %in% c("dx","diagnosis") ~ "Dx",
    z %in% c("12","12 wks","12 weeks","12wks") ~ "12 Wks",
    z %in% c("44","44 wks","44 weeks","44wks") ~ "44 Wks",
    TRUE ~ as.character(x)
  )
}

normalize_gatepath <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[\\u00A0\\p{Zs}]+", " ") %>%
    stringr::str_replace_all("\\s*,\\s*", ", ") %>%
    stringr::str_replace_all("\\s*\\|\\s*", " | ") %>%
    stringr::str_squish()
}

# Extract suffix AFTER the "CD45, Live Dead subset/" node; if not present, use entire path before metric.
suffix_after_ld <- function(gp_norm) {
  gp_norm <- normalize_gatepath(gp_norm)
  # remove trailing metric " | Freq. of ... (%)"
  core <- str_replace(gp_norm, "\\s*\\|\\s*Freq\\..*$", "")
  # split at the LD node; keep the right side if present
  parts <- str_split(core, "CD45,\\s*Live\\s*Dead\\s*subset/", n = 2, simplify = TRUE)
  suffix <- if (ncol(parts) == 2 && nzchar(parts[2])) parts[2] else core
  str_squish(suffix)
}

# Compact, unique feature name derived from suffix (plus short hash guard)
suffix_to_feat <- function(sfx) {
  base <- sfx %>%
    str_replace_all("[/\\s]+", "_") %>%
    str_replace_all("[^A-Za-z0-9_]+", "") %>%
    str_replace("^_+|_+$", "")
  paste0(base, "_", substr(digest(sfx, algo = "md5"), 1, 6))
}

## ---- Normalize common columns and derive suffix ----
cytokine_long <- cytokine_long %>%
  mutate(
    PID        = as.character(PID),
    PID_Type   = ifelse(!is.na(PID_Type), PID_Type, `PID Type`),
    Time_Point = norm_tp(ifelse(!is.na(Time_Point), Time_Point, `Time Point`)),
    GatePath_norm = normalize_gatepath(GatePath),
    GateSuffix    = suffix_after_ld(GatePath_norm)
  )
# ---- Harmonize column names (create helper fields without spaces) ----
cytokine_long <- cytokine_long %>%
  dplyr::mutate(
    PID_Type   = dplyr::coalesce(PID_Type, `PID Type`),
    Time_Point = dplyr::coalesce(Time_Point, `Time Point`)
  )
str(cytokine_long)
## --- Recompute a robust suffix (after CD45, Live Dead subset) ---
suffix_after_ld_safe <- function(gp_norm) {
  core <- sub("\\s*\\|\\s*Freq\\..*$", "", as.character(gp_norm), perl = TRUE)  # drop metric tail
  # try after LD node
  suf  <- sub("^.*CD45,\\s*Live\\s*Dead\\s*subset/\\s*", "", core, perl = TRUE)
  # fallback if LD node not present
  if (identical(suf, core)) {
    suf <- sub("^.*Single\\s*Cells/\\s*", "", core, perl = TRUE)
  }
  trimws(suf)
}

# make sure helper columns exist and recompute suffix
cytokine_long <- cytokine_long %>%
  dplyr::mutate(
    PID_Type     = dplyr::coalesce(PID_Type, `PID Type`),
    Time_Point   = dplyr::coalesce(Time_Point, `Time Point`),
    GatePath_norm= normalize_gatepath(GatePath),
    GateSuffix   = suffix_after_ld_safe(GatePath_norm)
  )
## =========================
## A) Infant outcome (12 Wks, BCG) — suffix-only match
## =========================
## --- Outcome: CD8 Proliferating with tolerant match (gd TCR or gd TCR-) ---
rx_suffix <- "gd TCR-/CD8/Proliferating$"

infant_outcome <- cytokine_long %>%
  dplyr::filter(
    PID_Type   == "Infant",
    Time_Point %in% c("12 Wks","44 Wks"),     # keep both for paired analyses
    Condition  == "BCG",
    stringr::str_detect(GateSuffix, rx_suffix),
    stringr::str_detect(GatePath_norm, "\\|\\s*Freq\\.\\s*of\\s*Parent")
  ) %>%
  dplyr::arrange(PID, Time_Point, dplyr::desc(GateDepth)) %>%
  dplyr::distinct(PID, Time_Point, .keep_all = TRUE) %>%
  dplyr::transmute(
    Infant_PID    = PID,
    Time_Point    = Time_Point,
    Maternal_PID  = as.character(`Maternal PID`),
    Batch_infant  = as.character(Batch),
    Infant_IGRA   = `Infant IGRA`,
    Maternal_IGRA = `Maternal IGRA`,
    y_abs   = as.numeric(Value_RawNormalized),
    y_delta = as.numeric(Delta_for_model)
  )

#🔹 y_abs — Absolute frequency of proliferating CD8⁺ T cells (after batch normalization)
#🔹 y_delta — Δ Stimulated–Baseline response (BCG minus MED)

# --- Compute noise threshold per batch using controls from the same gate ---
ctrl_noise <- cytokine_long %>%
  dplyr::filter(
    PID_Type == "Batch Control",
    Condition == "MED",
    stringr::str_detect(GateSuffix, "(?i)^gd\\s*TCR-?/CD8/Proliferating$"),
    stringr::str_detect(GatePath_norm, "\\|\\s*Freq\\.\\s*of\\s*Parent")
  ) %>%
  dplyr::group_by(Batch) %>%
  dplyr::summarise(
    noise95 = suppressWarnings(stats::quantile(Value, 0.95, na.rm = TRUE)),  # 95th percentile
    .groups = "drop"
  ) %>%
  dplyr::mutate(Batch = as.character(Batch))

# --- Join to infants and define y_bin ---
infant_outcome <- infant_outcome %>%
  dplyr::left_join(ctrl_noise, by = c("Batch_infant" = "Batch")) %>%
  dplyr::mutate(
    noise95 = dplyr::coalesce(noise95, 1.0),      # fallback if no control in batch
    y_bin = dplyr::case_when(
      # true responders: above batch noise + small margin (0.5%) and positive delta
      is.finite(y_abs) & is.finite(y_delta) &
        (y_abs >= (noise95 + 0.5) & y_delta >= 0.5) ~ 1L,
      # conservative fallback if absolute ≥ 2%
      is.finite(y_abs) & y_abs >= 2 ~ 1L,
      TRUE ~ 0L
    )
  )
# --- 1) Remove infants with no linked maternal PID ---
infant_outcome <- infant_outcome %>%
  dplyr::filter(!is.na(Maternal_PID) & Maternal_PID != "")

# --- 2) Fix Maternal IGRA values ("INT" -> "Pos") ---
infant_outcome <- infant_outcome %>%
  dplyr::mutate(
    Maternal_IGRA = dplyr::case_when(
      stringr::str_to_upper(trimws(Maternal_IGRA)) == "INT" ~ "Pos",
      TRUE ~ Maternal_IGRA
    )
  )

# --- 3) Define binary proliferator call based on BCG–MED delta ---
infant_outcome <- infant_outcome %>%
  dplyr::mutate(
    y_bin = dplyr::case_when(
      is.finite(y_delta) & y_delta > 5 ~ 1L,    # responder: BCG-induced Δ > 5%
      is.finite(y_delta) & y_delta <= 5 ~ 0L,   # non-responder: Δ ≤ 5%
      TRUE ~ NA_integer_
    )
  )
infant_outcome_12w <- infant_outcome %>%
  filter(Time_Point == "12 Wks") %>%
  select(-Time_Point)   # drop timepoint since all are 12W now

############# Maternal Features ################
maternal_plasma <- plasma_wide %>%
  dplyr::filter(Group == "Mother", TP == "Entry") %>%
  dplyr::mutate(PID = as.character(PID)) %>%
  dplyr::select(-Group, -TP)

maternal_plasma <- maternal_plasma %>%
  rename(Maternal_PID = PID)

# --- helpers ---
normalize_gatepath <- function(x) {
  x %>% stringr::str_replace_all("[\\u00A0\\p{Zs}]+", " ") %>% trimws()
}

suffix_after_ld_safe <- function(gp_norm) {
  core <- sub("\\s*\\|\\s*Freq\\..*$", "", as.character(gp_norm), perl = TRUE)
  suf  <- sub("^.*CD45,\\s*Live\\s*Dead\\s*subset/\\s*", "", core, perl = TRUE)
  if (identical(suf, core)) {
    suf <- sub("^.*Single\\s*Cells/\\s*", "", core, perl = TRUE)
  }
  trimws(suf)
}

format_suffix_for_feature <- function(suf) {
  s <- suf
  s <- gsub(",", "", s)                          # drop commas
  s <- gsub("\\s*/\\s*", "/", s)                 # tidy slashes
  s <- gsub("\\s+", " ", s)                      # collapse multiple spaces
  
  # preserve gdTCR vs gdTCR-
  s <- gsub("(?i)\\bgd\\s*TCR-\\b", "gdTCR-", s, perl = TRUE)  # negative gate
  s <- gsub("(?i)\\bgd\\s*TCR\\b",  "gdTCR",  s, perl = TRUE)  # positive gate
  
  s <- gsub("\\s+", "", s)                      # remove all spaces
  s <- gsub("/", "_", s)                        # slashes -> underscores
  s
}

baseline_lvls <- c("MED","Media","Unstim","UNSTIM","None","Baseline","CTRL","Control")

maternal_cytokine <- cytokine_long %>%
  mutate(
    PID_Type      = coalesce(PID_Type, `PID Type`),
    Time_Point    = coalesce(Time_Point, `Time Point`),
    GatePath_norm = normalize_gatepath(GatePath),
    GateSuffix    = suffix_after_ld_safe(GatePath_norm),
    freq_code     = case_when(
      Metric == "Freq. of Parent"      ~ "P",
      Metric == "Freq. of Grandparent" ~ "GP",
      TRUE                             ~ "UNK"
    )
  ) %>%
  filter(
    PID_Type   == "Mother",
    Time_Point == "Entry",
    !is.na(Condition),
    !Condition %in% baseline_lvls   # ⬅️ drop MED / baseline; keep all other stims
  ) %>%
  mutate(
    suf_fmt = format_suffix_for_feature(GateSuffix),
    feat    = paste0(Condition, "_", suf_fmt, "_", freq_code)
  ) %>%
  group_by(PID, feat) %>%
  summarise(
    value = dplyr::first(Delta_for_model),   # one ΔStim per mother×condition×gate
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = feat,
    values_from = value,
    values_fill = NA_real_
  )

master <- infant_outcome_12w %>%
  left_join(maternal_cytokine, by = c("Maternal_PID" = "PID"))

master <- master %>%
  left_join(maternal_plasma, by = "Maternal_PID")

ncol(master)
## ----------------------------------------------
## Multi-layer LASSO models for infant proliferation
## ----------------------------------------------


## 0) Output directory
out_root <- "C:/Users/ammas/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/Model_Results"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

## 1) Helper: identify plasma biomarker columns
##    (explicit list based on your plasma_wide structure)
plasma_markers <- c(
  "CD14","CD163","CD25","CRP","D-Dimer","GM-CSF","ICAM-1","iFABP",
  "IFN-a2","IFN-g","IL-10","IL-13","IL-17A","IL-18","IL-1a","IL-1b",
  "IL-1RA","IL-2","IL-22","IL-27","IL-4","IL-5","IL-6","IL-8","IP-10",
  "LBP","MCP-1","MMP-9","Neopterin","TNFa","TNFb","TNF-RI","TNF-RII",
  "VCAM-1"
)
## 2) Helper: general modeling function for one layer
run_layer_model <- function(
    df,
    layer_name,
    predictor_cols,
    y_col = "y_bin",
    id_col = "Infant_PID",
    out_root = out_root,
    missing_prop_max = 0.5,
    nzv_threshold = 1e-6
) {
  message("=== Layer: ", layer_name, " ===")
  
  layer_dir <- file.path(out_root, layer_name)
  dir.create(layer_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Keep only rows with non-missing outcome
  df_use <- df %>%
    filter(!is.na(.data[[y_col]]))
  
  # Keep only predictor cols that actually exist
  pred_cols <- intersect(predictor_cols, names(df_use))
  if (length(pred_cols) == 0) {
    warning("No predictor columns found for layer ", layer_name)
    return(invisible(NULL))
  }
  
  # For flow layers: restrict to rows that have at least one non-NA
  # in this layer (so infants without this stim are dropped)
  row_has_data <- df_use %>%
    select(all_of(pred_cols)) %>%
    mutate(row_na = apply(., 1, function(x) all(is.na(x)))) %>%
    pull(row_na)
  
  df_use <- df_use[!row_has_data, , drop = FALSE]
  
  if (nrow(df_use) < 20) {
    warning("Too few samples (", nrow(df_use), ") for layer ", layer_name, " – skipping.")
    return(invisible(NULL))
  }
  
  # Subset to predictors + outcome
  X_raw <- df_use %>%
    select(all_of(pred_cols))
  y <- df_use[[y_col]]
  ids <- df_use[[id_col]]
  
  # 2a) Drop predictors with too many missing values
  keep_by_missing <- colMeans(is.na(X_raw)) < missing_prop_max
  X_raw <- X_raw[, keep_by_missing, drop = FALSE]
  dropped_missing <- pred_cols[!keep_by_missing]
  
  # 2b) Drop near-zero variance predictors
  var_vec <- sapply(X_raw, function(z) stats::var(z, na.rm = TRUE))
  keep_by_var <- ifelse(is.na(var_vec), FALSE, var_vec > nzv_threshold)
  X_raw <- X_raw[, keep_by_var, drop = FALSE]
  dropped_nzv <- names(var_vec)[!keep_by_var]
  
  if (ncol(X_raw) == 0) {
    warning("All predictors dropped for layer ", layer_name)
    return(invisible(NULL))
  }
  
  # 2c) Median impute NAs
  X_imp <- X_raw %>%
    mutate(across(everything(), ~ {
      v <- .
      med <- median(v, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      v[is.na(v)] <- med
      v
    }))
  
  # 2d) Standardize
  X_scaled <- scale(as.matrix(X_imp))
  colnames(X_scaled) <- colnames(X_imp)
  
  # 3) Fit LASSO logistic with CV
  set.seed(1234)
  cvfit <- cv.glmnet(
    x = X_scaled,
    y = y,
    family = "binomial",
    alpha = 1,           # LASSO
    nfolds = min(10, length(unique(ids))), # up to 10 folds
    type.measure = "deviance"
  )
  
  # Choose lambda.1se (more conservative, better interpretability)
  lambda_use <- cvfit$lambda.1se
  
  coef_mat <- coef(cvfit, s = lambda_use)
  
  # coef_mat is a sparse matrix; turn into a clean tibble
  coef_df <- tibble::tibble(
    Feature     = rownames(coef_mat),
    Coefficient = as.numeric(coef_mat)
  ) %>%
    dplyr::mutate(
      OddsRatio = exp(Coefficient),
      Layer     = layer_name
    )
  
  
  # Separate intercept vs predictors
  coef_int <- coef_df %>% filter(Feature == "(Intercept)")
  coef_pred <- coef_df %>%
    filter(Feature != "(Intercept)", Coefficient != 0) %>%
    arrange(desc(abs(Coefficient)))
  
  # If no predictors selected, still save info
  coef_path <- file.path(layer_dir, paste0("coefficients_", layer_name, ".csv"))
  write_csv(coef_df, coef_path)
  
  # 4) ROC curve (using linear predictor at lambda_use)
  lp <- predict(cvfit, newx = X_scaled, s = lambda_use, type = "link")
  prob <- as.numeric(1 / (1 + exp(-lp)))
  roc_obj <- tryCatch(
    pROC::roc(response = y, predictor = prob),
    error = function(e) NULL
  )
  
  auc_val <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA_real_
  
  # Save a small summary text file
  summary_tbl <- tibble::tibble(
    Layer = layer_name,
    n_samples = length(y),
    n_predictors_input = length(pred_cols),
    n_dropped_missing = length(dropped_missing),
    n_dropped_nzv = length(dropped_nzv),
    n_predictors_final = ncol(X_scaled),
    n_selected_nonzero = nrow(coef_pred),
    lambda_use = lambda_use,
    AUC = auc_val
  )
  write_csv(summary_tbl, file.path(layer_dir, paste0("summary_", layer_name, ".csv")))
  
  # 5) Plots
  # 5a) Coefficient barplot
  if (nrow(coef_pred) > 0) {
    coef_pred$Feature <- factor(
      coef_pred$Feature,
      levels = coef_pred$Feature[order(coef_pred$Coefficient)]
    )
    
    p_coef <- ggplot(coef_pred, aes(x = Feature, y = Coefficient)) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      labs(
        title = paste0("LASSO coefficients — ", layer_name),
        subtitle = paste0("Non-zero predictors at lambda.1se (AUC ≈ ",
                          ifelse(is.na(auc_val), "NA", round(auc_val, 3)), ")"),
        x = NULL,
        y = "Log-odds coefficient (higher → higher log-odds of responder)"
      )
    
    ggsave(
      filename = file.path(layer_dir, paste0("coef_plot_", layer_name, ".png")),
      plot = p_coef,
      width = 7, height = max(3, 0.2 * nrow(coef_pred) + 2),
      dpi = 300, bg = "white"
    )
  }
  
  # 5b) ROC plot
  if (!is.null(roc_obj)) {
    roc_df <- tibble::tibble(
      TPR = rev(roc_obj$sensitivities),
      FPR = rev(1 - roc_obj$specificities)
    )
    
    p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
      geom_line() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
      theme_bw() +
      coord_equal() +
      labs(
        title = paste0("ROC — ", layer_name),
        subtitle = paste0("AUC = ", round(auc_val, 3)),
        x = "1 - Specificity",
        y = "Sensitivity"
      )
    
    ggsave(
      filename = file.path(layer_dir, paste0("roc_", layer_name, ".png")),
      plot = p_roc,
      width = 5, height = 5,
      dpi = 300, bg = "white"
    )
  }
  
  # Return key objects for inspection in R
  list(
    layer = layer_name,
    summary = summary_tbl,
    coef_full = coef_df,
    coef_selected = coef_pred,
    cvfit = cvfit,
    auc = auc_val
  )
}

## 3) Define predictor sets per layer

# Plasma: maternal plasma biomarkers only
plasma_cols <- intersect(plasma_markers, names(master))

# Flow layers by prefix
bcg_cols   <- names(master)[startsWith(names(master), "BCG_")]
dosr_cols  <- names(master)[startsWith(names(master), "DosR_")]
e6c10_cols <- names(master)[startsWith(names(master), "E6C10_")]
gag_cols   <- names(master)[startsWith(names(master), "GAG_")]

## 4) Run models per layer

results <- list()

# 4.1 Plasma
results$Plasma <- run_layer_model(
  df = master,
  layer_name = "Plasma",
  predictor_cols = plasma_cols,
  out_root = out_root
)

# 4.2 BCG
results$BCG <- run_layer_model(
  df = master,
  layer_name = "BCG",
  predictor_cols = bcg_cols,
  out_root = out_root
)

# 4.3 DosR
results$DosR <- run_layer_model(
  df = master,
  layer_name = "DosR",
  predictor_cols = dosr_cols,
  out_root = out_root
)

# 4.4 E6C10
results$E6C10 <- run_layer_model(
  df = master,
  layer_name = "E6C10",
  predictor_cols = e6c10_cols,
  out_root = out_root
)

# 4.5 GAG
results$GAG <- run_layer_model(
  df = master,
  layer_name = "GAG",
  predictor_cols = gag_cols,
  out_root = out_root
)

## 5) Save a combined summary of all layers
all_summaries <- purrr::map_dfr(
  results[!sapply(results, is.null)],
  ~ .x$summary
)
write_csv(all_summaries, file.path(out_root, "summary_all_layers.csv"))
