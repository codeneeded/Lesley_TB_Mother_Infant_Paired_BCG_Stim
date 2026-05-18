## ============================================================================
## Project: Flow Cytometry and Plasma Biomarkers on Mothers and Infants (P1078)
## Script:  4_LASSO_Models.R
## Purpose: Unified LASSO logistic regression modeling against two outcomes,
##          using maternal and infant flow + plasma feature layers.
##
## OUTCOME A — Infant 12wk BCG response
##   Definition: BCG−MED Δ in (gd TCR−/CD8/Proliferating, Freq. of Parent) > 5 pp
##   Predictors: maternal features ONLY (plasma at Entry, flow per stim at Entry).
##               Infant features are excluded because the outcome IS infant flow
##               at 12wk; predicting it from contemporaneous infant features
##               would be circular.
##
## OUTCOME B — Infant IGRA status (Pos vs Neg)
##   Definition: y = 1 if Infant IGRA == "Pos", 0 if "Neg", excluded otherwise.
##   Predictors: ALL available feature layers — maternal flow/plasma at Entry,
##               infant flow/plasma at 12wk and 44wk, plus a combined all-flow
##               layer. IGRA is measured independently of flow/plasma so
##               cross-source prediction is meaningful as a biomarker
##               discovery exercise (not a causal claim).
##
## Modeling choices (applied to every layer):
##   - LASSO logistic regression with stratified k-fold CV (folds balanced by
##     class), default k = 10 capped at min(n_pos, n_neg).
##   - cv.glmnet with type.measure = "auc" so the reported CV_AUC_cvglmnet
##     value is honestly cross-validated.
##   - Out-of-fold predictions at lambda.1se are concatenated for an OOF ROC
##     curve and OOF AUC (also honestly held out).
##   - Predictors with > 50% missingness or near-zero variance are dropped;
##     remaining NAs are median-imputed; features are standardised before fit.
##
## Justification for the 5pp threshold (Outcome A): unchanged from prior
## version. The 95th percentile of MED-control proliferation is < 0.5%, so
## 5 pp is ~10× the noise floor; abs ≥ 2% would classify all 35 infants as
## responders because infant baseline CD8 proliferation is already > 2%.
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(readr); library(purrr)
  library(digest); library(ggplot2); library(glmnet); library(pROC); library(tibble)
  library(qs2)
})

## ---- Paths ----------------------------------------------------------------
proj_root <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/Flow Cytometry and Plasma Biomarkers on Mothers and Infants (P1078)"
saved_dir <- file.path(proj_root, "saved_R_dat")
out_root  <- file.path(proj_root, "Results", "LASSO_Models")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

## ---- Load preprocessed inputs ---------------------------------------------
cyto_qs2 <- file.path(saved_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.qs2")
cyto_csv <- file.path(saved_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.csv")
cytokine_long <- if (file.exists(cyto_qs2)) qs2::qs_read(cyto_qs2) else read_csv(cyto_csv, show_col_types = FALSE)

plasma_qs2 <- file.path(saved_dir, "Plasma_Biomarker_wide_scaled.qs2")
plasma_csv <- file.path(saved_dir, "Plasma_Biomarker_wide_scaled.csv")
plasma_wide <- if (file.exists(plasma_qs2)) qs2::qs_read(plasma_qs2) else read_csv(plasma_csv, show_col_types = FALSE)

## ---- Helpers --------------------------------------------------------------
norm_tp <- function(x) {
  z <- tolower(trimws(as.character(x)))
  case_when(
    z %in% c("entry") ~ "Entry",
    z %in% c("dx","diagnosis") ~ "Dx",
    z %in% c("12","12 wks","12 weeks","12wks","12wk") ~ "12 Wks",
    z %in% c("44","44 wks","44 weeks","44wks","44wk") ~ "44 Wks",
    TRUE ~ as.character(x)
  )
}

normalize_gatepath <- function(x) {
  x %>% as.character() %>%
    str_replace_all("[\\u00A0\\p{Zs}]+", " ") %>%
    str_replace_all("\\s*,\\s*", ", ") %>%
    str_replace_all("\\s*\\|\\s*", " | ") %>%
    str_squish()
}

suffix_after_ld_safe <- function(gp_norm) {
  core <- sub("\\s*\\|\\s*Freq\\..*$", "", as.character(gp_norm), perl = TRUE)
  suf  <- sub("^.*CD45,\\s*Live\\s*Dead\\s*subset/\\s*", "", core, perl = TRUE)
  if (identical(suf, core)) suf <- sub("^.*Single\\s*Cells/\\s*", "", core, perl = TRUE)
  trimws(suf)
}

format_suffix_for_feature <- function(suf) {
  s <- suf
  s <- gsub(",", "", s)
  s <- gsub("\\s*/\\s*", "/", s)
  s <- gsub("\\s+", " ", s)
  s <- gsub("(?i)\\bgd\\s*TCR-\\b", "gdTCR-", s, perl = TRUE)
  s <- gsub("(?i)\\bgd\\s*TCR\\b",  "gdTCR",  s, perl = TRUE)
  s <- gsub("\\s+", "", s)
  s <- gsub("/", "_", s)
  s
}

# Stratified fold-id constructor (balanced positives/negatives per fold)
make_stratified_folds <- function(y, k, seed = 1234) {
  set.seed(seed)
  fid <- integer(length(y))
  for (cls in unique(y)) {
    idx <- which(y == cls)
    fid[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  fid
}

## ---- Normalize cytokine columns & derive suffix ---------------------------
cytokine_long <- cytokine_long %>%
  mutate(
    PID           = as.character(PID),
    PID_Type      = coalesce(PID_Type, `PID Type`),
    Time_Point    = norm_tp(coalesce(Time_Point, `Time Point`)),
    GatePath_norm = normalize_gatepath(GatePath),
    GateSuffix    = suffix_after_ld_safe(GatePath_norm)
  )

## ============================================================================
## A) BUILD FEATURE WIDE-TABLES (one row per PID; columns per feature)
## ============================================================================

# All proliferating cytokine features (excludes Non-Proliferating gates).
# Uses Delta_for_model (Stim − MED) so features are stimulation-specific.
baseline_lvls <- c("MED","Media","Unstim","UNSTIM","None","Baseline","CTRL","Control")

flow_long_prolif <- cytokine_long %>%
  filter(
    !is.na(Condition), !Condition %in% baseline_lvls,
    Metric %in% c("Freq. of Parent", "Freq. of Grandparent"),
    str_detect(GateSuffix, "Proliferating"),
    !str_detect(GateSuffix, "Non-Proliferating"),
    !is.na(Delta_for_model)
  ) %>%
  mutate(
    freq_code = case_when(
      Metric == "Freq. of Parent"      ~ "P",
      Metric == "Freq. of Grandparent" ~ "GP",
      TRUE                              ~ "UNK"
    ),
    suf_fmt = format_suffix_for_feature(GateSuffix),
    feat_core = paste0(Condition, "_", suf_fmt, "_", freq_code)
  )

# Helper to make a wide table for one (PID_Type, Time_Point) slice
make_flow_wide <- function(df_long, ptype, tp, prefix) {
  d <- df_long %>%
    filter(PID_Type == ptype, Time_Point == tp)
  
  # Detect duplicate (PID, feat_core) keys before collapse; warn but proceed
  dup <- d %>% group_by(PID, feat_core) %>% filter(n() > 1) %>% ungroup()
  if (nrow(dup) > 0) {
    dup_summary <- dup %>% count(PID, feat_core, name = "n_rows") %>% arrange(desc(n_rows))
    warning(sprintf(
      "[%s|%s] %d duplicate (PID, feat) rows across %d keys; keeping first(). Top:\n%s",
      ptype, tp, nrow(dup), nrow(dup_summary),
      paste(utils::capture.output(print(utils::head(dup_summary, 6))), collapse = "\n")
    ))
  }
  
  d %>%
    group_by(PID, feat_core) %>%
    summarise(value = dplyr::first(Delta_for_model), .groups = "drop") %>%
    mutate(feat = paste0(prefix, "_", feat_core)) %>%
    select(PID, feat, value) %>%
    pivot_wider(names_from = feat, values_from = value, values_fill = NA_real_)
}

# Maternal flow @ Entry
maternal_flow_entry <- make_flow_wide(flow_long_prolif, "Mother", "Entry", "MatFlow") %>%
  rename(Maternal_PID = PID)

# Infant flow @ 12 wks
infant_flow_12wk <- make_flow_wide(flow_long_prolif, "Infant", "12 Wks", "Inf12Flow") %>%
  rename(Infant_PID = PID)

# Infant flow @ 44 wks
infant_flow_44wk <- make_flow_wide(flow_long_prolif, "Infant", "44 Wks", "Inf44Flow") %>%
  rename(Infant_PID = PID)

# Plasma wide is already z-scored within TP from Script 3.
# Split by Group × TP and prefix to disambiguate maternal/infant timepoints.
plasma_marker_cols <- setdiff(names(plasma_wide), c("PID","Group","TP"))

maternal_plasma_entry <- plasma_wide %>%
  filter(Group == "Mother", TP == "Entry") %>%
  mutate(PID = as.character(PID)) %>%
  select(PID, all_of(plasma_marker_cols)) %>%
  rename_with(~ paste0("MatPlasma_", .x), all_of(plasma_marker_cols)) %>%
  rename(Maternal_PID = PID)

infant_plasma_12wk <- plasma_wide %>%
  filter(Group == "Infant", TP == "12 Wks") %>%
  mutate(PID = as.character(PID)) %>%
  select(PID, all_of(plasma_marker_cols)) %>%
  rename_with(~ paste0("Inf12Plasma_", .x), all_of(plasma_marker_cols)) %>%
  rename(Infant_PID = PID)

infant_plasma_44wk <- plasma_wide %>%
  filter(Group == "Infant", TP == "44 Wks") %>%
  mutate(PID = as.character(PID)) %>%
  select(PID, all_of(plasma_marker_cols)) %>%
  rename_with(~ paste0("Inf44Plasma_", .x), all_of(plasma_marker_cols)) %>%
  rename(Infant_PID = PID)

message(sprintf(
  "Feature tables built — MatFlow:%d cols, MatPlasma:%d cols, Inf12Flow:%d cols, Inf12Plasma:%d cols, Inf44Flow:%d cols, Inf44Plasma:%d cols.",
  ncol(maternal_flow_entry)-1, ncol(maternal_plasma_entry)-1,
  ncol(infant_flow_12wk)-1, ncol(infant_plasma_12wk)-1,
  ncol(infant_flow_44wk)-1, ncol(infant_plasma_44wk)-1
))

## ============================================================================
## B) BUILD OUTCOMES
## ============================================================================

## Outcome A — Infant 12wk BCG response
##   y_bin = 1 iff BCG−MED Δ in (gd TCR−/CD8/Proliferating, Freq. of Parent) > 5 pp
rx_suffix <- "gd TCR-/CD8/Proliferating$"

infant_outcome_bcg <- cytokine_long %>%
  filter(
    PID_Type   == "Infant",
    Time_Point == "12 Wks",
    Condition  == "BCG",
    str_detect(GateSuffix, rx_suffix),
    str_detect(GatePath_norm, "\\|\\s*Freq\\.\\s*of\\s*Parent")
  ) %>%
  arrange(PID, desc(GateDepth)) %>%
  distinct(PID, .keep_all = TRUE) %>%
  transmute(
    Infant_PID    = PID,
    Maternal_PID  = as.character(`Maternal PID`),
    Batch_infant  = as.character(Batch),
    Infant_IGRA   = `Infant IGRA`,
    Maternal_IGRA = `Maternal IGRA`,
    y_abs         = as.numeric(Value_RawNormalized),
    y_delta       = as.numeric(Delta_for_model)
  ) %>%
  filter(!is.na(Maternal_PID) & Maternal_PID != "") %>%
  mutate(
    # IGRA semantics: INT → Pos; ATB stays its own category (not merged with Pos).
    Maternal_IGRA = case_when(
      str_to_upper(trimws(Maternal_IGRA)) == "INT" ~ "Pos",
      TRUE                                          ~ Maternal_IGRA
    ),
    y_bin = case_when(
      is.finite(y_delta) & y_delta >  5 ~ 1L,
      is.finite(y_delta) & y_delta <= 5 ~ 0L,
      TRUE                               ~ NA_integer_
    )
  )

message(sprintf(
  "Outcome A (12wk BCG response): n=%d total | y=1: %d | y=0: %d | NA: %d",
  nrow(infant_outcome_bcg),
  sum(infant_outcome_bcg$y_bin == 1, na.rm = TRUE),
  sum(infant_outcome_bcg$y_bin == 0, na.rm = TRUE),
  sum(is.na(infant_outcome_bcg$y_bin))
))

## Outcome B — Infant IGRA status (Pos vs Neg)
##   Pos → 1, Neg → 0, everything else excluded.
infant_outcome_igra <- cytokine_long %>%
  filter(PID_Type == "Infant") %>%
  distinct(
    Infant_PID    = PID,
    Maternal_PID  = `Maternal PID`,
    Infant_IGRA   = `Infant IGRA`,
    Maternal_IGRA = `Maternal IGRA`
  ) %>%
  mutate(
    Infant_PID    = as.character(Infant_PID),
    Maternal_PID  = as.character(Maternal_PID),
    Infant_IGRA   = na_if(trimws(as.character(Infant_IGRA)), ""),
    Maternal_IGRA = case_when(
      str_to_upper(trimws(Maternal_IGRA)) == "INT" ~ "Pos",
      TRUE                                          ~ Maternal_IGRA
    ),
    y_bin = case_when(
      Infant_IGRA == "Pos" ~ 1L,
      Infant_IGRA == "Neg" ~ 0L,
      TRUE                 ~ NA_integer_
    )
  ) %>%
  filter(!is.na(y_bin)) %>%
  # one row per infant (some infants appear at multiple TPs; IGRA is per-infant)
  group_by(Infant_PID) %>%
  arrange(desc(!is.na(Maternal_PID))) %>%
  slice(1) %>%
  ungroup()

message(sprintf(
  "Outcome B (Infant IGRA): n=%d | y=1 (Pos): %d | y=0 (Neg): %d",
  nrow(infant_outcome_igra),
  sum(infant_outcome_igra$y_bin == 1),
  sum(infant_outcome_igra$y_bin == 0)
))

## ============================================================================
## C) BUILD MASTER TABLES (outcomes + features)
## ============================================================================

## Master A — Outcome A with maternal-only features
master_bcg <- infant_outcome_bcg %>%
  left_join(maternal_flow_entry,   by = "Maternal_PID") %>%
  left_join(maternal_plasma_entry, by = "Maternal_PID")

message("master_bcg: ", nrow(master_bcg), " rows, ", ncol(master_bcg), " cols")

## Master B — Outcome B with all features (maternal + infant 12wk + infant 44wk)
master_igra <- infant_outcome_igra %>%
  left_join(maternal_flow_entry,   by = "Maternal_PID") %>%
  left_join(maternal_plasma_entry, by = "Maternal_PID") %>%
  left_join(infant_flow_12wk,      by = "Infant_PID") %>%
  left_join(infant_plasma_12wk,    by = "Infant_PID") %>%
  left_join(infant_flow_44wk,      by = "Infant_PID") %>%
  left_join(infant_plasma_44wk,    by = "Infant_PID")

message("master_igra: ", nrow(master_igra), " rows, ", ncol(master_igra), " cols")

## ============================================================================
## D) MODEL-RUNNER FUNCTION (one layer at a time)
## ============================================================================

run_layer_model <- function(
    df, layer_name, predictor_cols,
    out_root,
    y_col = "y_bin",
    id_col = NULL,
    missing_prop_max = 0.5,
    nzv_threshold = 1e-6,
    k_folds = 10
) {
  message("=== Layer: ", layer_name, " ===")
  layer_dir <- file.path(out_root, layer_name)
  dir.create(layer_dir, recursive = TRUE, showWarnings = FALSE)
  
  df_use <- df %>% filter(!is.na(.data[[y_col]]))
  pred_cols <- intersect(predictor_cols, names(df_use))
  if (length(pred_cols) == 0) {
    warning("No predictor columns found for layer ", layer_name); return(invisible(NULL))
  }
  
  # Drop rows where every predictor is NA
  row_all_na <- df_use %>% select(all_of(pred_cols)) %>%
    apply(1, function(x) all(is.na(x)))
  df_use <- df_use[!row_all_na, , drop = FALSE]
  if (nrow(df_use) < 20) {
    warning("Too few samples (", nrow(df_use), ") for layer ", layer_name, " — skipping.")
    return(invisible(NULL))
  }
  
  X_raw <- df_use %>% select(all_of(pred_cols))
  y     <- df_use[[y_col]]
  
  # Drop high-missing and low-variance predictors
  keep_by_missing <- colMeans(is.na(X_raw)) < missing_prop_max
  X_raw <- X_raw[, keep_by_missing, drop = FALSE]
  dropped_missing <- pred_cols[!keep_by_missing]
  
  var_vec <- sapply(X_raw, function(z) stats::var(z, na.rm = TRUE))
  keep_by_var <- ifelse(is.na(var_vec), FALSE, var_vec > nzv_threshold)
  X_raw <- X_raw[, keep_by_var, drop = FALSE]
  dropped_nzv <- names(var_vec)[!keep_by_var]
  if (ncol(X_raw) == 0) {
    warning("All predictors dropped for layer ", layer_name); return(invisible(NULL))
  }
  
  # Median impute, then standardize
  X_imp <- X_raw %>% mutate(across(everything(), ~ {
    v <- .
    med <- stats::median(v, na.rm = TRUE); if (!is.finite(med)) med <- 0
    v[is.na(v)] <- med; v
  }))
  X_scaled <- scale(as.matrix(X_imp)); colnames(X_scaled) <- colnames(X_imp)
  
  # Stratified folds
  n_pos <- sum(y == 1); n_neg <- sum(y == 0)
  if (n_pos < 3 || n_neg < 3) {
    warning("Layer ", layer_name, ": only ", n_pos, " positives and ", n_neg,
            " negatives — need ≥3 of each for stratified CV. Skipping.")
    return(invisible(NULL))
  }
  k_use   <- min(k_folds, n_pos, n_neg)
  foldids <- make_stratified_folds(y, k = k_use, seed = 1234)
  
  # Fit cv.glmnet with AUC as the CV measure
  cvfit <- tryCatch(
    cv.glmnet(x = X_scaled, y = y, family = "binomial", alpha = 1,
              foldid = foldids, type.measure = "auc"),
    error = function(e) {
      warning("cv.glmnet AUC failed for ", layer_name, ": ", conditionMessage(e),
              " — falling back to deviance.")
      cv.glmnet(x = X_scaled, y = y, family = "binomial", alpha = 1,
                foldid = foldids, type.measure = "deviance")
    }
  )
  lambda_use <- cvfit$lambda.1se
  
  cv_auc_at_l1se <- if (!is.null(cvfit$name) && grepl("AUC", cvfit$name, ignore.case = TRUE)) {
    as.numeric(cvfit$cvm[which(cvfit$lambda == lambda_use)])
  } else NA_real_
  
  # Out-of-fold predictions for honest ROC
  oof_prob <- rep(NA_real_, length(y))
  for (k in sort(unique(foldids))) {
    tr <- which(foldids != k); te <- which(foldids == k)
    if (length(unique(y[tr])) < 2) next
    fit_k <- glmnet(x = X_scaled[tr, , drop = FALSE], y = y[tr],
                    family = "binomial", alpha = 1, lambda = lambda_use)
    lp_k <- predict(fit_k, newx = X_scaled[te, , drop = FALSE], s = lambda_use, type = "link")
    oof_prob[te] <- as.numeric(1 / (1 + exp(-lp_k)))
  }
  
  roc_obj <- tryCatch(
    pROC::roc(response = y[!is.na(oof_prob)], predictor = oof_prob[!is.na(oof_prob)], quiet = TRUE),
    error = function(e) NULL
  )
  oof_auc <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA_real_
  
  coef_mat <- coef(cvfit, s = lambda_use)
  coef_df <- tibble(Feature = rownames(coef_mat),
                    Coefficient = as.numeric(coef_mat)) %>%
    mutate(OddsRatio = exp(Coefficient), Layer = layer_name)
  coef_pred <- coef_df %>%
    filter(Feature != "(Intercept)", Coefficient != 0) %>%
    arrange(desc(abs(Coefficient)))
  
  write_csv(coef_df, file.path(layer_dir, paste0("coefficients_", layer_name, ".csv")))
  
  summary_tbl <- tibble(
    Layer              = layer_name,
    n_samples          = length(y),
    n_pos              = n_pos,
    n_neg              = n_neg,
    n_folds            = k_use,
    n_predictors_input = length(pred_cols),
    n_dropped_missing  = length(dropped_missing),
    n_dropped_nzv      = length(dropped_nzv),
    n_predictors_final = ncol(X_scaled),
    n_selected_nonzero = nrow(coef_pred),
    lambda_use         = lambda_use,
    CV_AUC_cvglmnet    = cv_auc_at_l1se,
    OOF_AUC_roc        = oof_auc
  )
  write_csv(summary_tbl, file.path(layer_dir, paste0("summary_", layer_name, ".csv")))
  
  # Coefficient barplot
  if (nrow(coef_pred) > 0) {
    coef_pred$Feature <- factor(coef_pred$Feature,
                                levels = coef_pred$Feature[order(coef_pred$Coefficient)])
    auc_label <- if (is.finite(oof_auc)) round(oof_auc, 3) else "NA"
    p_coef <- ggplot(coef_pred, aes(x = Feature, y = Coefficient)) +
      geom_col() + coord_flip() + theme_bw() +
      labs(title = paste0("LASSO coefficients — ", layer_name),
           subtitle = paste0("Non-zero predictors at lambda.1se (OOF AUC = ", auc_label, ")"),
           x = NULL,
           y = "Log-odds coefficient (higher → higher log-odds of y=1)")
    ggsave(file.path(layer_dir, paste0("coef_plot_", layer_name, ".png")),
           p_coef, width = 7, height = max(3, 0.2 * nrow(coef_pred) + 2),
           dpi = 300, bg = "white")
  }
  
  # ROC plot
  if (!is.null(roc_obj)) {
    roc_df <- tibble(TPR = rev(roc_obj$sensitivities),
                     FPR = rev(1 - roc_obj$specificities))
    p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
      geom_line() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
      theme_bw() + coord_equal() +
      labs(title = paste0("Out-of-fold ROC — ", layer_name),
           subtitle = paste0("OOF AUC = ", round(oof_auc, 3),
                             " | CV AUC (cv.glmnet) = ",
                             ifelse(is.finite(cv_auc_at_l1se), round(cv_auc_at_l1se, 3), "NA")),
           x = "1 - Specificity", y = "Sensitivity")
    ggsave(file.path(layer_dir, paste0("roc_", layer_name, ".png")),
           p_roc, width = 5, height = 5, dpi = 300, bg = "white")
  }
  
  list(layer = layer_name, summary = summary_tbl,
       coef_full = coef_df, coef_selected = coef_pred,
       cvfit = cvfit, oof_prob = oof_prob, foldids = foldids,
       cv_auc = cv_auc_at_l1se, oof_auc = oof_auc)
}

## ============================================================================
## E) PREDICTOR SET DEFINITIONS
## ============================================================================
nm <- function(df) names(df)
cols_with <- function(df, prefix) nm(df)[startsWith(nm(df), prefix)]

# Outcome A — maternal only
A_plasma_cols   <- cols_with(master_bcg, "MatPlasma_")
A_mat_bcg_cols  <- cols_with(master_bcg, "MatFlow_BCG_")
A_mat_dosr_cols <- cols_with(master_bcg, "MatFlow_DosR_")
A_mat_e6c10_cols<- cols_with(master_bcg, "MatFlow_E6C10_")
A_mat_gag_cols  <- cols_with(master_bcg, "MatFlow_GAG_")

# Outcome B — maternal + infant
B_mat_plasma     <- cols_with(master_igra, "MatPlasma_")
B_inf12_plasma   <- cols_with(master_igra, "Inf12Plasma_")
B_inf44_plasma   <- cols_with(master_igra, "Inf44Plasma_")

B_mat_bcg        <- cols_with(master_igra, "MatFlow_BCG_")
B_mat_dosr       <- cols_with(master_igra, "MatFlow_DosR_")
B_mat_e6c10      <- cols_with(master_igra, "MatFlow_E6C10_")
B_mat_gag        <- cols_with(master_igra, "MatFlow_GAG_")

B_inf12_bcg      <- cols_with(master_igra, "Inf12Flow_BCG_")
B_inf12_dosr     <- cols_with(master_igra, "Inf12Flow_DosR_")
B_inf12_e6c10    <- cols_with(master_igra, "Inf12Flow_E6C10_")
B_inf12_gag      <- cols_with(master_igra, "Inf12Flow_GAG_")

B_inf44_bcg      <- cols_with(master_igra, "Inf44Flow_BCG_")
B_inf44_dosr     <- cols_with(master_igra, "Inf44Flow_DosR_")
B_inf44_e6c10    <- cols_with(master_igra, "Inf44Flow_E6C10_")
B_inf44_gag      <- cols_with(master_igra, "Inf44Flow_GAG_")

B_all_flow       <- c(cols_with(master_igra, "MatFlow_"),
                      cols_with(master_igra, "Inf12Flow_"),
                      cols_with(master_igra, "Inf44Flow_"))
B_all_plasma     <- c(cols_with(master_igra, "MatPlasma_"),
                      cols_with(master_igra, "Inf12Plasma_"),
                      cols_with(master_igra, "Inf44Plasma_"))

## ============================================================================
## F) RUN MODELS
## ============================================================================

## -------- Outcome A: Infant 12wk BCG response (maternal predictors) --------
A_root <- file.path(out_root, "Infant_BCG_Response_12wk")
dir.create(A_root, recursive = TRUE, showWarnings = FALSE)

resA <- list()
resA$Plasma_Maternal_Entry      <- run_layer_model(master_bcg, "Plasma_Maternal_Entry",      A_plasma_cols,   out_root = A_root)
resA$Flow_Maternal_Entry_BCG    <- run_layer_model(master_bcg, "Flow_Maternal_Entry_BCG",    A_mat_bcg_cols,  out_root = A_root)
resA$Flow_Maternal_Entry_DosR   <- run_layer_model(master_bcg, "Flow_Maternal_Entry_DosR",   A_mat_dosr_cols, out_root = A_root)
resA$Flow_Maternal_Entry_E6C10  <- run_layer_model(master_bcg, "Flow_Maternal_Entry_E6C10",  A_mat_e6c10_cols,out_root = A_root)
resA$Flow_Maternal_Entry_GAG    <- run_layer_model(master_bcg, "Flow_Maternal_Entry_GAG",    A_mat_gag_cols,  out_root = A_root)

summaryA <- map_dfr(resA[!sapply(resA, is.null)], ~ .x$summary)
write_csv(summaryA, file.path(A_root, "summary_all_layers.csv"))

qs2::qs_save(master_bcg, file.path(A_root, "master_analysis_ready.qs2"))
write_csv(
  master_bcg %>% mutate(across(where(is.logical), as.integer),
                        across(where(is.factor),  as.character)),
  file.path(A_root, "master_analysis_ready.csv")
)

## -------- Outcome B: Infant IGRA status (all predictors) --------
B_root <- file.path(out_root, "Infant_IGRA_Status")
dir.create(B_root, recursive = TRUE, showWarnings = FALSE)

resB <- list()
# Plasma layers
resB$Plasma_Maternal_Entry      <- run_layer_model(master_igra, "Plasma_Maternal_Entry",      B_mat_plasma,   out_root = B_root)
resB$Plasma_Infant_12wk         <- run_layer_model(master_igra, "Plasma_Infant_12wk",         B_inf12_plasma, out_root = B_root)
resB$Plasma_Infant_44wk         <- run_layer_model(master_igra, "Plasma_Infant_44wk",         B_inf44_plasma, out_root = B_root)
resB$Plasma_All_Combined        <- run_layer_model(master_igra, "Plasma_All_Combined",        B_all_plasma,   out_root = B_root)

# Maternal flow layers
resB$Flow_Maternal_Entry_BCG    <- run_layer_model(master_igra, "Flow_Maternal_Entry_BCG",    B_mat_bcg,      out_root = B_root)
resB$Flow_Maternal_Entry_DosR   <- run_layer_model(master_igra, "Flow_Maternal_Entry_DosR",   B_mat_dosr,     out_root = B_root)
resB$Flow_Maternal_Entry_E6C10  <- run_layer_model(master_igra, "Flow_Maternal_Entry_E6C10",  B_mat_e6c10,    out_root = B_root)
resB$Flow_Maternal_Entry_GAG    <- run_layer_model(master_igra, "Flow_Maternal_Entry_GAG",    B_mat_gag,      out_root = B_root)

# Infant 12wk flow layers
resB$Flow_Infant_12wk_BCG       <- run_layer_model(master_igra, "Flow_Infant_12wk_BCG",       B_inf12_bcg,    out_root = B_root)
resB$Flow_Infant_12wk_DosR      <- run_layer_model(master_igra, "Flow_Infant_12wk_DosR",      B_inf12_dosr,   out_root = B_root)
resB$Flow_Infant_12wk_E6C10     <- run_layer_model(master_igra, "Flow_Infant_12wk_E6C10",     B_inf12_e6c10,  out_root = B_root)
resB$Flow_Infant_12wk_GAG       <- run_layer_model(master_igra, "Flow_Infant_12wk_GAG",       B_inf12_gag,    out_root = B_root)

# Infant 44wk flow layers
resB$Flow_Infant_44wk_BCG       <- run_layer_model(master_igra, "Flow_Infant_44wk_BCG",       B_inf44_bcg,    out_root = B_root)
resB$Flow_Infant_44wk_DosR      <- run_layer_model(master_igra, "Flow_Infant_44wk_DosR",      B_inf44_dosr,   out_root = B_root)
resB$Flow_Infant_44wk_E6C10     <- run_layer_model(master_igra, "Flow_Infant_44wk_E6C10",     B_inf44_e6c10,  out_root = B_root)
resB$Flow_Infant_44wk_GAG       <- run_layer_model(master_igra, "Flow_Infant_44wk_GAG",       B_inf44_gag,    out_root = B_root)

# Combined-flow layer
resB$Flow_All_Proliferating     <- run_layer_model(master_igra, "Flow_All_Proliferating",     B_all_flow,     out_root = B_root)

summaryB <- map_dfr(resB[!sapply(resB, is.null)], ~ .x$summary)
write_csv(summaryB, file.path(B_root, "summary_all_layers.csv"))

qs2::qs_save(master_igra, file.path(B_root, "master_analysis_ready.qs2"))
write_csv(
  master_igra %>% mutate(across(where(is.logical), as.integer),
                         across(where(is.factor),  as.character)),
  file.path(B_root, "master_analysis_ready.csv")
)

## ---- Combined top-level summary across BOTH outcomes ----------------------
summary_combined <- bind_rows(
  summaryA %>% mutate(Outcome = "Infant_BCG_Response_12wk"),
  summaryB %>% mutate(Outcome = "Infant_IGRA_Status")
) %>%
  select(Outcome, everything())

write_csv(summary_combined, file.path(out_root, "summary_all_models.csv"))

message("Script 4 complete. Model results in: ", out_root)

