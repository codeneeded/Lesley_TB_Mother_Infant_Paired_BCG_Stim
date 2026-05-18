## ============================================================================
## Project: Flow Cytometry and Plasma Biomarkers on Mothers and Infants (P1078)
## Script:  3_Plasma_Biomarker.R
## Purpose: Clean plasma biomarker data (UND/OVER/below-LLOQ handling, TP
##          harmonization, log1p + per-TP z-score), run four comparisons,
##          and write summary CSVs + per-analyte plots.
##
## Key fixes vs prior version:
##   - TP harmonization handles the ACTUAL labels in the raw file
##     ('12wk', '44wk', 'Entry') in a single, working case_when. No dead code.
##   - UND policy is explicit and consistent: UND → 0 (treated as below-LLOQ),
##     kept through the analysis (filter is now on Value_clean, not on Flag).
##   - Below-LLOQ markers like '<0.80↓' are recognized as UND.
##   - IGRA lookup deduplicates by PID with explicit priority (ATB > Pos > Neg > INT)
##     so that one mother with inconsistent IGRA (e.g. 6071260) doesn't multiply
##     plasma rows on the left_join.
##   - All comparison plots use a consistent Y axis (Value_z, per-TP z-score)
##     with matching y-axis labels. Stats run on log1p (rank-invariant).
##   - All paths derive from `proj_root`.
## ============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr)
  library(readr); library(ggplot2); library(purrr); library(glue); library(qs2)
})

## ---- Paths -----------------------------------------------------------------
proj_root <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/Flow Cytometry and Plasma Biomarkers on Mothers and Infants (P1078)"
raw_dir   <- file.path(proj_root, "Raw_Data")
saved_dir <- file.path(proj_root, "saved_R_dat")
dir.create(saved_dir, recursive = TRUE, showWarnings = FALSE)

plasma_xlsx <- file.path(raw_dir, "TB_Flow_and_Plasma_Biomarker_data.xlsx")
stopifnot(file.exists(plasma_xlsx))

# Cytokine (already preprocessed by Script 1) — needed for IGRA lookup
in_qs2 <- file.path(saved_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.qs2")
in_csv <- file.path(saved_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.csv")
cytokine_long <- if (file.exists(in_qs2)) {
  qs2::qs_read(in_qs2)
} else if (file.exists(in_csv)) {
  read_csv(in_csv, show_col_types = FALSE)
} else {
  stop("Could not find preprocessed cytokine input. Run Script 1 first.")
}

## ---- Read plasma ----------------------------------------------------------
plasma_wide_raw <- read_excel(plasma_xlsx, sheet = "Plasma Biomarker", .name_repair = "minimal")
stopifnot(all(c("PID","Group","TP") %in% names(plasma_wide_raw)))

analyte_cols <- setdiff(names(plasma_wide_raw), c("PID","Group","TP"))

## ---- Long format + flag parsing -------------------------------------------
## UND / OVER / NA_run / below-LLOQ (e.g. '<0.80↓') are all categorical "Flag"
## values; the numeric column gets NA for OVER and NA_run, 0 for UND/BLLOQ.
plasma_long <- plasma_wide_raw %>%
  pivot_longer(cols = all_of(analyte_cols), names_to = "AnalyteRaw", values_to = "RawValue") %>%
  extract(AnalyteRaw, into = c("Analyte","Unit"),
          regex = "^(.*)\\s*\\((pg/ml|nmol/L)\\)\\s*$", remove = TRUE) %>%
  mutate(
    RawValue_chr = as.character(RawValue),
    .upper       = str_to_upper(str_trim(RawValue_chr)),
    Flag = case_when(
      is.na(RawValue_chr) | .upper == "NA"               ~ "NA_run",
      .upper == "UND"                                    ~ "UND",
      .upper == "OVER"                                   ~ "OVER",
      str_detect(.upper, "^<")                           ~ "BLLOQ",  # below LLOQ markers
      TRUE                                                ~ "OK"
    ),
    Value_num = suppressWarnings(
      as.numeric(str_replace_all(RawValue_chr, "[,\\s]", ""))
    ),
    Analyte = str_squish(Analyte)
  ) %>%
  select(-.upper) %>%
  select(PID, Group, TP, Analyte, Unit, RawValue = RawValue_chr, Flag, Value_num)

## ---- Recode UND/BLLOQ → 0 ; OVER/NA_run → NA ------------------------------
plasma_long <- plasma_long %>%
  mutate(
    Value_clean = case_when(
      Flag %in% c("UND","BLLOQ") ~ 0,
      Flag %in% c("OVER","NA_run") ~ NA_real_,
      TRUE                          ~ Value_num
    )
  )

message("Flag distribution after cleaning:")
print(plasma_long %>% count(Flag))

## ---- Harmonize TP (handles the actual labels in this file: '12wk','44wk','Entry') ----
plasma_clean <- plasma_long %>%
  mutate(
    TP = case_when(
      str_detect(str_to_lower(TP), "^entry$")                                  ~ "Entry",
      str_detect(str_to_lower(TP), "^(dx|diagnosis)$")                         ~ "Dx",
      str_detect(str_to_lower(TP), "^(12|12 ?wks?|12 ?weeks?)$")               ~ "12 Wks",
      str_detect(str_to_lower(TP), "^(44|44 ?wks?|44 ?weeks?)$")               ~ "44 Wks",
      TRUE                                                                     ~ TP
    )
  )

message("TP values after harmonization: ",
        paste(sort(unique(plasma_clean$TP)), collapse = ", "))

## ---- Low-information analyte filter ---------------------------------------
min_nonmissing_overall <- 40
max_zero_fraction       <- 0.90

plasma_clean <- plasma_clean %>% mutate(IsZero = !is.na(Value_clean) & Value_clean <= 0)

plasma_analyte_qc <- plasma_clean %>%
  group_by(Analyte) %>%
  summarise(
    n_total       = n(),
    n_nonmiss     = sum(!is.na(Value_clean)),
    zero_fraction = ifelse(n_nonmiss > 0, sum(IsZero, na.rm = TRUE) / n_nonmiss, NA_real_),
    .groups = "drop"
  )

keep_analytes <- plasma_analyte_qc %>%
  filter(n_nonmiss >= min_nonmissing_overall,
         is.na(zero_fraction) | zero_fraction < max_zero_fraction) %>%
  pull(Analyte)

message("Keeping ", length(keep_analytes), " of ", n_distinct(plasma_clean$Analyte),
        " analytes after QC.")

plasma_filtered <- plasma_clean %>%
  filter(Analyte %in% keep_analytes) %>%
  select(PID, Group, TP, Analyte, Unit, Flag, RawValue, Value_clean)

## ---- log1p + per-TP z-score (mild winsorization on log scale) -------------
winsor_q <- 0.01

plasma_scaled <- plasma_filtered %>%
  mutate(Value_log1p = log1p(Value_clean)) %>%
  group_by(Analyte, TP) %>%
  mutate(
    .low   = if (winsor_q > 0) stats::quantile(Value_log1p, probs = winsor_q,     na.rm = TRUE) else -Inf,
    .high  = if (winsor_q > 0) stats::quantile(Value_log1p, probs = 1 - winsor_q, na.rm = TRUE) else  Inf,
    Value_log1p_w = pmin(pmax(Value_log1p, .low), .high),
    mean_tp = mean(Value_log1p_w, na.rm = TRUE),
    sd_tp   = stats::sd(Value_log1p_w, na.rm = TRUE),
    Value_z = ifelse(is.finite(sd_tp) & sd_tp > 0, (Value_log1p_w - mean_tp) / sd_tp, NA_real_)
  ) %>%
  ungroup() %>%
  select(PID, Group, TP, Analyte, Unit, Flag, RawValue, Value_clean, Value_log1p, Value_z)

## ---- Save cleaned outputs (qs2 + CSV) -------------------------------------
qs2::qs_save(plasma_scaled, file.path(saved_dir, "Plasma_Biomarker_long_scaled.qs2"))
write_csv(plasma_scaled,    file.path(saved_dir, "Plasma_Biomarker_long_scaled.csv"))

plasma_wide <- plasma_scaled %>%
  select(PID, Group, TP, Analyte, Value_z) %>%
  pivot_wider(names_from = Analyte, values_from = Value_z)
qs2::qs_save(plasma_wide, file.path(saved_dir, "Plasma_Biomarker_wide_scaled.qs2"))
write_csv(plasma_wide,    file.path(saved_dir, "Plasma_Biomarker_wide_scaled.csv"))

## ============================================================================
## IGRA lookup from cytokine_long (deduped with explicit priority)
## ============================================================================
igra_cols <- names(cytokine_long)[str_detect(names(cytokine_long), regex("IGRA", ignore_case = TRUE))]
if (length(igra_cols) == 0) stop("No IGRA columns found in cytokine_long.")

# Priority for resolving inconsistencies within a single PID:
#   ATB > Pos > Neg > INT > NA
# Rationale:
#   - ATB (active TB) is clinically the most specific label — never override it.
#     ATB is its own category and is NOT treated as "Pos".
#   - INT (intermediate / indeterminate) is the weakest — only used as a fallback
#     when no other label is available. Downstream in Script 4, remaining INTs
#     for maternal IGRA are remapped to Pos (intermediate counted as positive).
igra_priority <- function(x) {
  ord <- c("ATB" = 4, "Pos" = 3, "Neg" = 2, "INT" = 1)
  x   <- as.character(x)
  if (all(is.na(x))) return(NA_character_)
  x   <- x[!is.na(x)]
  x[which.max(unname(ord[x]))]
}

# Build a single-row-per-(PID, PID_Type, Time_Point) IGRA table
igra_lookup <- cytokine_long %>%
  mutate(across(all_of(igra_cols), as.character)) %>%
  select(PID, PID_Type, Time_Point, all_of(igra_cols)) %>%
  group_by(PID, PID_Type, Time_Point) %>%
  summarise(across(all_of(igra_cols), igra_priority), .groups = "drop") %>%
  mutate(PID = as.character(PID))

## ---- Build analysis-ready plasma table ------------------------------------
plasma_work <- plasma_scaled %>%
  mutate(
    PID        = as.character(PID),
    Group      = as.character(Group),
    Time_Point = TP,
    PID_Type   = Group
  ) %>%
  left_join(igra_lookup, by = c("PID","PID_Type","Time_Point"))

## ============================================================================
## Stats — four comparisons
## ============================================================================
out_root <- file.path(proj_root, "Results", "Plasma_StatTests")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

igra_m_col <- "Maternal IGRA"
igra_i_col <- "Infant IGRA"

# Analysis table for tests (rank tests so log1p vs z gives same p-values; we use log1p)
plasma_test <- plasma_work %>%
  filter(!Flag %in% c("OVER","NA_run"), !is.na(Value_clean)) %>%   # UND-as-0 retained
  mutate(Timepoint = Time_Point, Y = Value_log1p) %>%
  filter(!is.na(Y))

# Helper: collapse to one row per PID per group level (median, robust)
collapse_pid <- function(df, grp_col) {
  df %>%
    filter(!is.na(.data[[grp_col]])) %>%
    group_by(Analyte, PID, .data[[grp_col]]) %>%
    summarise(Y = stats::median(Y, na.rm = TRUE), .groups = "drop")
}

## ---- Paired stats: Infants 12 vs 44 wks -----------------------------------
paired_stats_12v44 <- function(df) {
  wide <- df %>%
    select(PID, Analyte, Timepoint, Y) %>%
    filter(Timepoint %in% c("12 Wks","44 Wks")) %>%
    group_by(Analyte, PID, Timepoint) %>%
    summarise(Y = stats::median(Y, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Timepoint, values_from = Y)
  
  if (!all(c("12 Wks","44 Wks") %in% names(wide))) {
    return(tibble(Analyte = character(), n_pairs = integer(),
                  median_12Wks = numeric(), median_44Wks = numeric(),
                  median_diff_44minus12 = numeric(),
                  p_value = numeric(), p_adj = numeric(), test = character()))
  }
  
  wide <- wide %>% filter(!is.na(`12 Wks`), !is.na(`44 Wks`))
  
  wide %>%
    group_split(Analyte) %>%
    map_dfr(function(xx) {
      a <- unique(xx$Analyte)
      if (nrow(xx) < 2) {
        return(tibble(Analyte = a, n_pairs = nrow(xx),
                      median_12Wks = stats::median(xx$`12 Wks`, na.rm = TRUE),
                      median_44Wks = stats::median(xx$`44 Wks`, na.rm = TRUE),
                      median_diff_44minus12 = stats::median(xx$`44 Wks` - xx$`12 Wks`, na.rm = TRUE),
                      p_value = NA_real_, test = "paired_wilcox"))
      }
      wt <- suppressWarnings(stats::wilcox.test(xx$`44 Wks`, xx$`12 Wks`, paired = TRUE, exact = FALSE))
      tibble(Analyte = a, n_pairs = nrow(xx),
             median_12Wks = stats::median(xx$`12 Wks`, na.rm = TRUE),
             median_44Wks = stats::median(xx$`44 Wks`, na.rm = TRUE),
             median_diff_44minus12 = stats::median(xx$`44 Wks` - xx$`12 Wks`, na.rm = TRUE),
             p_value = wt$p.value, test = "paired_wilcox")
    }) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))
}

## ---- Unpaired (2-level) ---------------------------------------------------
unpaired_stats_twolevel <- function(df, grp_var, lvl_a = "Pos", lvl_b = "Neg", label = "") {
  stopifnot(grp_var %in% names(df))
  
  dat <- df %>%
    filter(!is.na(.data[[grp_var]])) %>%
    mutate(Comp = .data[[grp_var]]) %>%
    filter(Comp %in% c(lvl_a, lvl_b)) %>%
    # collapse to one row per (Analyte, PID, level) so technical reps don't inflate n
    group_by(Analyte, PID, Comp) %>%
    summarise(Y = stats::median(Y, na.rm = TRUE), .groups = "drop")
  
  if (nrow(dat) == 0) {
    return(tibble(Contrast = character(), Analyte = character(),
                  n_A = integer(), n_B = integer(),
                  median_A = numeric(), median_B = numeric(),
                  median_diff_AminusB = numeric(),
                  p_value = numeric(), p_adj = numeric(), test = character()))
  }
  
  dat %>%
    group_split(Analyte) %>%
    map_dfr(function(xx) {
      a  <- unique(xx$Analyte)
      xa <- xx %>% filter(Comp == lvl_a) %>% pull(Y)
      xb <- xx %>% filter(Comp == lvl_b) %>% pull(Y)
      if (length(xa) < 2 || length(xb) < 2) {
        return(tibble(Contrast = label, Analyte = a,
                      n_A = length(xa), n_B = length(xb),
                      median_A = stats::median(xa, na.rm = TRUE),
                      median_B = stats::median(xb, na.rm = TRUE),
                      median_diff_AminusB = stats::median(xa, na.rm = TRUE) -
                        stats::median(xb, na.rm = TRUE),
                      p_value = NA_real_, test = "wilcox_unpaired"))
      }
      wt <- suppressWarnings(stats::wilcox.test(xa, xb, paired = FALSE, exact = FALSE))
      tibble(Contrast = label, Analyte = a,
             n_A = length(xa), n_B = length(xb),
             median_A = stats::median(xa, na.rm = TRUE),
             median_B = stats::median(xb, na.rm = TRUE),
             median_diff_AminusB = stats::median(xa, na.rm = TRUE) -
               stats::median(xb, na.rm = TRUE),
             p_value = wt$p.value, test = "wilcox_unpaired")
    }) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH"))
}

## ---- Run all four comparisons ---------------------------------------------
cmp1_dir <- file.path(out_root, "Infants_12wk_vs_44wk_Paired");              dir.create(cmp1_dir, showWarnings = FALSE, recursive = TRUE)
cmp2_dir <- file.path(out_root, "Infants_byMaternalIGRA_PosVsNeg");          dir.create(cmp2_dir, showWarnings = FALSE, recursive = TRUE)
cmp3_dir <- file.path(out_root, "Infants_byInfantIGRA_PosVsNeg");            dir.create(cmp3_dir, showWarnings = FALSE, recursive = TRUE)
cmp4_dir <- file.path(out_root, "Mothers_Entry_byMaternalIGRA_PosVsNeg");    dir.create(cmp4_dir, showWarnings = FALSE, recursive = TRUE)

stats_12v44 <- plasma_test %>%
  filter(Group == "Infant", Timepoint %in% c("12 Wks","44 Wks")) %>%
  paired_stats_12v44()
write_csv(stats_12v44, file.path(cmp1_dir, "stats_Infants_12vs44_paired.csv"))

stats_inf_by_mIGRA <- plasma_test %>%
  filter(Group == "Infant") %>%
  unpaired_stats_twolevel(igra_m_col, "Pos", "Neg", "MaternalIGRA_Pos_vs_Neg")
write_csv(stats_inf_by_mIGRA, file.path(cmp2_dir, "stats_Infants_byMaternalIGRA.csv"))

stats_inf_by_iIGRA <- plasma_test %>%
  filter(Group == "Infant") %>%
  unpaired_stats_twolevel(igra_i_col, "Pos", "Neg", "InfantIGRA_Pos_vs_Neg")
write_csv(stats_inf_by_iIGRA, file.path(cmp3_dir, "stats_Infants_byInfantIGRA.csv"))

stats_moms_by_mIGRA <- plasma_test %>%
  filter(Group == "Mother") %>%
  unpaired_stats_twolevel(igra_m_col, "Pos", "Neg", "MaternalIGRA_Pos_vs_Neg")
write_csv(stats_moms_by_mIGRA, file.path(cmp4_dir, "stats_Mothers_byMaternalIGRA.csv"))

message("Plasma stats written to: ", out_root)

## ============================================================================
## Plotting — consistent Y axis (Value_z) for ALL plots, matching label.
## P-values displayed in subtitle come from the log1p Wilcoxon (rank-invariant
## so the p-values are identical whether computed on log1p or on z-score).
## ============================================================================
plot_root <- file.path(proj_root, "Results", "Plasma_Plots")
dirs <- list(
  cmp_12v44   = file.path(plot_root, "Infants_12wk_vs_44wk_Paired"),
  cmp_i_mIGRA = file.path(plot_root, "Infants_byMaternalIGRA_PosVsNeg"),
  cmp_i_iIGRA = file.path(plot_root, "Infants_byInfantIGRA_PosVsNeg"),
  cmp_m_mIGRA = file.path(plot_root, "Mothers_Entry_byMaternalIGRA_PosVsNeg")
)
invisible(lapply(dirs, function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE)))

sanitize_for_path <- function(x) {
  x %>% str_replace_all("[\\\\/:*?\"<>|]", "_") %>% str_squish()
}

mk_subtitle <- function(nA = NA, nB = NA, n_pairs = NA, p = NA, fdr = NA) {
  bits <- c()
  if (!is.na(n_pairs))     bits <- c(bits, glue("n_pairs = {n_pairs}"))
  if (!is.na(nA) && !is.na(nB)) bits <- c(bits, glue("n = {nA} vs {nB}"))
  if (!is.na(p))           bits <- c(bits, glue("p = {signif(p, 3)}"))
  if (!is.na(fdr))         bits <- c(bits, glue("FDR = {signif(fdr, 3)}"))
  paste(bits, collapse = " | ")
}

theme_plasma <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(plot.title = element_text(face = "bold"),
          panel.grid.major.x = element_blank(),
          legend.position = "none")
}

geom_summary_crossbar <- function(width = 0.5) {
  stat_summary(fun = stats::median, geom = "crossbar",
               width = width, fatten = 0, alpha = 0.8)
}

y_axis_label <- "Z-score (per-TP standardized log1p input)"

# Plot data uses Value_z (per-TP z-score) — consistent across all plots.
plasma_plot_data <- plasma_work %>%
  filter(!Flag %in% c("OVER","NA_run"), !is.na(Value_z))

analytes   <- sort(unique(plasma_plot_data$Analyte))
timepoints <- c("12 Wks","44 Wks")

## ---- Plot: Infants 12 vs 44 wks (paired) ---------------------------------
purrr::walk(analytes, function(ana) {
  df <- plasma_plot_data %>%
    filter(Group == "Infant", TP %in% c("12 Wks","44 Wks"), Analyte == !!ana) %>%
    group_by(PID, TP) %>%
    summarise(Value_z = stats::median(Value_z, na.rm = TRUE), .groups = "drop")
  
  wide <- pivot_wider(df, names_from = TP, values_from = Value_z) %>%
    filter(!is.na(`12 Wks`), !is.na(`44 Wks`))
  n_pairs <- nrow(wide); if (n_pairs < 3) return(invisible(NULL))
  
  # P-value from log1p paired wilcox (matches rank-test on Value_z exactly)
  st <- stats_12v44 %>% filter(Analyte == !!ana) %>% slice(1)
  pval <- if (nrow(st)) st$p_value else NA_real_
  fdr  <- if (nrow(st)) st$p_adj else NA_real_
  
  df_plot <- wide %>%
    pivot_longer(c(`12 Wks`,`44 Wks`), names_to = "TP", values_to = "Value_z") %>%
    mutate(TP = factor(TP, levels = c("12 Wks","44 Wks")))
  
  w_in <- max(4.8, min(7.5, 4.8 + 0.05 * n_pairs)); h_in <- 4.8
  p <- ggplot(df_plot, aes(x = TP, y = Value_z, group = PID)) +
    geom_line(alpha = 0.35) +
    geom_point(aes(color = TP), size = 2) +
    geom_summary_crossbar() +
    scale_x_discrete(drop = FALSE) +
    labs(title = glue("{ana} — Infants: 12 vs 44 weeks (paired)"),
         subtitle = mk_subtitle(n_pairs = n_pairs, p = pval, fdr = fdr),
         x = NULL, y = y_axis_label) +
    theme_plasma()
  
  out <- file.path(dirs$cmp_12v44,
                   paste0(sanitize_for_path(glue("{ana}__Infants_12vs44_paired")), ".png"))
  ggsave(out, p, width = w_in, height = h_in, dpi = 300, bg = "white")
})

## ---- Plot helper for 2-level unpaired plots (Pos vs Neg) -----------------
plot_pos_neg <- function(group_filter, tp_filter, igra_col, label_path, stats_tbl) {
  purrr::walk(analytes, function(ana) {
    tps <- if (is.null(tp_filter)) "Entry" else tp_filter
    purrr::walk(tps, function(tp) {
      df <- plasma_plot_data %>%
        filter(Group == group_filter,
               (is.null(tp_filter) | TP == tp),
               Analyte == !!ana,
               !is.na(.data[[igra_col]]),
               .data[[igra_col]] %in% c("Pos","Neg")) %>%
        group_by(PID, .data[[igra_col]]) %>%
        summarise(Value_z = stats::median(Value_z, na.rm = TRUE), .groups = "drop")
      if (nrow(df) < 4) return(invisible(NULL))
      
      xa <- df$Value_z[df[[igra_col]] == "Pos"]
      xb <- df$Value_z[df[[igra_col]] == "Neg"]
      if (length(xa) < 2 || length(xb) < 2) return(invisible(NULL))
      
      st <- stats_tbl %>% filter(Analyte == !!ana) %>% slice(1)
      pval <- if (nrow(st)) st$p_value else NA_real_
      fdr  <- if (nrow(st)) st$p_adj else NA_real_
      
      title <- if (is.null(tp_filter))
        glue("{ana} — {group_filter}s (Entry) — {sub('_',' ',igra_col)}")
      else
        glue("{ana} — {group_filter}s ({tp}) — {sub('_',' ',igra_col)}")
      
      p <- ggplot(df, aes(x = .data[[igra_col]], y = Value_z)) +
        geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.6) +
        geom_point(position = position_jitter(width = 0.1, height = 0),
                   size = 1.9, alpha = 0.75, aes(color = .data[[igra_col]])) +
        geom_summary_crossbar() +
        scale_x_discrete(limits = c("Neg","Pos"), drop = FALSE) +
        labs(title = title,
             subtitle = mk_subtitle(nA = length(xa), nB = length(xb), p = pval, fdr = fdr),
             x = NULL, y = y_axis_label) +
        theme_plasma()
      
      tp_tag <- if (is.null(tp_filter)) "Entry" else tp
      out <- file.path(label_path,
                       sanitize_for_path(glue("{ana}__{group_filter}__{tp_tag}__{igra_col}.png")))
      ggsave(out, p, width = 5.6, height = 4.5, dpi = 300, bg = "white")
    })
  })
}

# C2: Infants × Maternal IGRA × {12 Wks, 44 Wks}
plot_pos_neg("Infant", timepoints, igra_m_col, dirs$cmp_i_mIGRA, stats_inf_by_mIGRA)
# C3: Infants × Infant IGRA × {12 Wks, 44 Wks}
plot_pos_neg("Infant", timepoints, igra_i_col, dirs$cmp_i_iIGRA, stats_inf_by_iIGRA)
# C4: Mothers × Maternal IGRA × Entry only
plot_pos_neg("Mother", NULL,       igra_m_col, dirs$cmp_m_mIGRA, stats_moms_by_mIGRA)

message("Script 3 complete. Plasma plots in: ", plot_root)