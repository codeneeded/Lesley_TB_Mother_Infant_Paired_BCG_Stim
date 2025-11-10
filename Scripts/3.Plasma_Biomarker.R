# Core tidyverse + plotting
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(purrr)

# For safer stats & labeling later
library(broom)
library(glue)

# ---- Paths (edit if needed) ----
# Linux/mac path:
#saved_dir <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/saved_R_dat"

# Windows example (leave commented if using Linux):
 saved_dir <- "C:/Users/ammas/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/saved_R_dat"

### load data
 
### cytokine Data
in_rds <- file.path(saved_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.rds")
in_csv <- file.path(saved_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.csv")

# ---- Load (prefer RDS; fallback to CSV) ----
if (file.exists(in_rds)) {
  cytokine_long <- readRDS(in_rds)
} else if (file.exists(in_csv)) {
  cytokine_long <- readr::read_csv(in_csv, show_col_types = FALSE)
} else {
  stop("Could not find preprocessed input at:\n  ", in_rds, "\n  ", in_csv)
}

### Plasma Biomarker Data
# ---- Path to the Excel with both Flow & Plasma ----
# Linux/mac:
#plasma_xlsx <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/Raw_Data/TB Flow and Plasma Biomarker data.xlsx"
# Windows example:
plasma_xlsx <- "C:/Users/ammas/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/Raw_Data/TB Flow and Plasma Biomarker data.xlsx"

# ---- Read the "Plasma Biomarker" sheet ----
plasma_wide <- readxl::read_excel(
  path = plasma_xlsx,
  sheet = "Plasma Biomarker",
  .name_repair = "minimal"
)


# Identify analyte columns (everything after PID, Group, TP)
stopifnot(all(c("PID","Group","TP") %in% names(plasma_wide)))
analyte_cols <- setdiff(names(plasma_wide), c("PID","Group","TP"))

# Long format + basic flag parsing (UND/OVER/NA) and numeric extraction
plasma_long <- plasma_wide %>%
  tidyr::pivot_longer(
    cols = all_of(analyte_cols),
    names_to  = "AnalyteRaw",
    values_to = "RawValue"
  ) %>%
  # Extract analyte and unit from names like "IL-6 (pg/ml)" or "Neopterin (nmol/L)"
  tidyr::extract(AnalyteRaw, into = c("Analyte","Unit"),
                 regex = "^(.*)\\s*\\((pg/ml|nmol/L)\\)\\s*$", remove = TRUE) %>%
  dplyr::mutate(
    RawValue_chr = as.character(RawValue),
    Flag = dplyr::case_when(
      is.na(RawValue_chr)                              ~ "NA_run",
      stringr::str_to_upper(stringr::str_trim(RawValue_chr)) == "UND"  ~ "UND",
      stringr::str_to_upper(stringr::str_trim(RawValue_chr)) == "OVER" ~ "OVER",
      stringr::str_to_upper(stringr::str_trim(RawValue_chr)) == "NA"   ~ "NA_run",
      TRUE ~ "OK"
    ),
    # Coerce numeric (strip commas/spaces); keep NA if flagged
    Value_num = suppressWarnings(
      as.numeric(stringr::str_replace_all(RawValue_chr, "[,\\s]", ""))
    ),
    Value = dplyr::if_else(Flag %in% c("UND","OVER","NA_run"), NA_real_, Value_num)
  ) %>%
  dplyr::select(PID, Group, TP, Analyte, Unit, RawValue = RawValue_chr, Flag, Value)



levels(as.factor(plasma_long$Flag))
str(plasma_long)

############# 
#---- Step 4.0 (fix): trust existing Flag/Value; create Value_num & trim Analyte ----
##############
plasma_long_v1 <- plasma_long %>%
  mutate(
    # clean trailing spaces on analyte names
    Analyte   = stringr::str_squish(Analyte),
    
    # working numeric column:
    # UND -> 0; OVER/NA_run -> NA; OK -> keep Value
    Value_num = dplyr::case_when(
      Flag == "UND"    ~ 0,
      Flag == "OVER"   ~ NA_real_,
      Flag == "NA_run" ~ NA_real_,
      TRUE             ~ as.numeric(Value)
    )
  )

# quick sanity check — you should see counts, not an empty tibble
plasma_long_v1 %>%
  count(Flag)

plasma_long_v1 %>%
  summarise(
    n_rows = n(),
    n_miss = sum(is.na(Value_num)),
    n_zero = sum(Value_num == 0, na.rm = TRUE)
  )
str(plasma_long_v1)

## ---- Step 4.1: prune & standardize ----
plasma_clean <- plasma_long_v1 %>%
  # drop OVER entirely
  dplyr::filter(Flag != "OVER") %>%
  # harmonize TP labels to: Entry, Dx, 12 Wks, 44 Wks
  dplyr::mutate(
    TP = dplyr::case_when(
      stringr::str_to_lower(TP) %in% c("entry")                              ~ "Entry",
      stringr::str_to_lower(TP) %in% c("dx", "diagnosis")                    ~ "Dx",
      stringr::str_to_lower(TP) %in% c("12", "12wks", "12 weeks", "12 wks")  ~ "12 Wks",
      stringr::str_to_lower(TP) %in% c("44", "44wks", "44 weeks", "44 wks")  ~ "44 Wks",
      TRUE ~ TP
    ),
    # working numeric to use everywhere downstream
    Value_clean = Value_num
  ) %>%
  # keep tidy column order
  dplyr::select(PID, Group, TP, Analyte, Unit, Flag, RawValue, Value, Value_clean)

# quick sanity check
plasma_clean %>% count(Flag)
levels(as.factor(plasma_clean$TP))

## ---- Step 4.2: harmonize TP + low-information screening ----

# thresholds you can tweak
min_nonmissing_overall <- 40   # keep analytes with at least this many non-NA values overall
max_zero_fraction       <- 0.90 # drop analytes with ≥90% zeros (after UND→0)

plasma_clean2 <- plasma_clean %>%
  dplyr::mutate(
    TP = dplyr::case_when(
      stringr::str_detect(stringr::str_to_lower(TP), "^entry$")                       ~ "Entry",
      stringr::str_detect(stringr::str_to_lower(TP), "^(dx|diagnosis)$")              ~ "Dx",
      stringr::str_detect(stringr::str_to_lower(TP), "^(12|12 ?wks?|12 ?weeks?)$")    ~ "12 Wks",
      stringr::str_detect(stringr::str_to_lower(TP), "^(44|44 ?wks?|44 ?weeks?)$")    ~ "44 Wks",
      TRUE ~ TP
    ),
    IsZero = !is.na(Value_clean) & Value_clean <= 0
  )

# per-analyte QC stats
plasma_analyte_qc <- plasma_clean2 %>%
  dplyr::group_by(Analyte) %>%
  dplyr::summarise(
    n_total      = dplyr::n(),
    n_nonmiss    = sum(!is.na(Value_clean)),
    zero_fraction= ifelse(n_nonmiss > 0, sum(IsZero, na.rm = TRUE) / n_nonmiss, NA_real_),
    .groups = "drop"
  )

# filter analytes that pass basic information criteria
keep_analytes <- plasma_analyte_qc %>%
  dplyr::filter(n_nonmiss >= min_nonmissing_overall,
                is.na(zero_fraction) | zero_fraction < max_zero_fraction) %>%
  dplyr::pull(Analyte)

plasma_filtered <- plasma_clean2 %>%
  dplyr::filter(Analyte %in% keep_analytes) %>%
  dplyr::select(PID, Group, TP, Analyte, Unit, Flag, RawValue, Value_clean)

# quick sanity checks
plasma_analyte_qc %>% dplyr::arrange(zero_fraction) %>% head(10)
plasma_analyte_qc %>% dplyr::arrange(dplyr::desc(zero_fraction)) %>% head(10)
levels(as.factor(plasma_filtered$TP))
length(unique(plasma_filtered$Analyte))

## ---- Step 5: transform + scale within timepoint (TP) ----
## log1p to stabilize heavy tails, then z-score within TP for each analyte
## optional: mild winsorization to reduce extreme outliers

winsor_q <- 0.01  # set to 0 to disable winsorization

plasma_scaled <- plasma_filtered %>%
  dplyr::mutate(
    # transform
    Value_log1p = log1p(Value_clean)
  ) %>%
  dplyr::group_by(Analyte, TP) %>%
  dplyr::mutate(
    # optional winsorization on the log scale
    .low  = if (winsor_q > 0) quantile(Value_log1p, probs = winsor_q, na.rm = TRUE) else -Inf,
    .high = if (winsor_q > 0) quantile(Value_log1p, probs = 1 - winsor_q, na.rm = TRUE) else  Inf,
    Value_log1p_w = pmin(pmax(Value_log1p, .low), .high),
    
    # z-score within TP (analyte-wise)
    mean_tp = mean(Value_log1p_w, na.rm = TRUE),
    sd_tp   = stats::sd(Value_log1p_w,  na.rm = TRUE),
    Value_z = ifelse(is.finite(sd_tp) & sd_tp > 0, (Value_log1p_w - mean_tp) / sd_tp, NA_real_)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(PID, Group, TP, Analyte, Unit, Flag, RawValue,
                Value_clean, Value_log1p, Value_z)

# quick checks
plasma_scaled %>% dplyr::count(TP)
plasma_scaled %>% dplyr::group_by(Analyte, TP) %>% dplyr::summarise(mu = mean(Value_z, na.rm=TRUE), sd = sd(Value_z, na.rm=TRUE), .groups="drop") %>% head()

## ---- Step 6: save cleaned plasma data (Windows paths) ----

out_dir <- "C:/Users/ammas/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/saved_R_dat"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# long format (clean + scaled)
saveRDS(plasma_scaled, file.path(out_dir, "Plasma_Biomarker_long_scaled.rds"))
readr::write_csv(plasma_scaled, file.path(out_dir, "Plasma_Biomarker_long_scaled.csv"))

# also make a wide matrix (for correlation / integration)
plasma_wide <- plasma_scaled %>%
  dplyr::select(PID, Group, TP, Analyte, Value_z) %>%
  tidyr::pivot_wider(names_from = Analyte, values_from = Value_z)

saveRDS(plasma_wide, file.path(out_dir, "Plasma_Biomarker_wide_scaled.rds"))
readr::write_csv(plasma_wide, file.path(out_dir, "Plasma_Biomarker_wide_scaled.csv"))

# sanity check
dim(plasma_wide)
head(colnames(plasma_wide))

####################################################
## ---- Step A: output dirs + working plasma table ----


# 1) Where to save things (Windows paths)
out_root <- "C:/Users/ammas/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/Results_Plasma"
subdirs  <- c(
  "Infants_12vs44",
  "Infants_ByMaternalIGRA",
  "Infants_ByInfantIGRA",
  "Mothers_ByMaternalIGRA"
)
for (sd in subdirs) dir.create(file.path(out_root, sd, "csv"), recursive = TRUE, showWarnings = FALSE)
for (sd in subdirs) dir.create(file.path(out_root, sd, "plots"), recursive = TRUE, showWarnings = FALSE)

# 2) Small helper to make safe filenames
sanitize_for_path <- function(x) {
  x %>%
    stringr::str_replace_all("[/:*?\"<>|]+", "_") %>%
    stringr::str_squish()
}

# 3) Build an IGRA lookup from cytokine_long (if needed)
#    (We’ll robustly find any columns containing 'IGRA')
igra_cols <- names(cytokine_long)[stringr::str_detect(names(cytokine_long), regex("IGRA", ignore_case = TRUE))]
if (length(igra_cols) == 0) stop("No IGRA columns found in cytokine_long. Please confirm IGRA column names.")

igra_lookup <- cytokine_long %>%
  dplyr::select(PID, PID_Type, Time_Point, dplyr::all_of(igra_cols)) %>%
  dplyr::distinct()

# 4) Choose your plasma base table (prefer your most processed one if it exists)
plasma_base <- plasma_scaled

# 5) Harmonize keys and attach IGRA if missing
plasma_work <- plasma_base %>%
  dplyr::mutate(
    PID       = as.character(PID),
    Group     = as.character(Group),
    TP        = as.character(TP),
    Time_Point = TP,         # align name to cytokine_long
    PID_Type   = Group       # 'Mother' / 'Infant'
  )

missing_igra <- setdiff(igra_cols, names(plasma_work))
if (length(missing_igra) > 0) {
  plasma_work <- plasma_work %>%
    dplyr::left_join(igra_lookup, by = c("PID","PID_Type","Time_Point"))
}

# 6) Quick sanity
message("plasma_work rows: ", nrow(plasma_work))
message("Columns now include IGRA fields: ", paste(intersect(names(plasma_work), igra_cols), collapse = ", "))
str(plasma_work)

## ---------------------------
## Step B — Plasma: stats for 4 comparisons (save CSVs)
## ---------------------------

## 0) Config / paths (Windows)
out_root <- "C:/Users/ammas/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/Plasma_Results"
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

## 1) Column hygiene
tp_col <- if ("Time_Point" %in% names(plasma_work)) "Time_Point" else if ("TP" %in% names(plasma_work)) "TP" else stop("No timepoint column found (Time_Point/TP)")
igra_m_col <- "Maternal IGRA"
igra_i_col <- "Infant IGRA"

## 2) Analysis-ready table
plasma_test <- plasma_work %>%
  filter(Flag == "OK", !is.na(Value_clean)) %>%     # keep only valid numeric runs
  mutate(
    Timepoint = .data[[tp_col]],
    Group     = .data[["Group"]],
    Analyte   = .data[["Analyte"]],
    Y         = .data[["Value_log1p"]]             # use log1p for tests
  ) %>%
  filter(!is.na(Y))

## ---------- helper: paired stats (Infants 12 vs 44 wks) WITH MEDIANS ----------
paired_stats_12v44 <- function(df) {
  # df is filtered to Group==Infant & Timepoint in c("12 Wks","44 Wks")
  w <- df %>%
    dplyr::select(PID, Analyte, Timepoint, Y) %>%
    dplyr::filter(Timepoint %in% c("12 Wks","44 Wks")) %>%
    dplyr::group_by(Analyte, PID) %>%
    dplyr::filter(dplyr::n_distinct(Timepoint) == 2) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = Timepoint, values_from = Y)
  
  if (!all(c("12 Wks","44 Wks") %in% names(w))) {
    return(tibble::tibble(
      Analyte = character(), n_pairs = integer(),
      median_12Wks = numeric(), median_44Wks = numeric(),
      median_diff_44minus12 = numeric(),
      p_value = numeric(), p_adj = numeric(), test = character()
    ))
  }
  
  w <- w %>% dplyr::filter(!is.na(`12 Wks`), !is.na(`44 Wks`))
  
  by_analyte <- w %>%
    dplyr::group_split(Analyte, .keep = TRUE) %>%
    purrr::map_dfr(function(xx) {
      a <- unique(xx$Analyte)
      if (nrow(xx) < 2) {
        return(tibble::tibble(
          Analyte = a, n_pairs = nrow(xx),
          median_12Wks = stats::median(xx$`12 Wks`, na.rm = TRUE),
          median_44Wks = stats::median(xx$`44 Wks`, na.rm = TRUE),
          median_diff_44minus12 = stats::median(xx$`44 Wks` - xx$`12 Wks`, na.rm = TRUE),
          p_value = NA_real_, p_adj = NA_real_, test = "paired_wilcox"
        ))
      }
      wt <- suppressWarnings(stats::wilcox.test(xx$`44 Wks`, xx$`12 Wks`, paired = TRUE, exact = FALSE))
      tibble::tibble(
        Analyte = a,
        n_pairs = nrow(xx),
        median_12Wks = stats::median(xx$`12 Wks`, na.rm = TRUE),
        median_44Wks = stats::median(xx$`44 Wks`, na.rm = TRUE),
        median_diff_44minus12 = stats::median(xx$`44 Wks` - xx$`12 Wks`, na.rm = TRUE),
        p_value = wt$p.value,
        test = "paired_wilcox"
      )
    }) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"))
  
  by_analyte
}

## ---------- helper: unpaired stats (2-level) WITH MEDIANS ----------
unpaired_stats_twolevel <- function(df, grp_var, lvl_a = "Pos", lvl_b = "Neg", label = "") {
  stopifnot(grp_var %in% names(df))
  
  dat <- df %>%
    dplyr::filter(!is.na(.data[[grp_var]])) %>%
    dplyr::mutate(Comp = .data[[grp_var]]) %>%
    dplyr::filter(Comp %in% c(lvl_a, lvl_b))
  
  if (nrow(dat) == 0) {
    return(tibble::tibble(
      Contrast = character(), Analyte = character(),
      n_A = integer(), n_B = integer(),
      median_A = numeric(), median_B = numeric(),
      median_diff_AminusB = numeric(),
      p_value = numeric(), p_adj = numeric(), test = character()
    ))
  }
  
  by_analyte <- dat %>%
    dplyr::group_split(Analyte, .keep = TRUE) %>%
    purrr::map_dfr(function(xx) {
      a <- unique(xx$Analyte)
      xa <- xx %>% dplyr::filter(Comp == lvl_a) %>% dplyr::pull(Y)
      xb <- xx %>% dplyr::filter(Comp == lvl_b) %>% dplyr::pull(Y)
      
      if (length(xa) < 2 || length(xb) < 2) {
        return(tibble::tibble(
          Contrast = label, Analyte = a,
          n_A = length(xa), n_B = length(xb),
          median_A = stats::median(xa, na.rm = TRUE),
          median_B = stats::median(xb, na.rm = TRUE),
          median_diff_AminusB = stats::median(xa, na.rm = TRUE) - stats::median(xb, na.rm = TRUE),
          p_value = NA_real_, p_adj = NA_real_, test = "wilcox_unpaired"
        ))
      }
      
      wt <- suppressWarnings(stats::wilcox.test(xa, xb, paired = FALSE, exact = FALSE))
      tibble::tibble(
        Contrast = label, Analyte = a,
        n_A = length(xa), n_B = length(xb),
        median_A = stats::median(xa, na.rm = TRUE),
        median_B = stats::median(xb, na.rm = TRUE),
        median_diff_AminusB = stats::median(xa, na.rm = TRUE) - stats::median(xb, na.rm = TRUE),
        p_value = wt$p.value, test = "wilcox_unpaired"
      )
    }) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"))
  
  by_analyte
}


## 3) Run the four comparisons

# 3.1 Infants: 12 vs 44 wks (paired)
cmp1_dir <- file.path(out_root, "Infants_12vs44_paired")
dir.create(cmp1_dir, showWarnings = FALSE, recursive = TRUE)

stats_12v44 <- plasma_test %>%
  filter(Group == "Infant", Timepoint %in% c("12 Wks","44 Wks")) %>%
  paired_stats_12v44()

write_csv(stats_12v44, file.path(cmp1_dir, "stats_Infants_12vs44_paired.csv"))

# 3.2 Infants: by Maternal IGRA (Pos vs Neg, unpaired)
cmp2_dir <- file.path(out_root, "Infants_byMaternalIGRA")
dir.create(cmp2_dir, showWarnings = FALSE, recursive = TRUE)

stats_inf_by_mIGRA <- plasma_test %>%
  filter(Group == "Infant") %>%
  unpaired_stats_twolevel(grp_var = igra_m_col, lvl_a = "Pos", lvl_b = "Neg", label = "MaternalIGRA_Pos_vs_Neg")

write_csv(stats_inf_by_mIGRA, file.path(cmp2_dir, "stats_Infants_byMaternalIGRA.csv"))

# 3.3 Infants: by Infant IGRA (Pos vs Neg, unpaired)
cmp3_dir <- file.path(out_root, "Infants_byInfantIGRA")
dir.create(cmp3_dir, showWarnings = FALSE, recursive = TRUE)

stats_inf_by_iIGRA <- plasma_test %>%
  filter(Group == "Infant") %>%
  unpaired_stats_twolevel(grp_var = igra_i_col, lvl_a = "Pos", lvl_b = "Neg", label = "InfantIGRA_Pos_vs_Neg")

write_csv(stats_inf_by_iIGRA, file.path(cmp3_dir, "stats_Infants_byInfantIGRA.csv"))

# 3.4 Mothers: by Maternal IGRA (Pos vs Neg, unpaired)
cmp4_dir <- file.path(out_root, "Mothers_byMaternalIGRA")
dir.create(cmp4_dir, showWarnings = FALSE, recursive = TRUE)

stats_moms_by_mIGRA <- plasma_test %>%
  filter(Group == "Mother") %>%
  unpaired_stats_twolevel(grp_var = igra_m_col, lvl_a = "Pos", lvl_b = "Neg", label = "MaternalIGRA_Pos_vs_Neg")

write_csv(stats_moms_by_mIGRA, file.path(cmp4_dir, "stats_Mothers_byMaternalIGRA.csv"))

## 4) Quick completion message
message("Step B complete. CSVs written to:\n",
        " - ", cmp1_dir, "\n",
        " - ", cmp2_dir, "\n",
        " - ", cmp3_dir, "\n",
        " - ", cmp4_dir)

# ---- Step C1: plotting setup ----

# Root for plots (Windows)
plot_root <- "C:/Users/ammas/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/Plasma_Plots"
if (!dir.exists(plot_root)) dir.create(plot_root, recursive = TRUE)

# Subfolders per comparison
dir_cmp_12v44   <- file.path(plot_root, "Infants_12vs44_paired")
dir_cmp_i_mIGRA <- file.path(plot_root, "Infants_byMaternalIGRA")
dir_cmp_i_iIGRA <- file.path(plot_root, "Infants_byInfantIGRA")
dir_cmp_m_mIGRA <- file.path(plot_root, "Mothers_byMaternalIGRA")
for (d in c(dir_cmp_12v44, dir_cmp_i_mIGRA, dir_cmp_i_iIGRA, dir_cmp_m_mIGRA)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Helper: safe filenames
sanitize_for_path <- function(x) {
  x %>%
    str_replace_all("[\\/:*?\"<>|]", "_") %>%
    str_squish()
}

# Helper: compact subtitle with n and (optional) stats
mk_subtitle <- function(n_a = NA, n_b = NA, n_pairs = NA, p = NA, fdr = NA) {
  bits <- c()
  if (!is.na(n_pairs)) bits <- c(bits, glue("n_pairs = {n_pairs}"))
  if (!is.na(n_a) && !is.na(n_b)) bits <- c(bits, glue("n = {n_a} vs {n_b}"))
  if (!is.na(p))   bits <- c(bits, glue("p = {signif(p, 3)}"))
  if (!is.na(fdr)) bits <- c(bits, glue("FDR = {signif(fdr, 3)}"))
  paste(bits, collapse = " | ")
}

# A clean default theme
theme_plasma <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.x = element_blank(),
      legend.position = "none"
    )
}

# Geoms we’ll re-use
geom_summary_crossbar <- function(width = 0.5) {
  stat_summary(fun = median, geom = "crossbar", width = width, fatten = 0, alpha = 0.8)
}

# We will plot on Z-scores for comparability across analytes and timepoints
plasma_plot_data <- plasma_work %>%
  filter(Flag == "OK", !is.na(Value_z))
# ---- Step C2: Infants 12 vs 44 weeks (paired) ----
analytes <- sort(unique(plasma_plot_data$Analyte))

purrr::walk(analytes, function(ana) {
  df <- plasma_plot_data %>%
    dplyr::filter(
      Group == "Infant",
      TP %in% c("12 Wks", "44 Wks"),
      Analyte == !!ana
    ) %>%
    dplyr::select(PID, TP, Analyte, Value_z) %>%
    dplyr::distinct()

  # Build paired matrix
  wide <- tidyr::pivot_wider(df, names_from = TP, values_from = Value_z)
  wide <- dplyr::filter(wide, !is.na(`12 Wks`), !is.na(`44 Wks`))

  n_pairs <- nrow(wide)
  if (n_pairs < 3) {
    message(glue::glue("Skip {ana}: not enough pairs (n_pairs={n_pairs})"))
    return(invisible(NULL))
  }

  # Paired Wilcoxon (44 - 12)
  wt <- suppressWarnings(stats::wilcox.test(wide$`44 Wks`, wide$`12 Wks`, paired = TRUE, exact = FALSE))
  med_delta <- stats::median(wide$`44 Wks` - wide$`12 Wks`, na.rm = TRUE)

  # Rebuild long for plotting
  df_plot <- wide %>%
    tidyr::pivot_longer(cols = c("12 Wks", "44 Wks"), names_to = "TP", values_to = "Value_z") %>%
    dplyr::mutate(TP = factor(TP, levels = c("12 Wks","44 Wks")))

  # Dynamic width by number of pairs
  w_in <- max(4.8, min(7.5, 4.8 + 0.05 * n_pairs))
  h_in <- 4.8

  p <- ggplot(df_plot, aes(x = TP, y = Value_z, group = PID)) +
    geom_line(alpha = 0.35) +
    geom_point(aes(color = TP), size = 2) +
    geom_summary_crossbar() +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = glue::glue("{ana} — Infants: 12 vs 44 weeks (paired)"),
      subtitle = mk_subtitle(n_pairs = n_pairs, p = wt$p.value),
      x = NULL,
      y = "Z-score (standardized within timepoint)"
    ) +
    theme_plasma() +
    theme(
      plot.title = element_text(size = 13, face = "bold", margin = margin(b = 6)),
      plot.subtitle = element_text(size = 10, color = "gray25")
    )
  

  out_name <- file.path(
    dir_cmp_12v44,
    paste0(sanitize_for_path(glue::glue("{ana}__Infants_12vs44_paired")), ".png")
  )
  ggsave(out_name, p, width = w_in, height = h_in, dpi = 300, bg = "white")
})
#### C2b) Infants — Maternal IGRA (Pos vs Neg), unpaired, plotted per analyte × timepoint
# ---- paths ----
cmp2_plotdir <- file.path(out_root, "Infants_byMaternalIGRA", "plots")
if (!dir.exists(cmp2_plotdir)) dir.create(cmp2_plotdir, recursive = TRUE)

# ---- helper (local) ----
mk_subtitle <- function(nA, nB, p, fdr = NA_real_) {
  glue::glue("n_Pos = {nA} | n_Neg = {nB} | p = {signif(p, 3)}{ifelse(is.na(fdr),'', paste0(' | FDR = ', signif(fdr,3)))}")
}

# ---- plotting loop ----
analytes <- sort(unique(plasma_test$Analyte))
timepoints <- c("12 Wks","44 Wks")

purrr::walk(analytes, function(ana) {
  purrr::walk(timepoints, function(tp) {
    df <- plasma_test %>%
      dplyr::filter(Group == "Infant",
                    Timepoint == tp,
                    Analyte == ana,
                    !is.na(.data[[igra_m_col]]),
                    .data[[igra_m_col]] %in% c("Pos","Neg"))
    
    if (nrow(df) < 4) return(invisible(NULL))
    
    # Wilcoxon (unpaired)
    xa <- df %>% dplyr::filter(.data[[igra_m_col]] == "Pos") %>% dplyr::pull(Y)
    xb <- df %>% dplyr::filter(.data[[igra_m_col]] == "Neg") %>% dplyr::pull(Y)
    if (length(xa) < 2 || length(xb) < 2) return(invisible(NULL))
    wt <- suppressWarnings(wilcox.test(xa, xb, paired = FALSE, exact = FALSE))
    
    # dynamic width by sample size
    nA <- length(xa); nB <- length(xb)
    w_in <- max(4.8, min(7, 4.5 + 0.02 * (nA + nB)))
    h_in <- 4.5
    
    p <- ggplot(df, aes(x = .data[[igra_m_col]], y = Y)) +
      geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.6) +
      ggplot2::geom_point(
        position = position_jitter(width = 0.1, height = 0),
        size = 1.9, alpha = 0.75, aes(color = .data[[igra_m_col]])
      ) +
      geom_summary_crossbar() +
      scale_x_discrete(limits = c("Neg","Pos"), drop = FALSE) +
      labs(
        title = glue::glue("{ana} — Infants ({tp}) — Maternal IGRA"),
        subtitle = mk_subtitle(nA, nB, wt$p.value),
        x = NULL,
        y = "log1p(pg/ml) (per-TP standardized input)"
      ) +
      theme_plasma() +
      theme(
        plot.title = element_text(size = 13, face = "bold", margin = margin(b = 6)),
        plot.subtitle = element_text(size = 10, color = "gray25"),
        legend.position = "none"
      )
    
    fname <- file.path(
      cmp2_plotdir,
      sanitize_for_path(glue::glue("{ana}__Infants__{tp}__MaternalIGRA.png"))
    )
    ggsave(fname, p, width = w_in, height = h_in, dpi = 300, bg = "white")
  })
})

### C2c) Infants — Infant IGRA (Pos vs Neg), unpaired, plotted per analyte × timepoint
# ---- paths ----
cmp3_plotdir <- file.path(out_root, "Infants_byInfantIGRA", "plots")
if (!dir.exists(cmp3_plotdir)) dir.create(cmp3_plotdir, recursive = TRUE)

purrr::walk(analytes, function(ana) {
  purrr::walk(timepoints, function(tp) {
    df <- plasma_test %>%
      dplyr::filter(Group == "Infant",
                    Timepoint == tp,
                    Analyte == ana,
                    !is.na(.data[[igra_i_col]]),
                    .data[[igra_i_col]] %in% c("Pos","Neg"))
    
    if (nrow(df) < 4) return(invisible(NULL))
    
    xa <- df %>% dplyr::filter(.data[[igra_i_col]] == "Pos") %>% dplyr::pull(Y)
    xb <- df %>% dplyr::filter(.data[[igra_i_col]] == "Neg") %>% dplyr::pull(Y)
    if (length(xa) < 2 || length(xb) < 2) return(invisible(NULL))
    wt <- suppressWarnings(wilcox.test(xa, xb, paired = FALSE, exact = FALSE))
    
    nA <- length(xa); nB <- length(xb)
    w_in <- max(4.8, min(7, 4.5 + 0.02 * (nA + nB)))
    h_in <- 4.5
    
    p <- ggplot(df, aes(x = .data[[igra_i_col]], y = Y)) +
      geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.6) +
      ggplot2::geom_point(
        position = position_jitter(width = 0.1, height = 0),
        size = 1.9, alpha = 0.75, aes(color = .data[[igra_i_col]])
      ) +
      geom_summary_crossbar() +
      scale_x_discrete(limits = c("Neg","Pos"), drop = FALSE) +
      labs(
        title = glue::glue("{ana} — Infants ({tp}) — Infant IGRA"),
        subtitle = mk_subtitle(nA, nB, wt$p.value),
        x = NULL,
        y = "log1p(pg/ml) (per-TP standardized input)"
      ) +
      theme_plasma() +
      theme(
        plot.title = element_text(size = 13, face = "bold", margin = margin(b = 6)),
        plot.subtitle = element_text(size = 10, color = "gray25"),
        legend.position = "none"
      )
    
    fname <- file.path(
      cmp3_plotdir,
      sanitize_for_path(glue::glue("{ana}__Infants__{tp}__InfantIGRA.png"))
    )
    ggsave(fname, p, width = w_in, height = h_in, dpi = 300, bg = "white")
  })
})


### C2d) Mothers — Maternal IGRA (Pos vs Neg), unpaired, plotted per analyte (Entry only)
# ---- paths ----
cmp4_plotdir <- file.path(out_root, "Mothers_byMaternalIGRA", "plots")
if (!dir.exists(cmp4_plotdir)) dir.create(cmp4_plotdir, recursive = TRUE)

purrr::walk(analytes, function(ana) {
  df <- plasma_test %>%
    dplyr::filter(Group == "Mother",
                  Timepoint == "Entry",
                  Analyte == ana,
                  !is.na(.data[[igra_m_col]]),
                  .data[[igra_m_col]] %in% c("Pos","Neg"))
  
  if (nrow(df) < 4) return(invisible(NULL))
  
  xa <- df %>% dplyr::filter(.data[[igra_m_col]] == "Pos") %>% dplyr::pull(Y)
  xb <- df %>% dplyr::filter(.data[[igra_m_col]] == "Neg") %>% dplyr::pull(Y)
  if (length(xa) < 2 || length(xb) < 2) return(invisible(NULL))
  wt <- suppressWarnings(wilcox.test(xa, xb, paired = FALSE, exact = FALSE))
  
  nA <- length(xa); nB <- length(xb)
  w_in <- max(4.8, min(7, 4.5 + 0.02 * (nA + nB)))
  h_in <- 4.5
  
  p <- ggplot(df, aes(x = .data[[igra_m_col]], y = Y)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.6) +
    ggplot2::geom_point(
      position = position_jitter(width = 0.1, height = 0),
      size = 1.9, alpha = 0.75, aes(color = .data[[igra_m_col]])
    ) +
    geom_summary_crossbar() +
    scale_x_discrete(limits = c("Neg","Pos"), drop = FALSE) +
    labs(
      title = glue::glue("{ana} — Mothers (Entry) — Maternal IGRA"),
      subtitle = mk_subtitle(nA, nB, wt$p.value),
      x = NULL,
      y = "log1p(pg/ml) (per-TP standardized input)"
    ) +
    theme_plasma() +
    theme(
      plot.title = element_text(size = 13, face = "bold", margin = margin(b = 6)),
      plot.subtitle = element_text(size = 10, color = "gray25"),
      legend.position = "none"
    )
  
  fname <- file.path(
    cmp4_plotdir,
    sanitize_for_path(glue::glue("{ana}__Mothers__Entry__MaternalIGRA.png"))
  )
  ggsave(fname, p, width = w_in, height = h_in, dpi = 300, bg = "white")
})
