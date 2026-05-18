## ============================================================================
## Project: Flow Cytometry and Plasma Biomarkers on Mothers and Infants (P1078)
## Script:  1_Data_Cleaning_Norm.R
## Purpose: Load TB Flow Cytokine Excel, parse gate paths, baseline-subtract
##          to MED, batch-normalize on a logit scale, and save cleaned long
##          tables for downstream analyses.
##
## Outputs: saved_R_dat/TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.{qs2,csv}
##
## Key fixes vs prior version:
##   - Metadata vs flow columns identified by NAME, not by hard-coded position.
##   - Constants `ctrl_min_pct`, `logit_clamp`, `shift_cap` defined once with
##     unambiguous names (no more redefinition of `min_ctrl_pct` / `tiny_eps`).
##   - MED baseline uses median (robust to batch-control duplicates like HC-SP
##     which has 10 MED rows because it runs once per batch).
##   - Delta is NA when MED is missing (no fallback to raw % — those are
##     semantically different quantities).
##   - `Value_Normalized` renamed to `Delta_unscaled` so the variable name
##     no longer suggests normalization that isn't applied.
## ============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr)
  library(readr); library(purrr); library(qs2)
})

## ---- Paths ----------------------------------------------------------------
proj_root <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/Flow Cytometry and Plasma Biomarkers on Mothers and Infants (P1078)"
raw_dir   <- file.path(proj_root, "Raw_Data")
out_dir   <- file.path(proj_root, "saved_R_dat")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

file_xlsx <- file.path(raw_dir, "TB_Flow_Cytokine_DATA_102725.xlsx")
stopifnot(file.exists(file_xlsx))

## ---- Read & clean character missings --------------------------------------
cytokine_data_raw <- read_excel(file_xlsx, .name_repair = "minimal") %>%
  mutate(across(where(is.character), ~ {
    v    <- trimws(.x)
    miss <- tolower(v) %in% c("na","n/a","n.a.","nan","null","missing","--","-")
    v[miss | v == ""] <- NA_character_
    v
  })) %>%
  type_convert()

## ---- Identify metadata vs flow columns by NAME ----------------------------
meta_names <- c("File","Batch","PID","PID Type","Time Point",
                "Maternal IGRA","Maternal PID","Infant IGRA","Infant PID","Condition")

missing_meta <- setdiff(meta_names, colnames(cytokine_data_raw))
if (length(missing_meta) > 0) {
  stop("Expected metadata columns missing from input: ",
       paste(missing_meta, collapse = ", "))
}
flow_cols <- setdiff(colnames(cytokine_data_raw), meta_names)
message("Detected ", length(flow_cols), " flow gate columns and ",
        length(meta_names), " metadata columns.")

## ---- Pivot to long --------------------------------------------------------
cytokine_long <- cytokine_data_raw %>%
  pivot_longer(cols = all_of(flow_cols),
               names_to  = "GatePath",
               values_to = "Value") %>%
  mutate(Value = if (is.numeric(Value)) Value else suppressWarnings(
    as.numeric(str_replace_all(as.character(Value), "[^0-9eE\\.\\-+]", ""))
  ))

## ---------------------------------------------------------------------------
## Step 4 — Gate-path parsing (Prolif status, T1..Tn, parent/grandparent, etc.)
## ---------------------------------------------------------------------------
cytokine_long <- cytokine_long %>%
  mutate(
    GatePath_norm = str_replace_all(GatePath, "[\\u00A0\\p{Zs}]+", " "),
    CorePath      = str_replace(GatePath_norm, "\\s*\\|\\s*Freq\\..*$", "")
  ) %>%
  mutate(
    parsed = map(CorePath, function(cp) {
      if (is.na(cp)) return(list(
        Prolif_Status = NA_character_, T_levels = character(0),
        ParentGate = NA_character_, GrandparentGate = NA_character_,
        Lineage = NA_character_, Compartment = NA_character_
      ))
      
      segs <- str_squish(str_replace_all(str_split(cp, "/", simplify = TRUE), ",", " "))
      segs <- segs[segs != ""]
      
      idx_p <- which(segs %in% c("Proliferating","Non-Proliferating"))
      has_p <- length(idx_p) > 0
      idx   <- if (has_p) idx_p[1] else NA_integer_
      
      T_levels <- if (has_p && idx < length(segs)) {
        str_split(paste(segs[(idx+1):length(segs)], collapse = "/"), "\\s*/\\s*")[[1]] %>%
          str_squish()
      } else character(0)
      k <- length(T_levels)
      
      if (has_p) {
        Prolif_Status <- segs[idx]
        PreProlif     <- if (idx > 1) segs[idx - 1] else NA_character_
        if (k >= 2) {
          ParentGate      <- T_levels[k - 1]
          GrandparentGate <- if (k >= 3) T_levels[k - 2] else Prolif_Status
        } else if (k == 1) {
          ParentGate      <- Prolif_Status
          GrandparentGate <- PreProlif
        } else {
          ParentGate      <- Prolif_Status
          GrandparentGate <- PreProlif
        }
      } else {
        Prolif_Status <- NA_character_
        if (length(segs) >= 2) {
          ParentGate      <- segs[length(segs) - 1]
          GrandparentGate <- segs[length(segs) - 2]
        } else {
          ParentGate      <- if (length(segs) >= 1) segs[length(segs)] else NA_character_
          GrandparentGate <- NA_character_
        }
      }
      
      Lineage <- dplyr::case_when(
        any(segs == "gd TCR-") ~ "gd TCR-",
        any(segs == "gd TCR")  ~ "gd TCR",
        TRUE                   ~ NA_character_
      )
      Compartment <- dplyr::case_when(
        any(segs == "CD4")    ~ "CD4",
        any(segs == "CD8")    ~ "CD8",
        any(segs == "gd TCR") ~ "gdTCR",
        TRUE                  ~ NA_character_
      )
      
      list(Prolif_Status = Prolif_Status, T_levels = T_levels,
           ParentGate = ParentGate, GrandparentGate = GrandparentGate,
           Lineage = Lineage, Compartment = Compartment)
    })
  ) %>%
  unnest_wider(parsed)

# Expand T_levels into fixed T1..Tn columns
max_levels <- cytokine_long %>%
  mutate(n_levels = map_int(T_levels, length)) %>%
  pull(n_levels) %>% max(na.rm = TRUE)

cytokine_long <- cytokine_long %>%
  mutate(T_levels = map(T_levels, ~ c(.x, rep(NA_character_, max_levels - length(.x))))) %>%
  unnest_wider(T_levels, names_sep = "") %>%
  rename_with(~ paste0("T", seq_along(.)), dplyr::starts_with("T_levels")) %>%
  mutate(across(starts_with("T"), ~ ifelse(is.na(.x), NA_character_, str_squish(.x)))) %>%
  mutate(
    Metric = case_when(
      str_detect(GatePath_norm, "\\|\\s*Freq\\.\\s*of\\s*Parent")      ~ "Freq. of Parent",
      str_detect(GatePath_norm, "\\|\\s*Freq\\.\\s*of\\s*Grandparent") ~ "Freq. of Grandparent",
      TRUE                                                              ~ NA_character_
    ),
    DenominatorGate = case_when(
      Metric == "Freq. of Parent"      ~ ParentGate,
      Metric == "Freq. of Grandparent" ~ GrandparentGate,
      TRUE                              ~ NA_character_
    ),
    GateDepth = str_count(CorePath, "/")
  ) %>%
  select(-GatePath_norm, -CorePath)

## ---------------------------------------------------------------------------
## Step 5 — Grouping keys, drop File column
## ---------------------------------------------------------------------------
cytokine_long <- cytokine_long %>%
  mutate(
    PID_Type   = `PID Type`,
    Time_Point = `Time Point`,
    Batch      = as.character(Batch),
    FamilyID   = dplyr::coalesce(as.character(`Maternal PID`), as.character(PID))
  ) %>%
  select(-File)

## ---------------------------------------------------------------------------
## Step 6 — Baseline subtraction (Stim − MED) within PID × GatePath × TP
##   Uses MEDIAN so duplicate MED rows (e.g., HC-SP batch control runs once
##   per batch → 10 MED rows for that PID/TP combo) don't get pulled by outliers.
## ---------------------------------------------------------------------------
med_table <- cytokine_long %>%
  dplyr::filter(Condition == "MED") %>%
  dplyr::group_by(PID, PID_Type, Time_Point, GatePath) %>%
  dplyr::summarise(MED_Value = stats::median(Value, na.rm = TRUE), .groups = "drop")

cytokine_long <- cytokine_long %>%
  dplyr::left_join(med_table, by = c("PID","PID_Type","Time_Point","GatePath")) %>%
  dplyr::mutate(
    Value_BaselineAdj = dplyr::if_else(!is.na(MED_Value), Value - MED_Value, NA_real_),
    MED_Missing       = is.na(MED_Value)
  )

# Zero-like flags (raw + MED), purely diagnostic
zerolike_eps <- 0.01
cytokine_long <- cytokine_long %>%
  dplyr::mutate(
    ZeroLike_Raw = !is.na(Value)     & Value     < zerolike_eps,
    ZeroLike_MED = !is.na(MED_Value) & MED_Value < zerolike_eps
  )

## ---------------------------------------------------------------------------
## Step 7 — Batch normalization on the logit scale
##   - Reference per Batch × GatePath is the median of MED-condition batch
##     controls (HC-SP), but only when that median is ≥ ctrl_min_pct and ≤ 99.9%.
##   - Shift is logit(batch_ref) − logit(reference_batch_ref), capped at ±shift_cap.
##   - Applied to absolute values only; deltas are left untouched.
## ---------------------------------------------------------------------------
ctrl_min_pct <- 1.0     # ignore controls <1% when computing reference medians
logit_clamp  <- 1e-4    # clamp for logit transform
shift_cap    <- 2.0     # cap absolute shift on logit scale (~7.4× odds)

logit_safe   <- function(p) qlogis(pmin(pmax(p, logit_clamp), 1 - logit_clamp))
invlogit_pct <- function(x) plogis(x) * 100

# 7.1 — Per Batch × GatePath median of MED controls (filtered)
ctrl_meds_med <- cytokine_long %>%
  dplyr::filter(PID_Type == "Batch Control", Condition == "MED") %>%
  dplyr::mutate(Value_use = dplyr::if_else(Value >= ctrl_min_pct, Value, NA_real_)) %>%
  dplyr::group_by(Batch, GatePath) %>%
  dplyr::summarise(ctrl_med_med = stats::median(Value_use, na.rm = TRUE), .groups = "drop")

# Reference batch: prefer "Batch 1", else first alphabetical
REF_BATCH <- if ("Batch 1" %in% unique(ctrl_meds_med$Batch)) "Batch 1" else sort(unique(ctrl_meds_med$Batch))[1]
message("Using reference batch for normalization: ", REF_BATCH)

ref_ctrl_med <- ctrl_meds_med %>%
  dplyr::filter(Batch == REF_BATCH) %>%
  dplyr::select(GatePath, ref_ctrl_med = ctrl_med_med)

# 7.2 — Compute capped logit shifts per Batch × GatePath
ctrl_shifts <- ctrl_meds_med %>%
  dplyr::left_join(ref_ctrl_med, by = "GatePath") %>%
  dplyr::mutate(
    ok_batch  = !is.na(ctrl_med_med) & ctrl_med_med >= ctrl_min_pct & ctrl_med_med <= 99.9,
    ok_ref    = !is.na(ref_ctrl_med) & ref_ctrl_med  >= ctrl_min_pct & ref_ctrl_med  <= 99.9,
    p_b       = ctrl_med_med / 100,
    p_r       = ref_ctrl_med / 100,
    logit_b   = dplyr::if_else(ok_batch, logit_safe(p_b), NA_real_),
    logit_r   = dplyr::if_else(ok_ref,   logit_safe(p_r), NA_real_),
    raw_shift = logit_b - logit_r,
    shift     = dplyr::case_when(
      is.finite(raw_shift) ~ pmax(pmin(raw_shift, shift_cap), -shift_cap),
      TRUE                 ~ NA_real_
    ),
    Shift_Reason = dplyr::case_when(
      ok_batch & ok_ref ~ "shifted:MED logit",
      !ok_ref           ~ "skip:ref<1% or >99%",
      !ok_batch         ~ "skip:batch<1% or >99%",
      TRUE              ~ "skip:other"
    )
  ) %>%
  dplyr::select(Batch, GatePath, shift, Shift_Reason)

# 7.3 — Apply: absolute values get the logit shift; deltas pass through
cytokine_long <- cytokine_long %>%
  dplyr::left_join(ctrl_shifts, by = c("Batch", "GatePath")) %>%
  dplyr::mutate(
    Value_RawNormalized = dplyr::case_when(
      is.na(Value)               ~ NA_real_,
      !is.finite(shift)          ~ Value,                        # no shift info: keep raw
      Value < ctrl_min_pct       ~ Value,                        # don't "correct" noise
      TRUE ~ {
        p <- pmin(pmax(Value / 100, logit_clamp), 1 - logit_clamp)
        invlogit_pct(logit_safe(p) - shift)
      }
    ),
    Delta_unscaled = Value_BaselineAdj,   # NA when MED is missing (no fallback)
    Shift_Used     = is.finite(shift) & Value >= ctrl_min_pct
  )

## ---------------------------------------------------------------------------
## Step 7d — Diagnostics
## ---------------------------------------------------------------------------
summary_shift <- cytokine_long %>%
  mutate(is_mid = !is.na(Value) & Value > 1 & Value < 99) %>%
  summarise(
    n_total       = n(),
    n_mid         = sum(is_mid),
    n_shift_used  = sum(is_mid & Shift_Used, na.rm = TRUE),
    frac_shifted  = n_shift_used / pmax(n_mid, 1)
  )
print(summary_shift)

ref_saturated <- ref_ctrl_med %>%
  mutate(ref_is_sat = ref_ctrl_med >= 99 | ref_ctrl_med < 1) %>%
  arrange(desc(ref_is_sat), desc(ref_ctrl_med)) %>%
  head(15)
message("Reference-batch saturation check (top 15):")
print(ref_saturated)

# Delta missingness audit
delta_audit <- cytokine_long %>%
  dplyr::filter(!Condition %in% c("MED"), PID_Type != "Batch Control") %>%
  dplyr::summarise(
    n_rows         = n(),
    n_delta_NA     = sum(is.na(Delta_unscaled)),
    n_MED_missing  = sum(MED_Missing, na.rm = TRUE)
  )
message("Delta missingness in stim conditions (mothers/infants):")
print(delta_audit)

## ---------------------------------------------------------------------------
## Step 8 — Transforms for modeling
## ---------------------------------------------------------------------------
logit_eps <- 1e-4
cytokine_long <- cytokine_long %>%
  dplyr::mutate(
    prop_abs    = pmin(pmax(Value_RawNormalized / 100, logit_eps), 1 - logit_eps),
    Value_Logit = qlogis(prop_abs)
  ) %>%
  dplyr::select(-prop_abs) %>%
  dplyr::mutate(Delta_for_model = Delta_unscaled)  # NA-preserving alias

## ---------------------------------------------------------------------------
## Step 9 — Save (qs2 + CSV)
## ---------------------------------------------------------------------------
qs2::qs_save(
  cytokine_long,
  file.path(out_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.qs2")
)

cytokine_long_csv <- cytokine_long %>%
  dplyr::mutate(dplyr::across(where(is.list), ~ sapply(.x, toString)))
readr::write_csv(
  cytokine_long_csv,
  file.path(out_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.csv")
)

message("Done. Wrote cleaned cytokine data to: ", out_dir)
