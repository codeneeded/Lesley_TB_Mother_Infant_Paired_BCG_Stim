# install.packages(c("readxl","dplyr","tidyr","stringr","readr")) # if needed
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
#######3 Set directories
data_dir  <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/Raw_Data"
file_xlsx <- file.path(data_dir, "TB_Flow_Cytokine DATA_102725.xlsx")

############ Read and Clean dataset
# Keep names exactly as in Excel; fix 'NA' style strings only in character columns
cytokine_data_raw <- read_excel(file_xlsx, .name_repair = "minimal") %>%
  mutate(
    across(
      where(is.character),
      ~{
        v <- trimws(.x)
        # treat common typed-missing variants as NA
        miss <- tolower(v) %in% c("na","n/a","n.a.","nan","null","missing","--","-")
        v[miss | v == ""] <- NA_character_
        v
      }
    )
  ) %>%
  # Re-parse character columns to numeric/dates where possible
  type_convert()

# Indices: first 10 are metadata; 11+ are flow frequency readouts
meta_cols <- colnames(cytokine_data_raw)[1:10]
flow_cols <- colnames(cytokine_data_raw)[11:length(colnames(cytokine_data_raw))]

cytokine_long <- cytokine_data_raw %>%
  pivot_longer(
    cols = all_of(flow_cols),
    names_to = "GatePath",
    values_to = "Value"
  )

# If any readouts arrived as character like "<0.01" or "3%", strip non-numeric and coerce
cytokine_long <- cytokine_long %>%
  mutate(
    Value = if (is.numeric(Value)) Value else suppressWarnings(
      as.numeric(str_replace_all(as.character(Value), "[^0-9eE\\.\\-+]", ""))
    )
  )
# ---------------------------
# Step 4 — Robust gate parsing
# ---------------------------

# ---- Reparse Step 4 with correct Parent/Grandparent relative to final gate ----
cytokine_long <- cytokine_long %>%
  mutate(
    GatePath_norm = str_replace_all(GatePath, "[\\u00A0\\p{Zs}]+", " "),
    CorePath = str_replace(GatePath_norm, "\\s*\\|\\s*Freq\\..*$", "")
  ) %>%
  mutate(
    parsed = map(CorePath, function(cp) {
      if (is.na(cp)) return(list(
        Prolif_Status = NA_character_, T_levels = character(0),
        ParentGate = NA_character_, GrandparentGate = NA_character_,
        Lineage = NA_character_, Compartment = NA_character_
      ))
      
      # split on "/" into segments; within segments, treat commas as spaces
      segs <- str_squish(str_replace_all(str_split(cp, "/", simplify = TRUE), ",", " "))
      segs <- segs[segs != ""]
      
      # locate proliferating status
      idx_p <- which(segs %in% c("Proliferating","Non-Proliferating"))
      has_p <- length(idx_p) > 0
      idx <- if (has_p) idx_p[1] else NA_integer_
      
      # levels after prolif (T1..Tk)
      T_levels <- if (has_p && idx < length(segs)) {
        str_split(paste(segs[(idx+1):length(segs)], collapse = "/"), "\\s*/\\s*")[[1]] %>% str_squish()
      } else character(0)
      
      k <- length(T_levels)
      
      # structural parent/grandparent relative to the *final* gate
      if (has_p) {
        Prolif_Status <- segs[idx]
        # the segment before Prolif (e.g., CD4) if it exists
        PreProlif <- if (idx > 1) segs[idx - 1] else NA_character_
        
        if (k >= 2) {
          ParentGate      <- T_levels[k - 1]                    # T(k-1)
          GrandparentGate <- if (k >= 3) T_levels[k - 2] else Prolif_Status
        } else if (k == 1) {
          ParentGate      <- Prolif_Status                      # immediate parent of T1
          GrandparentGate <- PreProlif                          # e.g., CD4
        } else { # k == 0 (no gate after Prolif)
          ParentGate      <- Prolif_Status
          GrandparentGate <- PreProlif
        }
      } else {
        # no Prolif segment: use last two segments of the path
        Prolif_Status <- NA_character_
        if (length(segs) >= 2) {
          ParentGate      <- segs[length(segs) - 1]
          GrandparentGate <- segs[length(segs) - 2]
        } else {
          ParentGate      <- if (length(segs) >= 1) segs[length(segs)] else NA_character_
          GrandparentGate <- NA_character_
        }
      }
      
      # lineage / compartment hints from the whole path
      Lineage <- dplyr::case_when(
        any(segs == "gd TCR-") ~ "gd TCR-",
        any(segs == "gd TCR")  ~ "gd TCR",
        TRUE ~ NA_character_
      )
      Compartment <- dplyr::case_when(
        any(segs == "CD4")    ~ "CD4",
        any(segs == "CD8")    ~ "CD8",
        any(segs == "gd TCR") ~ "gdTCR",
        TRUE ~ NA_character_
      )
      
      list(
        Prolif_Status   = Prolif_Status,
        T_levels        = T_levels,
        ParentGate      = ParentGate,
        GrandparentGate = GrandparentGate,
        Lineage         = Lineage,
        Compartment     = Compartment
      )
    })
  ) %>%
  unnest_wider(parsed)

# build T1..Tn
max_levels <- cytokine_long %>%
  mutate(n_levels = map_int(T_levels, length)) %>%
  pull(n_levels) %>% max(na.rm = TRUE)

cytokine_long <- cytokine_long %>%
  mutate(T_levels = map(T_levels, ~ c(.x, rep(NA_character_, max_levels - length(.x))))) %>%
  unnest_wider(T_levels, names_sep = "") %>%
  rename_with(~ paste0("T", seq_along(.)), dplyr::starts_with("T_levels")) %>%
  mutate(across(starts_with("T"), ~ ifelse(is.na(.x), NA_character_, str_squish(.x)))) %>%
  # metric & explicit denominator gate
  mutate(
    Metric = case_when(
      str_detect(GatePath_norm, "\\|\\s*Freq\\.\\s*of\\s*Parent") ~ "Freq. of Parent",
      str_detect(GatePath_norm, "\\|\\s*Freq\\.\\s*of\\s*Grandparent") ~ "Freq. of Grandparent",
      TRUE ~ NA_character_
    ),
    DenominatorGate = case_when(
      Metric == "Freq. of Parent"      ~ ParentGate,
      Metric == "Freq. of Grandparent" ~ GrandparentGate,
      TRUE ~ NA_character_
    ),
    GateDepth = str_count(CorePath, "/")
  ) %>%
  select(-GatePath_norm, -CorePath)
# ---------------------------
# Step 5 — Build grouping keys, drop File
# ---------------------------
cytokine_long <- cytokine_long %>%
  mutate(
    PID_Type   = `PID Type`,
    Time_Point = `Time Point`,
    Batch      = as.character(Batch),
    # Prefer Maternal PID (character) when present, else PID
    FamilyID   = dplyr::coalesce(as.character(`Maternal PID`), as.character(PID))
  ) %>%
  select(-File)

# ---------------------------
# Step 6 — Baseline adjust to MED within PID × Gate (ΔStim = Stim − MED)
# (Uses GatePath so the exact gate/denominator is respected)
# ---------------------------
med_table <- cytokine_long %>%
  filter(Condition == "MED") %>%
  group_by(PID, PID_Type, Time_Point, GatePath) %>%
  summarize(MED_Value = mean(Value, na.rm = TRUE), .groups = "drop")

cytokine_long <- cytokine_long %>%
  left_join(med_table, by = c("PID","PID_Type","Time_Point","GatePath")) %>%
  mutate(Value_BaselineAdj = ifelse(!is.na(MED_Value), Value - MED_Value, NA_real_))

# ---------------------------
# Step 7 — Batch normalization using HC-SP @ MED ONLY
# One scale factor per Batch × GatePath; apply to all conditions in that batch
# ---------------------------

# Control summary (Batch Control, MED, raw Value) Compute the mean of the batch control (HC-SP, MED) per batch × gate.
ctrl_means_med <- cytokine_long %>%
  filter(PID_Type == "Batch Control", Condition == "MED") %>%
  group_by(Batch, GatePath) %>%
  summarize(ctrl_mean_med = mean(Value, na.rm = TRUE), .groups = "drop")

# Choose a reference batch robustly: prefer "Batch 1" if present; else the first in sorted order
REF_BATCH <- if ("Batch 1" %in% unique(ctrl_means_med$Batch)) "Batch 1" else sort(unique(ctrl_means_med$Batch))[1]

ref_ctrl_med <- ctrl_means_med %>%
  filter(Batch == REF_BATCH) %>%
  select(GatePath, ref_ctrl_med = ctrl_mean_med)

ctrl_factors_med <- ctrl_means_med %>%
  left_join(ref_ctrl_med, by = "GatePath") %>%
  mutate(
    scale_factor = dplyr::case_when(
      !is.na(ctrl_mean_med) & !is.na(ref_ctrl_med) & ref_ctrl_med != 0 ~ ctrl_mean_med / ref_ctrl_med,
      TRUE ~ 1
    )
  ) %>%
  select(Batch, GatePath, scale_factor)

# ---------------------------
# Step 7a — Also normalize the raw Value for absolute-condition plots
# ---------------------------

cytokine_long <- cytokine_long %>%
  left_join(ctrl_factors_med, by = c("Batch", "GatePath")) %>%
  mutate(
    scale_factor = ifelse(is.na(scale_factor) | scale_factor == 0, 1, scale_factor),
    Value_RawNormalized = Value / scale_factor      # raw signal normalized across batches
  ) %>%
  select(-scale_factor)
# Apply scaling to the working signal (prefer ΔStim; fall back to raw Value)
cytokine_long <- cytokine_long %>%
  left_join(ctrl_factors_med, by = c("Batch","GatePath")) %>%
  mutate(
    scale_factor      = ifelse(is.na(scale_factor) | scale_factor == 0, 1, scale_factor),
    Working           = dplyr::coalesce(Value_BaselineAdj, Value),
    Value_Normalized  = Working / scale_factor
  ) %>%
  select(-scale_factor, -Working)
#If you batch-normalize first, you’d be scaling both MED and Stim by possibly different factors, which could distort within-sample differences.
# ---------------------------
# Step 8 — Transform for modeling (logit on % → proportion)
#Convert a percentage to a proportion and then apply a logit transform (log-odds) so you can use models that assume (approximately) normal residuals.
# ---------------------------
# ---- Step 8a: absolute percentages → logit (for level-based modeling/plots)
eps <- 1e-4
cytokine_long <- cytokine_long %>%
  mutate(
    prop_abs = Value_RawNormalized / 100,                 # needs the raw normalized column
    prop_abs = pmin(pmax(prop_abs, eps), 1 - eps),
    Value_Logit = qlogis(prop_abs)
  ) %>%
  select(-prop_abs)
#This is your batch-corrected absolute percentage of cytokine+ cells,
#converted from % → proportion → logit (log-odds).
#
#---- Step 8b: delta (Stim - MED, batch-corrected) → keep on pp scale
cytokine_long <- cytokine_long %>%
  mutate(
    Delta_for_model = Value_Normalized   # this is your Δ; can be negative; no logit
  )
#induced response (stimulus-specific change)
# ---------------------------
# Step 9 — Save processed cytokine data (RDS + CSV)
# ---------------------------

out_dir <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/saved_R_dat"

# Ensure the folder exists
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save full object (preserves lists and data types)
saveRDS(
  cytokine_long,
  file.path(out_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.rds")
)

# Flatten list-columns for CSV export
cytokine_long_csv <- cytokine_long %>%
  mutate(across(where(is.list), ~ sapply(.x, toString)))

# Save as CSV
readr::write_csv(
  cytokine_long_csv,
  file.path(out_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.csv")
)

