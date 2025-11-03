# install.packages(c("readxl","dplyr","tidyr","stringr","readr")) # if needed
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
#######3 Set directories
#data_dir  <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/Raw_Data"
data_dir  <- "C:/Users/ammas/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/Raw_Data"

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
#xs <- cytokine_long
#cytokine_long <- xs
# ---------------------------
# Step 6 — Baseline adjust to MED within PID × Gate (ΔStim = Stim − MED)
# ---------------------------
med_table <- cytokine_long %>%
  dplyr::filter(Condition == "MED") %>%
  dplyr::group_by(PID, PID_Type, Time_Point, GatePath) %>%
  dplyr::summarise(MED_Value = mean(Value, na.rm = TRUE), .groups = "drop")

cytokine_long <- cytokine_long %>%
  dplyr::left_join(med_table, by = c("PID","PID_Type","Time_Point","GatePath")) %>%
  dplyr::mutate(
    Value_BaselineAdj = dplyr::if_else(!is.na(MED_Value), Value - MED_Value, NA_real_)
  )

# Helpful flags about zeros/near-zeros in raw and MED
min_ctrl_pct <- 0.5      # threshold for safe scaling (in %)
tiny_eps     <- 0.01     # anything below this we treat as "below detection" for flags

cytokine_long <- cytokine_long %>%
  dplyr::mutate(
    ZeroLike_Raw = !is.na(Value)     & Value     < tiny_eps,
    ZeroLike_MED = !is.na(MED_Value) & MED_Value < tiny_eps
  )

# ---------------------------
# Step 7 — Safer batch "normalization" using HC-SP @ MED ONLY
# Do it on the logit scale (additive shift), then invert back to percentage.
# ---------------------------

# ---- Parameters ----
min_ctrl_pct <- 1.0      # ignore controls <1% when computing medians (unchanged)
tiny_eps     <- 1e-4     # clamp for logit
shift_cap    <- 2.0      # cap absolute shift on logit scale (~7.4× odds)

# Helper: logit / invlogit that clamp to avoid Inf
logit_safe   <- function(p) qlogis(pmin(pmax(p, tiny_eps), 1 - tiny_eps))
invlogit_pct <- function(x) plogis(x) * 100

# 7.1 MED controls per Batch×GatePath — compute median on % scale, but only if >=1%
ctrl_meds_med <- cytokine_long %>%
  dplyr::filter(PID_Type == "Batch Control", Condition == "MED") %>%
  dplyr::mutate(Value_use = dplyr::if_else(Value >= min_ctrl_pct, Value, NA_real_)) %>%
  dplyr::group_by(Batch, GatePath) %>%
  dplyr::summarise(ctrl_med_med = median(Value_use, na.rm = TRUE), .groups = "drop")

# Pick reference batch robustly (keep your rule)
REF_BATCH <- if ("Batch 1" %in% unique(ctrl_meds_med$Batch)) "Batch 1" else sort(unique(ctrl_meds_med$Batch))[1]

ref_ctrl_med <- ctrl_meds_med %>%
  dplyr::filter(Batch == REF_BATCH) %>%
  dplyr::select(GatePath, ref_ctrl_med = ctrl_med_med)

# 7.2 Compute logit shift = logit(batch_med) - logit(ref_med), with guards and caps
ctrl_shifts <- ctrl_meds_med %>%
  dplyr::left_join(ref_ctrl_med, by = "GatePath") %>%
  dplyr::mutate(
    ok_batch = !is.na(ctrl_med_med) & ctrl_med_med >= min_ctrl_pct & ctrl_med_med <= 99.9,
    ok_ref   = !is.na(ref_ctrl_med)  & ref_ctrl_med  >= min_ctrl_pct & ref_ctrl_med  <= 99.9,
    p_b      = ctrl_med_med / 100,
    p_r      = ref_ctrl_med / 100,
    logit_b  = dplyr::if_else(ok_batch, logit_safe(p_b), NA_real_),
    logit_r  = dplyr::if_else(ok_ref,   logit_safe(p_r), NA_real_),
    raw_shift = logit_b - logit_r,
    # cap extreme shifts to avoid over-correction
    shift     = dplyr::case_when(
      is.finite(raw_shift) ~ pmax(pmin(raw_shift, shift_cap), -shift_cap),
      TRUE                 ~ NA_real_
    ),
    Shift_Reason = dplyr::case_when(
      ok_batch & ok_ref                           ~ "shifted:MED logit",
      !ok_ref                                     ~ "skip:ref<1% or >99%",
      !ok_batch                                   ~ "skip:batch<1% or >99%",
      TRUE                                        ~ "skip:other"
    )
  ) %>%
  dplyr::select(Batch, GatePath, shift, Shift_Reason)

# 7.3 Apply the logit shift to *absolute* values only; leave ΔStim unscaled
cytokine_long <- cytokine_long %>%
  dplyr::left_join(ctrl_shifts, by = c("Batch", "GatePath")) %>%
  dplyr::mutate(
    # Absolute (% of parent/grandparent) — adjust on logit scale when we have a finite shift
    Value_RawNormalized = dplyr::case_when(
      is.na(Value)                 ~ NA_real_,
      !is.finite(shift)            ~ Value,           # no shift info: keep raw
      Value < min_ctrl_pct         ~ Value,           # don't "correct" noise <1%
      TRUE ~ {
        p <- pmin(pmax(Value / 100, tiny_eps), 1 - tiny_eps)
        invlogit_pct(logit_safe(p) - shift)
      }
    ),
    # ΔStim: prefer Δ (Stim−MED), keep as-is (safer than multiplicative scaling)
    Working_Delta    = dplyr::coalesce(Value_BaselineAdj, Value),
    Value_Normalized = Working_Delta,   # NO scaling for deltas
    Shift_Used       = is.finite(shift) & Value >= min_ctrl_pct
  )

# ---------------------------
# Step 7d — Quick diagnostics: confirm no more huge inflations
# ---------------------------
diag_targets <- c("IL-17a", "CXCR5\\+")  # keep your patterns
diagnostic_df <- cytokine_long %>%
  dplyr::mutate(GatePath_simple = stringr::str_replace_all(GatePath, ",", " ")) %>%
  dplyr::filter(
    stringr::str_detect(GatePath_simple, paste(diag_targets, collapse = "|")),
    !is.na(Value_RawNormalized) | !is.na(Value_Normalized)
  ) %>%
  dplyr::select(
    Batch, PID, Condition, Metric, GatePath,
    Value, MED_Value, Value_BaselineAdj,
    shift, Shift_Used, Shift_Reason,
    Value_RawNormalized, Value_Normalized
  ) %>%
  dplyr::arrange(dplyr::desc(dplyr::coalesce(Value_RawNormalized, Value_Normalized)))

print(head(diagnostic_df, 20))

### 1) See which rows really changed (non-100s)
tol <- 1e-6

changed_examples <- cytokine_long %>%
  filter(!is.na(Value), Value > 1, Value < 99) %>%          # ignore saturated & near-zero
  mutate(delta = Value_RawNormalized - Value) %>%
  filter(is.finite(delta), abs(delta) > 0.05) %>%            # show meaningful changes (>0.05%)
  arrange(desc(abs(delta))) %>%
  select(Batch, PID, Condition, GatePath, Value, Value_RawNormalized, delta, shift, Shift_Used) %>%
  head(25)

print(changed_examples)
summary_shift <- cytokine_long %>%
  mutate(is_mid = !is.na(Value) & Value > 1 & Value < 99) %>%
  summarise(
    n_total       = n(),
    n_mid         = sum(is_mid),
    n_shift_used  = sum(is_mid & Shift_Used, na.rm=TRUE),
    frac_shifted  = n_shift_used / n_mid
  )
print(summary_shift)

# ref_ctrl_med came from your Step 7 build
ref_saturated <- ref_ctrl_med %>%
  mutate(ref_is_sat = ref_ctrl_med >= 99 | ref_ctrl_med < 1) %>%
  arrange(desc(ref_is_sat), desc(ref_ctrl_med)) %>%
  head(30)

print(ref_saturated)

# ---------------------------
# Step 8 — Transforms for modeling
# ---------------------------
# 8a) Absolute levels on logit scale (requires proportions in (0,1); clamp)
eps <- 1e-4
cytokine_long <- cytokine_long %>%
  dplyr::mutate(
    prop_abs = Value_RawNormalized / 100,
    prop_abs = pmin(pmax(prop_abs, eps), 1 - eps),
    Value_Logit = qlogis(prop_abs)
  ) %>%
  dplyr::select(-prop_abs)

# 8b) ΔStim (Stim − MED), batch-normalized — keep on % (can be negative)
cytokine_long <- cytokine_long %>%
  dplyr::mutate(
    Delta_for_model = Value_Normalized
  )

# ---------------------------
# Step 9 — Save processed cytokine data (RDS + CSV)
# ---------------------------
#out_dir <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/saved_R_dat"
out_dir <- "C:/Users/ammas/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/saved_R_dat"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Save RDS (preserves list columns if any)
saveRDS(
  cytokine_long,
  file.path(out_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.rds")
)

# Flatten list-cols for CSV; keep flags for transparency
cytokine_long_csv <- cytokine_long %>%
  dplyr::mutate(dplyr::across(where(is.list), ~ sapply(.x, toString)))

readr::write_csv(
  cytokine_long_csv,
  file.path(out_dir, "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.csv")
)
