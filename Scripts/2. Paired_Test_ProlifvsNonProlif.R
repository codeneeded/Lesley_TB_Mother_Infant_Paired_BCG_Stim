#########################
# Libraries
#########################
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)
library(ggplot2)
library(glue)
library(ragg)
library(ggtext)


#########################
# Paths & I/O
#########################
base_dir <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim"
in_rds   <- file.path(base_dir, "saved_R_dat", "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.rds")
in_csv   <- file.path(base_dir, "saved_R_dat", "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.csv")

# Prefer RDS (preserves types); fallback to CSV if needed
if (file.exists(in_rds)) {
  cytokine_long <- readRDS(in_rds)
} else if (file.exists(in_csv)) {
  cytokine_long <- readr::read_csv(in_csv, show_col_types = FALSE)
} else {
  stop("Could not find preprocessed input at:\n  ", in_rds, "\n  ", in_csv)
}

# Results directories
res_root <- file.path(base_dir, "Results", "Prolif_vs_NonProlif")
if (!dir.exists(res_root)) dir.create(res_root, recursive = TRUE)

#########################
# Analysis settings
#########################
stim_keep <- c("BCG", "DosR", "E6C10", "GAG")
levels(as.factor(cytokine_long$`Time_Point`))
# Identify ONLY T-level columns like T1, T2, ... (avoid Time_Point etc.)
T_cols <- grep("^T[0-9]+$", names(cytokine_long), value = TRUE)

# Build the four groups explicitly
cytokine_an <- cytokine_long %>%
  # Exclude batch controls
  filter(!is.na(PID_Type), PID_Type != "Batch Control") %>%
  # Only non-MED conditions
  filter(Condition %in% stim_keep) %>%
  # We’re comparing absolute levels across batches (no baseline subtraction)
  # so we use the batch-normalized raw value
  mutate(Value_for_test = Value_RawNormalized) %>%
  # Gate consistency: use % of Parent for cytokine gates
  filter(Metric == "Freq. of Parent") %>%
  # Keep only rows with a Proliferation status
  filter(Prolif_Status %in% c("Proliferating","Non-Proliferating")) %>%
  mutate(
    # Build a combined phenotype name across T1..Tn ONLY
    Readout = if (length(T_cols) > 0) {
      pmap_chr(across(all_of(T_cols)), function(...) {
        toks <- c(...)
        toks <- toks[!is.na(toks) & toks != ""]
        if (length(toks) == 0) NA_character_ else str_squish(paste(toks, collapse = " / "))
      })
    } else {
      NA_character_
    }
  ) %>%
  filter(!is.na(Readout)) %>%
  # Define analysis groups
  mutate(
    Time_Point = str_squish(Time_Point),
    Time_Point_std = dplyr::recode(
      Time_Point,
      "12 Wks" = "12wks",
      "44 Wks" = "44wks",
      "Entry"  = "Entry",
      "Dx"     = "Dx",
      .default = Time_Point
    ),
    Group = dplyr::case_when(
      PID_Type == "Mother" & Time_Point_std == "Entry" ~ "Mothers Entry",
      PID_Type == "Mother" & Time_Point_std == "Dx"    ~ "Mothers Dx",
      PID_Type == "Infant" & Time_Point_std == "12wks" ~ "Infants 12wks",
      PID_Type == "Infant" & Time_Point_std == "44wks" ~ "Infants 44wks",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Group)) %>%
  mutate(Compartment = ifelse(is.na(Compartment), "Unknown", Compartment))


# Helper to compute paired Wilcoxon and effect sizes
paired_stats <- function(df) {
  # 1) Clean and standardize Prolif_Status
  df2 <- df %>%
    mutate(
      # normalize spaces
      Prolif_Status = stringr::str_squish(Prolif_Status),
      # normalize unicode dashes to ASCII hyphen
      Prolif_Status = stringr::str_replace_all(Prolif_Status, "[\u2010-\u2015\u2212]", "-"),
      # fix common variant
      Prolif_Status = dplyr::recode(Prolif_Status,
                                    "Non Proliferating" = "Non-Proliferating",
                                    .default = Prolif_Status),
      # force the two expected levels so pivot_wider creates both cols
      Prolif_Status = factor(Prolif_Status,
                             levels = c("Proliferating", "Non-Proliferating"))
    )
  
  # 2) Wide table with both columns guaranteed
  wide <- df2 %>%
    dplyr::select(PID, Prolif_Status, Value_for_test) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      names_from  = Prolif_Status,
      values_from = Value_for_test,
      values_fill = NA_real_
    )
  
  # 3) Keep complete pairs
  if (!("Proliferating" %in% names(wide))) wide$Proliferating <- NA_real_
  if (!("Non-Proliferating" %in% names(wide))) wide$`Non-Proliferating` <- NA_real_
  
  wide <- wide %>%
    dplyr::filter(!is.na(Proliferating), !is.na(`Non-Proliferating`))
  
  n_pairs <- nrow(wide)
  if (n_pairs < 1) {
    return(dplyr::tibble(
      n_pairs = 0,
      median_Prolif = NA_real_,
      median_NonProlif = NA_real_,
      median_diff_pp = NA_real_,
      wilcox_stat = NA_real_,
      p_value = NA_real_,
      effect_rbc = NA_real_
    ))
  }
  
  # 4) Paired Wilcoxon
  w <- suppressWarnings(stats::wilcox.test(
    x = wide$Proliferating,
    y = wide$`Non-Proliferating`,
    paired = TRUE,
    exact = FALSE,
    conf.int = FALSE
  ))
  V <- unname(w$statistic)
  
  # Rank-biserial correlation for paired signed-rank
  r_rb <- if (!is.null(V) && n_pairs > 0) {
    1 - (2 * as.numeric(V)) / (n_pairs * (n_pairs + 1))
  } else {
    NA_real_
  }
  
  dplyr::tibble(
    n_pairs = n_pairs,
    median_Prolif = stats::median(wide$Proliferating, na.rm = TRUE),
    median_NonProlif = stats::median(wide$`Non-Proliferating`, na.rm = TRUE),
    median_diff_pp = stats::median(wide$Proliferating - wide$`Non-Proliferating`, na.rm = TRUE),
    wilcox_stat = as.numeric(V),
    p_value = w$p.value,
    effect_rbc = r_rb
  )
}

#########################
# Run tests
#########################

# Compute stats for every combo
stats_all <- cytokine_an %>%
  group_by(Readout, Group, Condition, Compartment) %>%
  group_modify(~ paired_stats(.x)) %>%
  ungroup()

# Adjust p-values per Readout (FDR across Group×Condition×Compartment within that cytokine)
stats_all <- stats_all %>%
  group_by(Readout) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Save combined summary
write_csv(stats_all, file.path(res_root, "_ALL_readouts_summary.csv"))

#########################
# Save per-readout folders & files
#########################

### Paths
# Build the universe of readouts from the data (single + combinatorial)
readouts <- cytokine_an %>%
  dplyr::distinct(Readout) %>%
  dplyr::pull(Readout) %>%
  sort()
#
# Helper to make safe folder/file names from readout labels
sanitize_for_path <- function(x) {
  x |>
    stringr::str_replace_all("[/\\\\|?:*<>%]", "_") |>  # remove path-breaking chars
    stringr::str_replace_all("[\\s]+", " ") |>         # collapse whitespace
    stringr::str_trim()
}

# Write per-readout outputs
purrr::walk(readouts, function(rd) {
  safe_rd <- sanitize_for_path(rd)
  subdir  <- file.path(res_root, safe_rd)
  if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)
  
  # Per-readout stats
  stats_rd <- stats_all %>% dplyr::filter(Readout == rd)
  readr::write_csv(stats_rd, file.path(subdir, paste0(safe_rd, "_paired_stats.csv")))
  
  # PID-level paired differences for this readout
  diffs_rd <- pid_diffs %>% dplyr::filter(Readout == rd)
  readr::write_csv(diffs_rd, file.path(subdir, paste0(safe_rd, "_paired_differences.csv")))
})

# ========= PLOTTING =========
library(ggplot2)
library(glue)

# Helper again (if not already in your script)
sanitize_for_path <- function(x) {
  x |>
    stringr::str_replace_all("[/\\\\|?:*<>%]", "_") |>
    stringr::str_replace_all("[\\s]+", " ") |>
    stringr::str_trim()
}

# Use the same readouts vector we built from the data
if (!exists("readouts")) {
  readouts <- cytokine_an %>% dplyr::distinct(Readout) %>% dplyr::pull(Readout) %>% sort()
}

# Ensure status order for x-axis
cytokine_an <- cytokine_an %>%
  dplyr::mutate(Prolif_Status = factor(Prolif_Status, levels = c("Non-Proliferating","Proliferating")))

# Plot loop
purrr::walk(readouts, function(rd) {
  safe_rd <- sanitize_for_path(rd)
  subdir  <- file.path(res_root, safe_rd)
  if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)
  
  data_rd  <- cytokine_an %>% dplyr::filter(Readout == rd)
  combos   <- data_rd %>% dplyr::distinct(Group, Condition, Compartment)
  
  purrr::pwalk(
    combos,
    function(Group, Condition, Compartment) {
      df <- data_rd %>%
        dplyr::filter(Group == !!Group, Condition == !!Condition, Compartment == !!Compartment) %>%
        dplyr::select(PID, Prolif_Status, Value_for_test) %>%
        dplyr::distinct() %>%
        # --- normalize status & force both levels ---
        dplyr::mutate(
          Prolif_Status = stringr::str_squish(Prolif_Status),
          Prolif_Status = stringr::str_replace_all(Prolif_Status, "[\u2010-\u2015\u2212]", "-"),
          Prolif_Status = dplyr::recode(Prolif_Status,
                                        "Non Proliferating" = "Non-Proliferating",
                                        .default = Prolif_Status),
          Prolif_Status = factor(Prolif_Status, levels = c("Non-Proliferating","Proliferating"))
        )
      
      # Wide table with both columns guaranteed
      wide <- df %>%
        tidyr::pivot_wider(
          names_from  = Prolif_Status,
          values_from = Value_for_test,
          values_fill = NA_real_
        )
      
      # If any status column is missing, create it (rare but safe)
      if (!("Proliferating" %in% names(wide)))         wide$Proliferating <- NA_real_
      if (!("Non-Proliferating" %in% names(wide))) wide$`Non-Proliferating` <- NA_real_
      
      # keep complete pairs
      wide <- wide %>% dplyr::filter(!is.na(Proliferating), !is.na(`Non-Proliferating`))
      n_pairs <- nrow(wide)
      if (n_pairs < 1) return(invisible(NULL))
      
      # stats for subtitle (pull from computed stats_all if available)
      stat_row <- stats_all %>%
        dplyr::filter(Readout == !!rd, Group == !!Group, Condition == !!Condition, Compartment == !!Compartment) %>%
        dplyr::slice(1)
      
      if (nrow(stat_row) == 0) {
        w <- suppressWarnings(stats::wilcox.test(wide$Proliferating, wide$`Non-Proliferating`, paired = TRUE, exact = FALSE))
        median_diff <- stats::median(wide$Proliferating - wide$`Non-Proliferating`, na.rm = TRUE)
        pval <- w$p.value
        padj <- NA_real_
      } else {
        median_diff <- stat_row$median_diff_pp
        pval <- stat_row$p_value
        padj <- stat_row$p_adj
      }
      
      subtitle_txt <- glue::glue(
        "n_pairs = {n_pairs} | median Δ(P − NP) = {round(median_diff, 3)} | p = {signif(pval, 3)}{ifelse(is.na(padj),'', paste0(' | FDR = ', signif(padj,3)))}"
      )
      
      # long format (with normalized status) for plotting
      df_plot <- df
      
      # dynamic width
      w_in <- max(4.5, min(7, 4.5 + 0.05 * n_pairs))
      h_in <- 4.5
      
      p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = Prolif_Status, y = Value_for_test, group = PID)) +
        ggplot2::geom_line(alpha = 0.35) +
        ggplot2::geom_point(ggplot2::aes(color = Prolif_Status), size = 2) +
        ggplot2::stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 0, alpha = 0.7) +
        ggplot2::scale_x_discrete(drop = FALSE) +
        ggplot2::labs(
          title = glue::glue("{rd} — {Group} — {Condition} — {Compartment}"),
          subtitle = subtitle_txt,
          x = NULL,
          y = "% of Parent (batch-normalized)"
        ) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          legend.position = "none",
          plot.title = ggplot2::element_text(face = "bold"),
          panel.grid.major.x = ggplot2::element_blank()
        )
      
      fname <- file.path(
        subdir,
        paste0(
          sanitize_for_path(glue::glue("{rd}__{Group}__{Condition}__{Compartment}")),
          ".png"
        )
      )
      ggplot2::ggsave(filename = fname, plot = p, width = w_in, height = h_in, dpi = 300, bg = "white")
    }
  )
})

#Within Proliferating only compare:
#1.   Infant 12 vs 44wks (timepoint comparison, paired, within each stim (BCG and doser for infants))
#2 Infant - within each timepoint and condition- Maternal IGRA Pos vs Neg (Exclude ATB)
#3 Infants - within each timepoint and condition- Infant IGRA Pos vs Neg (this is NEG biased- can we control loopk at data)
#4. Mothers - within each conditon - IMaternal IGRA Pos vs NEg vs ATB
#column ids # PID Type  = tells you if its Batch Control Maternal or Infant, Time_Point has timepoint for Infants (mothers which is 12 or 44wks), mothers have no timepoint
#Condition has the stim (MED= media or no stim, other conditions are BCG, GAG, DoseR etc...) Maternal IGRA and Infant IGRA status is the IGRA status for infants and Maternal for mother
#for infants, the Maternal IGRA means the IGRA status of the infants mother and vice versa for the mother

# ================================
# PROLIFERATING-ONLY COMPARISONS
# ================================

# Output directory
# ----------------
comp_dir <- file.path(base_dir, "Results", "Prolif_Only", "Comparisons")
if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE)

# ---------------------------------------------------------
# Build proliferating-only analysis table with full context
# ---------------------------------------------------------
# - Keeps Metric ("Freq. of Parent"/"Freq. of Grandparent")
# - DenominatorGate (explicit ref population)
# - AnchorGate (gate immediately upstream of Proliferating branch)
#   -> If ParentGate == "Proliferating", AnchorGate = GrandparentGate (e.g., CD4/CD8)
#   -> else AnchorGate = ParentGate
# - Compartment is retained (falls back to AnchorGate if missing)
# - Readout is built from T1..Tn only (post-Proliferating path)
# - Value_for_test: absolute % (batch-normalized), clamped to [0, 100]
# ---------------------------------------------------------
# Build proliferating-only analysis table with full context
# ---------------------------------------------------------

prolif <- cytokine_long %>%
  filter(!is.na(PID_Type), PID_Type != "Batch Control") %>%
  mutate(Condition = str_squish(Condition)) %>%
  filter(Prolif_Status == "Proliferating",
         Metric %in% c("Freq. of Parent", "Freq. of Grandparent")) %>%
  # --- reparse path to get the segment immediately BEFORE "Proliferating" ---
  mutate(
    PathCore = str_replace(GatePath, "\\s*\\|\\s*Freq\\..*$", ""),
    PathCore = str_replace_all(PathCore, ",", " "),
    PathCore = str_squish(PathCore),
    LineageGate = map_chr(
      str_split(PathCore, "/"),
      function(segs) {
        segs <- segs[segs != ""]
        idx  <- which(segs == "Proliferating")
        if (length(idx) && idx[1] > 1) segs[idx[1] - 1] else NA_character_
      }
    ),
    # normalize common lineage labels
    LineageGate = case_when(
      LineageGate %in% c("gd TCR", "gdTCR") ~ "gd TCR",
      TRUE ~ LineageGate
    ),
    Time_Point  = str_squish(Time_Point),
    Time_Point_std = recode(
      Time_Point,
      "12 Wks" = "12wks",
      "44 Wks" = "44wks",
      "Entry"  = "Entry",
      "Dx"     = "Dx",
      .default = Time_Point
    ),
    # Keep Compartment; if missing, fall back to LineageGate
    Compartment = coalesce(Compartment, LineageGate),
    # Build Readout from T1..Tn (post-Proliferating)
    Readout = {
      T_cols <- grep("^T[0-9]+$", names(.), value = TRUE)
      if (length(T_cols) > 0) {
        pmap_chr(across(all_of(T_cols)), function(...) {
          toks <- c(...); toks <- toks[!is.na(toks) & toks != ""]
          if (length(toks) == 0) NA_character_ else str_squish(paste(toks, collapse = " / "))
        })
      } else NA_character_
    },
    Value_for_test = pmin(pmax(Value_RawNormalized, 0), 100),
    DenominatorLabel = case_when(
      Metric == "Freq. of Parent"      ~ paste0("Freq. of Parent: ", DenominatorGate),
      Metric == "Freq. of Grandparent" ~ paste0("Freq. of Grandparent: ", DenominatorGate),
      TRUE ~ Metric
    ),
    Group = case_when(
      PID_Type == "Mother" & Time_Point_std == "Entry" ~ "Mothers Entry",
      PID_Type == "Mother" & Time_Point_std == "Dx"    ~ "Mothers Dx",
      PID_Type == "Infant" & Time_Point_std == "12wks" ~ "Infants 12wks",
      PID_Type == "Infant" & Time_Point_std == "44wks" ~ "Infants 44wks",
      PID_Type == "Mother" & is.na(Time_Point_std)     ~ "Mothers",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Readout)) %>%
  rename(
    Maternal_IGRA = `Maternal IGRA`,
    Infant_IGRA   = `Infant IGRA`
  ) %>%
  select(
    PID, PID_Type, Group, Time_Point_std, Condition,
    Metric, DenominatorGate, DenominatorLabel,
    ParentGate, GrandparentGate,
    LineageGate, Compartment,   # <- NEW: lineage just before Proliferating
    Readout, Maternal_IGRA, Infant_IGRA, Value_for_test
  )

# -------------------------
# Effect size helper functions
# -------------------------
paired_rbc <- function(x, y) {
  x <- x[is.finite(x) & is.finite(y)]
  y <- y[is.finite(x) & is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  w <- suppressWarnings(stats::wilcox.test(x, y, paired = TRUE, exact = FALSE))
  V <- unname(w$statistic); n <- length(x)
  if (is.null(V) || n < 1) return(NA_real_)
  1 - (2 * as.numeric(V)) / (n * (n + 1))
}

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  diffs <- outer(x, y, "-")
  (sum(diffs > 0) - sum(diffs < 0)) / (length(x) * length(y))
}

# ==============================================
# 1) INFANTS: Paired 12wks vs 44wks within stimulus
# ==============================================
stim_pair <- c("BCG","DosR")

comp1 <- prolif %>%
  filter(PID_Type == "Infant", Condition %in% c("BCG","DosR")) %>%   # <- only DosR (and BCG)
  filter(Time_Point_std %in% c("12wks","44wks")) %>%
  group_by(Readout, Condition, Metric, DenominatorGate, DenominatorLabel, LineageGate) %>%
  group_modify(~{
    wide <- .x %>%
      select(PID, Time_Point_std, Value_for_test) %>%
      group_by(PID, Time_Point_std) %>%
      summarise(Value = median(Value_for_test, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Time_Point_std, values_from = Value) %>%
      filter(is.finite(`12wks`), is.finite(`44wks`))
    n_pairs <- nrow(wide)
    if (n_pairs < 1) {
      tibble(n_pairs=0, median_12w=NA_real_, median_44w=NA_real_,
             median_diff_pp=NA_real_, p_value=NA_real_, effect_rbc=NA_real_)
    } else {
      w <- suppressWarnings(wilcox.test(wide$`12wks`, wide$`44wks`, paired = TRUE, exact = FALSE))
      V <- unname(w$statistic)
      tibble(
        n_pairs        = n_pairs,
        median_12w     = median(wide$`12wks`, na.rm = TRUE),
        median_44w     = median(wide$`44wks`, na.rm = TRUE),
        median_diff_pp = median(wide$`44wks` - wide$`12wks`, na.rm = TRUE),
        p_value        = w$p.value,
        effect_rbc     = if (is.null(V)) NA_real_ else 1 - (2*as.numeric(V))/(n_pairs*(n_pairs+1))
      )
    }
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(Condition, Metric, DenominatorGate, LineageGate, desc(abs(median_diff_pp)))

write_csv(comp1, file.path(comp_dir, "Comp1_Infants_Paired_12w_vs_44w_byStim.csv"))


# ===========================================================
# 2) INFANTS: Maternal IGRA Pos vs Neg (by TP × stim)
# ===========================================================
comp2 <- prolif %>%
  filter(PID_Type == "Infant",
         Time_Point_std %in% c("12wks","44wks"),
         !is.na(Maternal_IGRA),
         Maternal_IGRA %in% c("Pos","Neg")) %>%
  group_by(Readout, Time_Point_std, Condition, Metric, DenominatorGate, DenominatorLabel, LineageGate) %>%
  group_modify(~{
    dat  <- .x %>% select(Maternal_IGRA, Value_for_test)
    gpos <- dat %>% filter(Maternal_IGRA == "Pos") %>% pull(Value_for_test)
    gneg <- dat %>% filter(Maternal_IGRA == "Neg") %>% pull(Value_for_test)
    if (!length(gpos) || !length(gneg)) {
      tibble(n_pos=length(gpos), n_neg=length(gneg),
             median_pos=NA_real_, median_neg=NA_real_,
             median_diff_pp=NA_real_, p_value=NA_real_, cliffs_delta=NA_real_)
    } else {
      w <- suppressWarnings(wilcox.test(gpos, gneg, paired = FALSE, exact = FALSE))
      tibble(
        n_pos = length(gpos), n_neg = length(gneg),
        median_pos = median(gpos, na.rm = TRUE),
        median_neg = median(gneg, na.rm = TRUE),
        median_diff_pp = median(gpos, na.rm = TRUE) - median(gneg, na.rm = TRUE),
        p_value = w$p.value,
        cliffs_delta = {
          diffs <- outer(gpos, gneg, "-")
          (sum(diffs > 0) - sum(diffs < 0)) / (length(gpos) * length(gneg))
        }
      )
    }
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(Time_Point_std, Condition, Metric, DenominatorGate, LineageGate, desc(abs(median_diff_pp)))


write_csv(comp2, file.path(comp_dir, "Comp2_Infants_MaternalIGRA_Pos_vs_Neg_byTP_byStim.csv"))

# ===========================================================
# 3) INFANTS: Infant IGRA Pos vs Neg (by TP × stim)
# ===========================================================
comp3 <- prolif %>%
  filter(PID_Type == "Infant",
         Time_Point_std %in% c("12wks","44wks"),
         !is.na(Infant_IGRA),
         Infant_IGRA %in% c("Pos","Neg")) %>%
  group_by(Readout, Time_Point_std, Condition, Metric, DenominatorGate, DenominatorLabel, LineageGate) %>%
  group_modify(~{
    dat  <- .x %>% select(Infant_IGRA, Value_for_test)
    gpos <- dat %>% filter(Infant_IGRA == "Pos") %>% pull(Value_for_test)
    gneg <- dat %>% filter(Infant_IGRA == "Neg") %>% pull(Value_for_test)
    if (!length(gpos) || !length(gneg)) {
      tibble(n_pos=length(gpos), n_neg=length(gneg),
             median_pos=NA_real_, median_neg=NA_real_,
             median_diff_pp=NA_real_, p_value=NA_real_, cliffs_delta=NA_real_)
    } else {
      w <- suppressWarnings(wilcox.test(gpos, gneg, paired = FALSE, exact = FALSE))
      tibble(
        n_pos = length(gpos), n_neg = length(gneg),
        median_pos = median(gpos, na.rm = TRUE),
        median_neg = median(gneg, na.rm = TRUE),
        median_diff_pp = median(gpos, na.rm = TRUE) - median(gneg, na.rm = TRUE),
        p_value = w$p.value,
        cliffs_delta = {
          diffs <- outer(gpos, gneg, "-")
          (sum(diffs > 0) - sum(diffs < 0)) / (length(gpos) * length(gneg))
        }
      )
    }
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    note_imbalance = case_when(
      n_pos < 5 & n_neg >= 10 ~ "POS small",
      n_neg < 5 & n_pos >= 10 ~ "NEG small",
      n_pos < 5 & n_neg < 5  ~ "both small",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(Time_Point_std, Condition, Metric, DenominatorGate, LineageGate, desc(abs(median_diff_pp)))


write_csv(comp3, file.path(comp_dir, "Comp3_Infants_InfantIGRA_Pos_vs_Neg_byTP_byStim.csv"))

# ===========================================================
# 4) MOTHERS: IGRA Pos vs Neg vs ATB (per stimulus)
# ===========================================================
comp4_overall <- prolif %>%
  filter(PID_Type == "Mother",
         !is.na(Maternal_IGRA),
         Maternal_IGRA %in% c("Pos","Neg","ATB")) %>%
  group_by(Readout, Condition, Metric, DenominatorGate, DenominatorLabel, LineageGate) %>%
  group_modify(~{
    dat <- .x %>% select(Maternal_IGRA, Value_for_test) %>% drop_na()
    if (nrow(dat) < 3 || length(unique(dat$Maternal_IGRA)) < 2) {
      tibble(n = nrow(dat), k = length(unique(dat$Maternal_IGRA)),
             kw_stat = NA_real_, kw_p = NA_real_)
    } else {
      kw <- suppressWarnings(kruskal.test(Value_for_test ~ Maternal_IGRA, data = dat))
      tibble(n = nrow(dat), k = length(unique(dat$Maternal_IGRA)),
             kw_stat = unname(kw$statistic), kw_p = kw$p.value)
    }
  }) %>%
  ungroup() %>%
  mutate(kw_p_adj = p.adjust(kw_p, method = "BH")) %>%
  arrange(Condition, Metric, DenominatorGate, LineageGate, kw_p)

comp4_pairwise <- prolif %>%
  filter(PID_Type == "Mother",
         !is.na(Maternal_IGRA),
         Maternal_IGRA %in% c("Pos","Neg","ATB")) %>%
  group_by(Readout, Condition, Metric, DenominatorGate, DenominatorLabel, LineageGate) %>%
  group_modify(~{
    dat <- .x %>% select(Maternal_IGRA, Value_for_test) %>% drop_na()
    if (length(unique(dat$Maternal_IGRA)) < 2) {
      tibble(comp = character(), p = numeric())
    } else {
      pw <- suppressWarnings(pairwise.wilcox.test(
        dat$Value_for_test, dat$Maternal_IGRA, p.adjust.method = "BH", exact = FALSE
      ))
      mat <- pw$p.value
      if (is.null(mat)) {
        tibble(comp = character(), p = numeric())
      } else {
        as_tibble(as.data.frame(as.table(mat)), .name_repair = "minimal") %>%
          rename(group1 = Var1, group2 = Var2, p = Freq) %>%
          filter(!is.na(p)) %>%
          mutate(comp = paste(group1, "vs", group2)) %>%
          select(comp, p)
      }
    }
  }) %>%
  ungroup() %>%
  arrange(Condition, Metric, DenominatorGate, LineageGate, comp)


write_csv(comp4_overall, file.path(comp_dir, "Comp4_Mothers_KruskalWallis_byStim.csv"))

write_csv(comp4_pairwise, file.path(comp_dir, "Comp4_Mothers_PairwiseWilcoxon_byStim.csv"))

# ================================
# PLOTTING SUITE — Proliferating Only
# ================================

# ----------------
# Output directories
# ----------------
plots_root <- file.path(base_dir, "Results", "Prolif_Only", "Plots")
dirs <- list(
  comp1_slopes   = file.path(plots_root, "Comp1_Slopes"),
  comp1_heatmap  = file.path(plots_root, "Comp1_Summary"),
  comp2_violin   = file.path(plots_root, "Comp2_Violin"),
  comp2_heatmap  = file.path(plots_root, "Comp2_Summary"),
  comp3_violin   = file.path(plots_root, "Comp3_Violin"),
  comp3_heatmap  = file.path(plots_root, "Comp3_Summary"),
  comp4_box      = file.path(plots_root, "Comp4_Box"),
  comp4_heatmap  = file.path(plots_root, "Comp4_Summary")
)
invisible(lapply(dirs, function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)))

# ----------------
# Helpers
# ----------------
sanitize_for_path <- function(x) {
  x |>
    stringr::str_replace_all("[/\\\\|?:*<>%]", "_") |>
    stringr::str_replace_all("[\\s]+", " ") |>
    stringr::str_trim()
}
wrap_str <- function(x, width = 70) stringr::str_wrap(x, width = width)

sig_stars <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*", ""))))
}

# Label that makes lineage + denominator explicit
mk_readout_label <- function(lineage, readout, denom_label) {
  lineage <- ifelse(is.na(lineage) | lineage == "", "Lineage?", lineage)
  rl <- ifelse(is.na(readout) | readout == "", "Readout?", readout)
  paste0(lineage, " :: ", rl, " [", denom_label, "]")
}
# -------- anti-clipping helpers --------
wrap_str <- function(x, width = 80) stringr::str_wrap(x, width = width)

count_lines <- function(x) {
  if (is.null(x) || is.na(x)) return(1L)
  length(unlist(strsplit(x, "\n", fixed = TRUE)))
}

# Given a (title, subtitle), return wrapped strings + suggested height bump
wrap_title_sub <- function(title, subtitle, title_width = 80, sub_width = 100,
                           base_height = 4.6, extra_per_line = 0.28) {
  tw <- wrap_str(title, title_width)
  sw <- wrap_str(subtitle, sub_width)
  n_title <- max(1L, count_lines(tw))
  n_sub   <- max(1L, count_lines(sw))
  # add height for extra lines beyond 1
  h <- base_height + extra_per_line * ((n_title - 1L) + (n_sub - 1L))
  list(title = tw, subtitle = sw, height = h)
}

# Use in theme() to keep generous margins + prevent cuts
anti_clip_theme <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) %+replace%
    ggplot2::theme(
      legend.position    = "none",
      plot.title.position= "plot",
      plot.title         = ggplot2::element_text(face = "bold", margin = ggplot2::margin(b = 8)),
      plot.subtitle      = ggplot2::element_text(margin = ggplot2::margin(b = 10)),
      plot.caption       = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(t = 20, r = 28, b = 20, l = 20)
    )
}

sig_stars <- function(p){
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*", ""))))
}

mk_readout_label <- function(lineage, readout, denom_label) {
  lineage <- ifelse(is.na(lineage) | lineage == "", "Lineage?", lineage)
  rl <- ifelse(is.na(readout) | readout == "", "Readout?", readout)
  paste0(lineage, " :: ", rl, " [", denom_label, "]")
}

# ======================================================================
# COMP 1 (Paired): Infant 12wks vs 44wks within stimulus (BCG, DosR only)
# ======================================================================

# 1A) SLOPE PLOTS — per (Readout × LineageGate × Condition)
# Uses the same pairing as comp1: we reconstruct pairs from `prolif` to draw lines.
plot_comp1_slope <- function(df, condition, readout, lineage, outdir) {
  dat <- df %>%
    dplyr::filter(PID_Type == "Infant",
                  Condition == condition,
                  Time_Point_std %in% c("12wks","44wks"),
                  Readout == readout,
                  LineageGate == lineage) %>%
    dplyr::select(PID, Time_Point_std, Value_for_test, DenominatorLabel) %>%
    dplyr::group_by(PID, Time_Point_std) %>%
    dplyr::summarise(Value = median(Value_for_test, na.rm = TRUE),
                     DenominatorLabel = dplyr::first(DenominatorLabel), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Time_Point_std, values_from = Value) %>%
    dplyr::filter(is.finite(`12wks`), is.finite(`44wks`)) %>%
    dplyr::mutate(diff = `44wks` - `12wks`)
  n_pairs <- nrow(dat)
  if (n_pairs < 1) return(invisible(NULL))
  
  w <- suppressWarnings(stats::wilcox.test(dat$`12wks`, dat$`44wks`, paired = TRUE, exact = FALSE))
  med12 <- median(dat$`12wks`, na.rm = TRUE)
  med44 <- median(dat$`44wks`, na.rm = TRUE)
  meddiff <- median(dat$diff, na.rm = TRUE)
  pval <- w$p.value
  
  denom_lbl <- unique(dat$DenominatorLabel)[1]
  title_txt <- glue::glue("{mk_readout_label(lineage, readout, denom_lbl)} — {condition}")
  subt_txt  <- glue::glue("n={n_pairs} | median 12w={round(med12,2)}, 44w={round(med44,2)} | median Δ(44−12)={round(meddiff,2)} | p={signif(pval,3)} {sig_stars(pval)}")
  
  long <- dat |>
    dplyr::select(PID, `12wks`, `44wks`) |>
    tidyr::pivot_longer(cols = c(`12wks`,`44wks`), names_to = "Time", values_to = "Value")
  long$Time <- factor(long$Time, levels = c("12wks","44wks"))
  
  # wrap + dynamic height and width
  wraps <- wrap_title_sub(title_txt, subt_txt, title_width = 90, sub_width = 110,
                          base_height = 4.8, extra_per_line = 0.28)
  w_in <- max(5.2, min(8.6, 5.2 + 0.05*n_pairs))
  h_in <- wraps$height
  
  p <- ggplot2::ggplot(long, ggplot2::aes(x = Time, y = Value, group = PID)) +
    ggplot2::geom_line(alpha = 0.38, linewidth = 0.55) +
    ggplot2::geom_point(ggplot2::aes(color = Time), size = 2, position = ggplot2::position_identity()) +
    ggplot2::stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 0, alpha = 0.92) +
    ggplot2::coord_cartesian(ylim = c(0, 100), clip = "off") +   # <- anti-clipping
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::labs(
      title    = wraps$title,
      subtitle = wraps$subtitle,
      x = NULL, y = "% of Parent (batch-normalized)"
    ) +
    anti_clip_theme()
  
  fname <- file.path(outdir, paste0(
    sanitize_for_path(glue::glue("{lineage}__{readout}__{condition}__slope")), ".png"))
  ggplot2::ggsave(fname, p, width = w_in, height = h_in, units = "in",
                  dpi = 300, bg = "white")
}

# Generate slope plots for top effects (or all). Here: top 25 by |median_diff_pp|
comp1_top <- comp1 %>%
  filter(Condition %in% c("BCG","DosR")) %>%
  mutate(Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel)) %>%
  arrange(Condition, desc(abs(median_diff_pp))) %>%
  group_by(Condition) %>%
  slice_head(n = 25) %>%
  ungroup()

# Draw slopes
invisible(pmap(
  list(comp1_top$Condition, comp1_top$Readout, comp1_top$LineageGate),
  ~ plot_comp1_slope(prolif, ..1, ..2, ..3, dirs$comp1_slopes)
))



# ======================================================
# COMP 2 (Unpaired): Maternal IGRA Pos vs Neg (Infants)
# ======================================================
# 2A) VIOLIN/BOX/JITTER — per (TP × Stim × Readout × Lineage)
plot_comp2_violin <- function(df, tp, stim, readout, lineage, outdir) {
  dat <- df %>%
    dplyr::filter(PID_Type == "Infant",
                  Time_Point_std == tp,
                  Condition == stim,
                  Readout == readout,
                  LineageGate == lineage,
                  !is.na(Maternal_IGRA),
                  Maternal_IGRA %in% c("Pos","Neg")) %>%
    dplyr::select(Maternal_IGRA, Value_for_test, DenominatorLabel)
  
  # need at least one value in each group
  if (nrow(dat) < 2 || length(unique(dat$Maternal_IGRA)) < 2) return(invisible(NULL))
  
  # --- vector interface (no formula), so we can set paired=FALSE explicitly ---
  g_pos <- dat$Value_for_test[dat$Maternal_IGRA == "Pos"]
  g_neg <- dat$Value_for_test[dat$Maternal_IGRA == "Neg"]
  g_pos <- g_pos[is.finite(g_pos)]
  g_neg <- g_neg[is.finite(g_neg)]
  if (length(g_pos) < 1L || length(g_neg) < 1L) return(invisible(NULL))
  
  w <- suppressWarnings(stats::wilcox.test(g_pos, g_neg, paired = FALSE, exact = FALSE))
  pval  <- w$p.value
  gpos  <- length(g_pos)
  gneg  <- length(g_neg)
  
  denom_lbl <- unique(dat$DenominatorLabel)[1]
  title_txt <- glue::glue("{mk_readout_label(lineage, readout, denom_lbl)} — {stim} @ {tp}")
  subt_txt  <- glue::glue("Pos n={gpos}, Neg n={gneg} | p={signif(pval,3)} {sig_stars(pval)}")
  
  wraps <- wrap_title_sub(title_txt, subt_txt, title_width = 90, sub_width = 110,
                          base_height = 4.8, extra_per_line = 0.28)
  
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = Maternal_IGRA, y = Value_for_test, fill = Maternal_IGRA)) +
    ggplot2::geom_violin(width = 0.9, alpha = 0.25, color = NA, trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.46, outlier.shape = NA, alpha = 0.95) +
    ggplot2::geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.7) +
    ggplot2::scale_y_continuous(limits = c(0,100), expand = ggplot2::expansion(mult = c(0.02, 0.06))) +
    ggplot2::coord_cartesian(clip = "off") +  # anti-clipping
    ggplot2::labs(title = wraps$title, subtitle = wraps$subtitle,
                  x = NULL, y = "% of Parent (batch-normalized)") +
    anti_clip_theme()
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  fname <- file.path(outdir, paste0(
    sanitize_for_path(glue::glue("{lineage}__{readout}__{stim}__{tp}__MaternalIGRA_violin")),".png"))
  ggplot2::ggsave(fname, p, width = 6.4, height = wraps$height, units = "in",
                  dpi = 300, bg = "white")
}


comp2_top <- comp2 %>%
  mutate(Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel)) %>%
  arrange(Time_Point_std, Condition, desc(abs(median_diff_pp))) %>%
  group_by(Time_Point_std, Condition) %>%
  slice_head(n = 25) %>%
  ungroup()

invisible(pmap(
  list(comp2_top$Time_Point_std, comp2_top$Condition, comp2_top$Readout, comp2_top$LineageGate),
  ~ plot_comp2_violin(prolif, ..1, ..2, ..3, ..4, dirs$comp2_violin)
))



# ======================================================
# COMP 3 (Unpaired): Infant IGRA Pos vs Neg (Infants)
# ======================================================

# 3A) VIOLIN/BOX/JITTER — per (TP × Stim × Readout × Lineage)
plot_comp3_violin <- function(df, tp, stim, readout, lineage, outdir) {
  dat <- df %>%
    dplyr::filter(PID_Type == "Infant",
                  Time_Point_std == tp,
                  Condition == stim,
                  Readout == readout,
                  LineageGate == lineage,
                  !is.na(Infant_IGRA),
                  Infant_IGRA %in% c("Pos","Neg")) %>%
    dplyr::select(Infant_IGRA, Value_for_test, DenominatorLabel)
  
  # need at least one value in each group
  if (nrow(dat) < 2 || length(unique(dat$Infant_IGRA)) < 2) return(invisible(NULL))
  
  # --- vector interface ---
  g_pos <- dat$Value_for_test[dat$Infant_IGRA == "Pos"]
  g_neg <- dat$Value_for_test[dat$Infant_IGRA == "Neg"]
  g_pos <- g_pos[is.finite(g_pos)]
  g_neg <- g_neg[is.finite(g_neg)]
  if (length(g_pos) < 1L || length(g_neg) < 1L) return(invisible(NULL))
  
  w <- suppressWarnings(stats::wilcox.test(g_pos, g_neg, paired = FALSE, exact = FALSE))
  pval  <- w$p.value
  gpos  <- length(g_pos)
  gneg  <- length(g_neg)
  
  denom_lbl <- unique(dat$DenominatorLabel)[1]
  title_txt <- glue::glue("{mk_readout_label(lineage, readout, denom_lbl)} — {stim} @ {tp}")
  subt_txt  <- glue::glue("Pos n={gpos}, Neg n={gneg} | p={signif(pval,3)} {sig_stars(pval)}")
  
  wraps <- wrap_title_sub(title_txt, subt_txt, title_width = 90, sub_width = 110,
                          base_height = 4.8, extra_per_line = 0.28)
  
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = Infant_IGRA, y = Value_for_test, fill = Infant_IGRA)) +
    ggplot2::geom_violin(width = 0.9, alpha = 0.25, color = NA, trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.46, outlier.shape = NA, alpha = 0.95) +
    ggplot2::geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.7) +
    ggplot2::scale_y_continuous(limits = c(0,100), expand = ggplot2::expansion(mult = c(0.02, 0.06))) +
    ggplot2::coord_cartesian(clip = "off") +  # anti-clipping
    ggplot2::labs(title = wraps$title, subtitle = wraps$subtitle,
                  x = NULL, y = "% of Parent (batch-normalized)") +
    anti_clip_theme()
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  fname <- file.path(outdir, paste0(
    sanitize_for_path(glue::glue("{lineage}__{readout}__{stim}__{tp}__InfantIGRA_violin")),".png"))
  ggplot2::ggsave(fname, p, width = 6.4, height = wraps$height, units = "in",
                  dpi = 300, bg = "white")
}


comp3_top <- comp3 %>%
  mutate(Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel)) %>%
  arrange(Time_Point_std, Condition, desc(abs(median_diff_pp))) %>%
  group_by(Time_Point_std, Condition) %>%
  slice_head(n = 25) %>%
  ungroup()

invisible(pmap(
  list(comp3_top$Time_Point_std, comp3_top$Condition, comp3_top$Readout, comp3_top$LineageGate),
  ~ plot_comp3_violin(prolif, ..1, ..2, ..3, ..4, dirs$comp3_violin)
))



# =========================
# COMP 4 (Mothers): Pos vs Neg vs ATB
# =========================

# 4A) BOX/JITTER — per (Stim × Readout × Lineage)
plot_comp4_box <- function(df, stim, readout, lineage, outdir) {
  dat <- df %>%
    dplyr::filter(PID_Type == "Mother",
                  Condition == stim,
                  Readout == readout,
                  LineageGate == lineage,
                  !is.na(Maternal_IGRA),
                  Maternal_IGRA %in% c("Pos","Neg","ATB")) %>%
    dplyr::select(Maternal_IGRA, Value_for_test, DenominatorLabel)
  
  # need at least 2 groups present & ≥3 total obs
  if (nrow(dat) < 3 || length(unique(dat$Maternal_IGRA)) < 2) return(invisible(NULL))
  
  # KW overall
  kw   <- suppressWarnings(stats::kruskal.test(Value_for_test ~ Maternal_IGRA, data = dat))
  p_kw <- kw$p.value
  
  # Counts and axis labels with n
  lvls <- c("Pos","Neg","ATB")
  cnt  <- dat %>%
    dplyr::count(Maternal_IGRA) %>%
    dplyr::right_join(tibble::tibble(Maternal_IGRA = lvls), by = "Maternal_IGRA") %>%
    dplyr::mutate(n = dplyr::coalesce(n, 0L))
  axis_labels <- setNames(paste0(cnt$Maternal_IGRA, " (n=", cnt$n, ")"), cnt$Maternal_IGRA)
  
  denom_lbl <- unique(dat$DenominatorLabel)[1]
  title_txt <- glue::glue("{mk_readout_label(lineage, readout, denom_lbl)} — {stim}")
  subt_txt  <- glue::glue("Kruskal–Wallis p={signif(p_kw,3)} {sig_stars(p_kw)}")
  
  wraps <- wrap_title_sub(title_txt, subt_txt,
                          title_width = 90, sub_width = 110,
                          base_height = 4.8, extra_per_line = 0.28)
  
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = Maternal_IGRA, y = Value_for_test, fill = Maternal_IGRA)) +
    ggplot2::geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.95) +
    ggplot2::geom_jitter(width = 0.13, height = 0, size = 1.6, alpha = 0.7) +
    ggplot2::scale_y_continuous(limits = c(0,100),
                                expand = ggplot2::expansion(mult = c(0.02, 0.06))) +
    ggplot2::scale_x_discrete(limits = lvls, labels = axis_labels, drop = FALSE) +
    ggplot2::coord_cartesian(clip = "off") +  # anti-clipping for long titles/subtitles
    ggplot2::labs(
      title = wraps$title,
      subtitle = wraps$subtitle,
      x = NULL,
      y = "% of Parent (batch-normalized)"
    ) +
    anti_clip_theme()
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  fname <- file.path(outdir, paste0(
    sanitize_for_path(glue::glue("{lineage}__{readout}__{stim}__Mothers_KW_box")), ".png"))
  ggplot2::ggsave(fname, p, width = 6.8, height = wraps$height, units = "in",
                  dpi = 300, bg = "white")
}

# Plot boxes for top effects by KW p (from comp4_overall)
if (exists("comp4_overall") && nrow(comp4_overall) > 0) {
  comp4_top <- comp4_overall %>%
    dplyr::filter(!is.na(kw_p)) %>%
    dplyr::mutate(Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel)) %>%
    dplyr::arrange(Condition, kw_p) %>%
    dplyr::group_by(Condition) %>%
    dplyr::slice_head(n = 25) %>%
    dplyr::ungroup()
  
  invisible(purrr::pmap(
    list(comp4_top$Condition, comp4_top$Readout, comp4_top$LineageGate),
    ~ plot_comp4_box(prolif, ..1, ..2, ..3, dirs$comp4_box)
  ))
}

######## Summary Heatmaps
# Single output folder for all summary heatmaps
dirs$summary_heatmaps <- file.path(plots_root, "Summary_Heatmaps")
if (!dir.exists(dirs$summary_heatmaps)) dir.create(dirs$summary_heatmaps, recursive = TRUE)

# Select union of (significant) and (top |effect|) per panel group
select_sig_plus_top <- function(df, group_vars, effect_col, p_col, top_n = 30, alpha = 0.05) {
  df <- df %>% dplyr::mutate(.is_sig = !!rlang::sym(p_col) < alpha)
  
  sig <- df %>% dplyr::filter(.is_sig)
  
  top <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::arrange(dplyr::desc(abs(!!rlang::sym(effect_col))), .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
  
  dplyr::bind_rows(sig, top) %>% dplyr::distinct()
}


# Compute a sensible height based on the max number of rows across facets
compute_heatmap_height <- function(df, facet_vars, row_label_col, base = 5.5, per_row = 0.18, min_h = 6) {
  if (length(facet_vars)) {
    facet_counts <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(facet_vars))) %>%
      dplyr::summarise(n_rows = dplyr::n_distinct(!!rlang::sym(row_label_col)), .groups = "drop")
    max_rows <- if (nrow(facet_counts)) max(facet_counts$n_rows) else dplyr::n_distinct(df[[row_label_col]])
  } else {
    max_rows <- dplyr::n_distinct(df[[row_label_col]])
  }
  max(min_h, base + per_row * max_rows)
}
### COMP 1 — Infants, paired 12→44 (per Condition)
if (nrow(comp1) > 0) {
  comp1_hm0 <- comp1 %>%
    dplyr::mutate(
      Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel),
      star          = sig_stars(p_value)   # raw p
    )
  
  comp1_hm <- select_sig_plus_top(
    df = comp1_hm0,
    group_vars = c("Condition"),
    effect_col = "median_diff_pp",
    p_col      = "p_value",
    top_n      = 30,
    alpha      = 0.05
  ) %>%
    dplyr::group_by(Condition) %>%
    dplyr::arrange(dplyr::desc(abs(median_diff_pp)), .by_group = TRUE) %>%
    dplyr::ungroup()
  
  h_in <- compute_heatmap_height(comp1_hm, facet_vars = c("Condition"), row_label_col = "Readout_label")
  
  p1 <- ggplot(comp1_hm,
               aes(x = "Effect", y = reorder(Readout_label, median_diff_pp), fill = median_diff_pp)) +
    geom_tile(color = "grey70", linewidth = 0.2) +
    geom_text(aes(label = ifelse(star == "", "", star)), size = 3.6, fontface = "bold") +
    scale_x_discrete(expand = expansion(add = 0.1)) +
    scale_fill_gradient2(
      low = "#1b9e77", mid = "grey92", high = "#d95f02", midpoint = 0,
      name = "Median Δ(44−12)"
    ) +
    facet_wrap(~Condition, scales = "free_y", ncol = 1) +
    labs(
      title    = wrap_str("Comp1: 12wks → 44wks paired change (Infants, Proliferating)", 95),
      subtitle = wrap_str("Tiles = median Δ(44−12). Text = significance stars (raw p).", 110),
      x = NULL, y = NULL
    ) +
    anti_clip_theme(12) +
    theme(
      strip.text   = element_text(face = "bold"),
      axis.text.y  = element_text(size = 8.5),
      legend.position = "right",
      legend.title    = element_text(size = 11, face = "bold"),
      legend.text     = element_text(size = 10)
    ) +
    guides(fill = guide_colorbar(barwidth = 1.5, barheight = 8))
  
  ggsave(file.path(dirs$summary_heatmaps, "Comp1_Heatmap_medianDiff_p.png"),
         p1, width = 8.0, height = h_in, dpi = 300, bg = "white")
}


#COMP 2 — Infants, Maternal IGRA Pos–Neg (per Timepoint × Condition)
if (nrow(comp2) > 0) {
  comp2_hm0 <- comp2 %>%
    dplyr::mutate(
      Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel),
      star          = sig_stars(p_value)   # raw p
    )
  
  comp2_hm <- select_sig_plus_top(
    df = comp2_hm0,
    group_vars = c("Time_Point_std", "Condition"),
    effect_col = "median_diff_pp",
    p_col      = "p_value",
    top_n      = 30,
    alpha      = 0.05
  ) %>%
    dplyr::group_by(Time_Point_std, Condition) %>%
    dplyr::arrange(dplyr::desc(abs(median_diff_pp)), .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Condition = factor(Condition, levels = c("BCG", "DosR")))
  
  h_in <- compute_heatmap_height(comp2_hm, facet_vars = c("Time_Point_std"), row_label_col = "Readout_label")
  
  p2 <- ggplot(comp2_hm,
               aes(x = Condition, y = reorder(Readout_label, median_diff_pp), fill = median_diff_pp)) +
    geom_tile(color = "grey70", linewidth = 0.2) +
    geom_text(aes(label = ifelse(star == "", "", star)), size = 3.6, fontface = "bold") +
    scale_x_discrete(drop = TRUE) +
    scale_fill_gradient2(
      low = "#4575b4", mid = "grey92", high = "#d73027", midpoint = 0,
      name = "Median (Pos−Neg)"
    ) +
    facet_grid(Time_Point_std ~ ., scales = "free_y") +
    labs(
      title    = wrap_str("Comp2: Maternal IGRA (Pos vs Neg), Infants (Proliferating)", 95),
      subtitle = wrap_str("Tiles = median (Pos−Neg). Text = significance stars (raw p).", 110),
      x = NULL, y = NULL
    ) +
    anti_clip_theme(12) +
    theme(
      axis.text.y     = element_text(size = 8.5),
      strip.text.y    = element_text(face = "bold", angle = 0),
      legend.position = "right",
      legend.title    = element_text(size = 11, face = "bold"),
      legend.text     = element_text(size = 10)
    ) +
    guides(fill = guide_colorbar(barwidth = 1.5, barheight = 8))
  
  ggsave(file.path(dirs$summary_heatmaps, "Comp2_Heatmap_medianDiff_p.png"),
         p2, width = 8.4, height = h_in, dpi = 300, bg = "white")
}


# COMP 3 — Infants, Infant IGRA Pos–Neg
if (nrow(comp3) > 0) {
  comp3_hm0 <- comp3 %>%
    dplyr::mutate(
      Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel),
      star          = sig_stars(p_value)   # raw p
    )
  
  comp3_hm <- select_sig_plus_top(
    df = comp3_hm0,
    group_vars = c("Time_Point_std", "Condition"),
    effect_col = "median_diff_pp",
    p_col      = "p_value",
    top_n      = 30,
    alpha      = 0.05
  ) %>%
    dplyr::group_by(Time_Point_std, Condition) %>%
    dplyr::arrange(dplyr::desc(abs(median_diff_pp)), .by_group = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Condition = factor(Condition, levels = c("BCG", "DosR")))
  
  h_in <- compute_heatmap_height(comp3_hm, facet_vars = c("Time_Point_std"), row_label_col = "Readout_label")
  
  p3 <- ggplot(comp3_hm,
               aes(x = Condition, y = reorder(Readout_label, median_diff_pp), fill = median_diff_pp)) +
    geom_tile(color = "grey70", linewidth = 0.2) +
    geom_text(aes(label = ifelse(star == "", "", star)), size = 3.6, fontface = "bold") +
    scale_x_discrete(drop = TRUE) +
    scale_fill_gradient2(
      low = "#4575b4", mid = "grey92", high = "#d73027", midpoint = 0,
      name = "Median (Pos−Neg)"
    ) +
    facet_grid(Time_Point_std ~ ., scales = "free_y") +
    labs(
      title    = wrap_str("Comp3: Infant IGRA (Pos vs Neg), Infants (Proliferating)", 95),
      subtitle = wrap_str("Tiles = median (Pos−Neg). Text = significance stars (raw p).", 110),
      x = NULL, y = NULL
    ) +
    anti_clip_theme(12) +
    theme(
      axis.text.y     = element_text(size = 8.5),
      strip.text.y    = element_text(face = "bold", angle = 0),
      legend.position = "right",
      legend.title    = element_text(size = 11, face = "bold"),
      legend.text     = element_text(size = 10)
    ) +
    guides(fill = guide_colorbar(barwidth = 1.5, barheight = 8))
  
  ggsave(file.path(dirs$summary_heatmaps, "Comp3_Heatmap_medianDiff_p.png"),
         p3, width = 8.4, height = h_in, dpi = 300, bg = "white")
}


#COMP 4 — Mothers, KW
if (exists("comp4_overall") && nrow(comp4_overall) > 0) {
  comp4_hm0 <- comp4_overall %>%
    dplyr::mutate(
      Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel),
      star          = sig_stars(kw_p)   # raw KW p
    )
  
  # keep sig + top (by smallest p) within each Condition block, then combine and plot with Condition on x
  comp4_hm <- select_sig_plus_top(
    df = comp4_hm0,
    group_vars = c("Condition"),
    effect_col = "kw_p",   # smaller is stronger; selection still uses p
    p_col      = "kw_p",
    top_n      = 30,
    alpha      = 0.05
  ) %>%
    dplyr::mutate(Condition = factor(Condition, levels = c("BCG", "DosR"))) %>%
    dplyr::group_by(Condition) %>%
    dplyr::arrange(kw_p, .by_group = TRUE) %>%
    dplyr::ungroup()
  
  # dynamic height based on total unique rows (union across conditions)
  h_in <- compute_heatmap_height(comp4_hm, facet_vars = character(0), row_label_col = "Readout_label")
  
  p4 <- ggplot(
    comp4_hm,
    aes(x = Condition,
        y = reorder(Readout_label, -kw_p),
        fill = -log10(kw_p))) +
    geom_tile(color = "grey70", linewidth = 0.2) +
    geom_text(aes(label = ifelse(star == "", "", star)), size = 3.6, fontface = "bold") +
    scale_x_discrete(drop = TRUE) +  # drop any stim not present
    scale_fill_viridis_c(
      name = expression(-log[10]~p),
      option = "plasma",
      na.value = "grey90"
    ) +
    labs(
      title    = wrap_str("Comp4: Mothers (Pos vs Neg vs ATB), Kruskal–Wallis", 95),
      subtitle = wrap_str("Tiles = −log10(p); text = significance stars (raw p).", 110),
      x = NULL, y = NULL
    ) +
    anti_clip_theme(12) +
    theme(
      axis.text.y     = element_text(size = 8.5),
      legend.position = "right",
      legend.title    = element_text(size = 11, face = "bold"),
      legend.text     = element_text(size = 10)
    ) +
    guides(fill = guide_colorbar(barwidth = 1.5, barheight = 8))
  
  ggsave(file.path(dirs$summary_heatmaps, "Comp4_Heatmap_KW_p.png"),
         p4, width = 8.6, height = h_in, dpi = 300, bg = "white")
}
