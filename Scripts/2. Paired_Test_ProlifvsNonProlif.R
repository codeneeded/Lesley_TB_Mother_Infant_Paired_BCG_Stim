#########################
# Libraries
#########################
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)

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

