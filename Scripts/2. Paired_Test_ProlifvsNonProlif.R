#########################
# Libraries
###########
##############
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)

#########################
# Paths & I/O
#########################
#base_dir <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim"
base_dir <- "C:/Users/ammas/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim"
getwd()
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

# Helper to make safe folder/file names from readout labels
sanitize_for_path <- function(x) {
  x |>
    stringr::str_replace_all("[/\\\\|?:*<>%]", "_") |>  # remove path-breaking chars
    stringr::str_replace_all("[\\s]+", " ") |>         # collapse whitespace
    stringr::str_trim()
}
# =========================
# Build PID-level differences (pid_diffs)
# =========================

# helper to normalize Prolif_Status exactly like in paired_stats()
normalize_prolif <- function(x) {
  x <- stringr::str_squish(x)
  x <- stringr::str_replace_all(x, "[\u2010-\u2015\u2212]", "-")  # unicode dashes -> "-"
  dplyr::recode(x, "Non Proliferating" = "Non-Proliferating", .default = x)
}

pid_diffs <- cytokine_an %>%
  dplyr::mutate(Prolif_Status = factor(normalize_prolif(Prolif_Status),
                                       levels = c("Proliferating","Non-Proliferating"))) %>%
  # If there are replicate rows per PID/Status, summarise first to a single value per side
  dplyr::group_by(PID, Group, Condition, Compartment, Readout, Prolif_Status) %>%
  dplyr::summarise(Value_side = median(Value_for_test, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = Prolif_Status,
    values_from = Value_side,
    values_fill = NA_real_
  ) %>%
  dplyr::mutate(
    pair_complete = !is.na(Proliferating) & !is.na(`Non-Proliferating`),
    diff_pp       = Proliferating - `Non-Proliferating`,
    abs_diff_pp   = abs(diff_pp),
    direction     = dplyr::case_when(
      is.na(diff_pp)               ~ NA_character_,
      diff_pp > 0                  ~ "Proliferating > Non-Proliferating",
      diff_pp < 0                  ~ "Proliferating < Non-Proliferating",
      TRUE                         ~ "Equal"
    )
  ) %>%
  dplyr::filter(pair_complete)

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
# ========= PLOTTING (aligned lines, anti-clipping, dynamic sizing) =========
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
      
      # ------- data prep -------
      df_raw <- data_rd %>%
        dplyr::filter(Group == !!Group, Condition == !!Condition, Compartment == !!Compartment) %>%
        dplyr::select(PID, Prolif_Status, Value_for_test) %>%
        dplyr::distinct() %>%
        dplyr::mutate(
          Prolif_Status = stringr::str_squish(Prolif_Status),
          Prolif_Status = stringr::str_replace_all(Prolif_Status, "[\u2010-\u2015\u2212]", "-"),
          Prolif_Status = dplyr::recode(Prolif_Status,
                                        "Non Proliferating" = "Non-Proliferating",
                                        .default = Prolif_Status),
          Prolif_Status = factor(Prolif_Status, levels = c("Non-Proliferating","Proliferating")),
          Value_for_test = pmin(pmax(Value_for_test, 0), 100)
        )
      
      # one value per PID × Status (ensures single point per side)
      df_plot <- df_raw %>%
        dplyr::group_by(PID, Prolif_Status) %>%
        dplyr::summarise(Value_for_test = median(Value_for_test, na.rm = TRUE), .groups = "drop")
      
      # wide for pairing & n_pairs; use df_plot so stats match what we draw
      wide <- df_plot %>%
        tidyr::pivot_wider(
          names_from  = Prolif_Status,
          values_from = Value_for_test,
          values_fill = NA_real_
        )
      
      if (!("Proliferating" %in% names(wide)))         wide$Proliferating <- NA_real_
      if (!("Non-Proliferating" %in% names(wide))) wide$`Non-Proliferating` <- NA_real_
      
      wide <- wide %>% dplyr::filter(!is.na(Proliferating), !is.na(`Non-Proliferating`))
      n_pairs <- nrow(wide)
      if (n_pairs < 1) return(invisible(NULL))
      
      # ------- stats for subtitle -------
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
      
      # ------- titles (wrapped) & sizing -------
      title_raw      <- glue("{rd} — {Group} — {Condition} — {Compartment}")
      title_wrapped  <- stringr::str_wrap(title_raw, width = 70)
      
      subtitle_txt   <- glue(
        "n={n_pairs} | median Δ(P−NP)={round(median_diff, 3)} | p={signif(pval, 3)}",
        if (!is.na(padj)) glue(" | FDR={signif(padj,3)}") else ""
      )
      subtitle_wrapped <- stringr::str_wrap(subtitle_txt, width = 85)
      
      n_title_lines <- max(1, length(unlist(strsplit(title_wrapped, "\n", fixed = TRUE))))
      n_sub_lines   <- max(1, length(unlist(strsplit(subtitle_wrapped, "\n", fixed = TRUE))))
      
      w_in <- max(4.8, min(8.5, 4.6 + 0.07 * n_pairs))
      h_in <- 4.6 + 0.22 * (n_title_lines - 1 + n_sub_lines - 1)
      
      # ------- plot (no jitter; perfect alignment) -------
      p <- ggplot(df_plot, aes(x = Prolif_Status, y = Value_for_test, group = PID)) +
        geom_line(alpha = 0.4, linewidth = 0.5, lineend = "round", position = position_identity()) +
        geom_point(aes(color = Prolif_Status), size = 2, position = position_identity()) +
        stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 0, alpha = 0.85) +
        scale_x_discrete(drop = FALSE) +
        coord_cartesian(ylim = c(0, 100), clip = "off") +
        labs(
          title = title_wrapped,
          subtitle = subtitle_wrapped,
          x = NULL,
          y = "% of Parent (batch-normalized)"
        ) +
        theme_bw(base_size = 12) +
        theme(
          legend.position = "none",
          plot.title.position = "plot",
          plot.title   = element_text(face = "bold", margin = margin(b = 6)),
          plot.subtitle= element_text(margin = margin(b = 8)),
          panel.grid.major.x = element_blank(),
          plot.margin = margin(t = 14, r = 18, b = 14, l = 18, unit = "pt")
        )
      
      # ------- save -------
      fname <- file.path(
        subdir,
        paste0(
          sanitize_for_path(glue("{rd}__{Group}__{Condition}__{Compartment}")),
          ".png"
        )
      )
      
      ggplot2::ggsave(
        filename = fname, plot = p,
        width = w_in, height = h_in, units = "in",
        dpi = 300, bg = "white"
      )
    }
  )
})

# ==== Overview Heatmap: Δ(Proliferating − Non-Proliferating) ====
library(ggplot2)
library(dplyr)
library(glue)

# Output dir for overview figs
overview_dir <- file.path(res_root, "Overviews")
if (!dir.exists(overview_dir)) dir.create(overview_dir, recursive = TRUE)

# Optional: lock a specific order for facets/conditions
# stats_all$Group     <- factor(stats_all$Group,     levels = c("Mothers Entry","Mothers Dx","Infants 12wks","Infants 44wks"))
# stats_all$Condition <- factor(stats_all$Condition, levels = c("BCG","DosR","E6C10","GAG"))

# Basic size heuristics (adjust multipliers if needed)
n_readouts <- stats_all %>% distinct(Readout) %>% nrow()
n_groups   <- stats_all %>% distinct(Group)   %>% nrow()
h_in_heat  <- max(6, 0.14 * n_readouts)   # more readouts -> taller
w_in_heat  <- 9

# Optional: add a significance flag
heat_df <- stats_all %>%
  mutate(sig = !is.na(p_adj) & p_adj < 0.05)

p_heat <- ggplot(
  heat_df,
  aes(x = Condition, y = Readout, fill = median_diff_pp)
) +
  geom_tile(color = "grey70", linewidth = 0.2) +
  # optional small dot for significant cells:
  geom_point(
    data = subset(heat_df, sig),
    aes(x = Condition, y = Readout),
    inherit.aes = FALSE, shape = 21, size = 1.5, stroke = 0.2, fill = "black", alpha = 0.7
  ) +
  facet_wrap(~ Group, ncol = 1, scales = "free_y") +
  scale_fill_viridis_c(option = "plasma", name = "Δ(P−NP)") +
  labs(
    title = "Proliferation bias by readout (Δ(Proliferating − Non-Proliferating))",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", margin = margin(b = 8)),
    plot.margin = margin(t = 12, r = 16, b = 12, l = 16)
  )

# Save PNG + PDF
ggsave(file.path(overview_dir, "Heatmap_DeltaProlif_vs_NonProlif.png"),
       p_heat, width = 20, height = 14, units = "in", dpi = 300, bg = "white")

