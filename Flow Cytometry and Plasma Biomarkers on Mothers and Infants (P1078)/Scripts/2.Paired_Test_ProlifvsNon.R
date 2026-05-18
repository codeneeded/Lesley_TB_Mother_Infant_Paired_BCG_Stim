## ============================================================================
## Project: Flow Cytometry and Plasma Biomarkers on Mothers and Infants (P1078)
## Script:  2_Paired_Test_ProlifvsNonProlif.R
## Purpose: Paired Wilcoxon tests of Proliferating vs Non-Proliferating frequencies,
##          and Prolif-only between-group comparisons (Comp1–Comp4).
##
## Key fixes vs prior version:
##   - Uses the project-root path layout (no more Linux/Windows mix).
##   - `pid_diffs` (long table of per-PID Prolif − NonProlif deltas) is now
##     actually built before the loop that writes it.
##   - Helper functions defined ONCE near the top, not redefined further down.
##   - Comp2/3/4 dedupe to one row per PID within each test group before
##     running the Wilcoxon, so technical replicates can't inflate n.
##   - Heatmap blocks check that input tables are not empty before plotting.
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(stringr)
  library(readr); library(ggplot2); library(glue); library(qs2)
})

## ---- Paths -----------------------------------------------------------------
proj_root <- "/home/akshay-iyer/Documents/Lesley_TB_Mother_Infant_Paired_BCG_Stim/Flow Cytometry and Plasma Biomarkers on Mothers and Infants (P1078)"
in_qs2    <- file.path(proj_root, "saved_R_dat", "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.qs2")
in_csv    <- file.path(proj_root, "saved_R_dat", "TB_Flow_Cytokine_DATA_102725_preprocessed_MEDnorm.csv")

cytokine_long <- if (file.exists(in_qs2)) {
  qs2::qs_read(in_qs2)
} else if (file.exists(in_csv)) {
  read_csv(in_csv, show_col_types = FALSE)
} else {
  stop("Could not find preprocessed input at:\n  ", in_qs2, "\n  ", in_csv)
}

res_root <- file.path(proj_root, "Results", "Flow_Prolif_vs_NonProlif_PairedTests")
dir.create(res_root, showWarnings = FALSE, recursive = TRUE)

## ---- Helpers (defined ONCE) -----------------------------------------------
sanitize_for_path <- function(x) {
  x |>
    stringr::str_replace_all("[/\\\\|?:*<>%\"]", "_") |>
    stringr::str_replace_all("[\\s]+", " ") |>
    stringr::str_trim()
}

wrap_str <- function(x, width = 80) stringr::str_wrap(x, width = width)

count_lines <- function(x) {
  if (is.null(x) || any(is.na(x))) return(1L)
  length(unlist(strsplit(x, "\n", fixed = TRUE)))
}

wrap_title_sub <- function(title, subtitle, title_width = 80, sub_width = 100,
                           base_height = 4.6, extra_per_line = 0.28) {
  tw <- wrap_str(title, title_width)
  sw <- wrap_str(subtitle, sub_width)
  n_title <- max(1L, count_lines(tw))
  n_sub   <- max(1L, count_lines(sw))
  h <- base_height + extra_per_line * ((n_title - 1L) + (n_sub - 1L))
  list(title = tw, subtitle = sw, height = h)
}

anti_clip_theme <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) %+replace%
    ggplot2::theme(
      legend.position    = "none",
      plot.title.position = "plot",
      plot.title         = ggplot2::element_text(face = "bold", margin = ggplot2::margin(b = 8)),
      plot.subtitle      = ggplot2::element_text(margin = ggplot2::margin(b = 10)),
      plot.caption       = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(t = 20, r = 28, b = 20, l = 20)
    )
}

sig_stars <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01,  "**",
                       ifelse(p < 0.05,  "*", ""))))
}

mk_readout_label <- function(lineage, readout, denom_label) {
  lineage <- ifelse(is.na(lineage) | lineage == "", "Lineage?", lineage)
  rl      <- ifelse(is.na(readout) | readout == "", "Readout?", readout)
  paste0(lineage, " :: ", rl, " [", denom_label, "]")
}

## ---- Analysis settings -----------------------------------------------------
stim_keep <- c("BCG", "DosR", "E6C10", "GAG")
T_cols    <- grep("^T[0-9]+$", names(cytokine_long), value = TRUE)

## ---- Build the four analysis groups ---------------------------------------
cytokine_an <- cytokine_long %>%
  filter(!is.na(PID_Type), PID_Type != "Batch Control") %>%
  filter(Condition %in% stim_keep) %>%
  mutate(Value_for_test = Value_RawNormalized) %>%
  filter(Metric == "Freq. of Parent") %>%
  filter(Prolif_Status %in% c("Proliferating","Non-Proliferating")) %>%
  mutate(
    Readout = if (length(T_cols) > 0) {
      pmap_chr(across(all_of(T_cols)), function(...) {
        toks <- c(...); toks <- toks[!is.na(toks) & toks != ""]
        if (length(toks) == 0) NA_character_ else str_squish(paste(toks, collapse = " / "))
      })
    } else NA_character_
  ) %>%
  filter(!is.na(Readout)) %>%
  mutate(
    Time_Point     = str_squish(Time_Point),
    Time_Point_std = dplyr::recode(Time_Point,
                                   "12 Wks" = "12wks", "44 Wks" = "44wks",
                                   "Entry"  = "Entry", "Dx"     = "Dx",
                                   .default = Time_Point),
    Group = dplyr::case_when(
      PID_Type == "Mother" & Time_Point_std == "Entry" ~ "Mothers Entry",
      PID_Type == "Mother" & Time_Point_std == "Dx"    ~ "Mothers Dx",
      PID_Type == "Infant" & Time_Point_std == "12wks" ~ "Infants 12wks",
      PID_Type == "Infant" & Time_Point_std == "44wks" ~ "Infants 44wks",
      TRUE ~ NA_character_
    ),
    Compartment = ifelse(is.na(Compartment), "Unknown", Compartment)
  ) %>%
  filter(!is.na(Group))

## ---- Paired Wilcoxon helper (Prolif vs Non-Prolif within PID) -------------
paired_stats <- function(df) {
  df2 <- df %>%
    mutate(
      Prolif_Status = stringr::str_squish(Prolif_Status),
      Prolif_Status = stringr::str_replace_all(Prolif_Status, "[\u2010-\u2015\u2212]", "-"),
      Prolif_Status = dplyr::recode(Prolif_Status,
                                    "Non Proliferating" = "Non-Proliferating",
                                    .default = Prolif_Status),
      Prolif_Status = factor(Prolif_Status, levels = c("Proliferating","Non-Proliferating"))
    )
  
  # collapse multiple rows per (PID, status) to a single value (median, robust)
  wide <- df2 %>%
    dplyr::group_by(PID, Prolif_Status) %>%
    dplyr::summarise(Value_for_test = stats::median(Value_for_test, na.rm = TRUE),
                     .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Prolif_Status, values_from = Value_for_test,
                       values_fill = NA_real_)
  
  if (!("Proliferating"     %in% names(wide))) wide$Proliferating <- NA_real_
  if (!("Non-Proliferating" %in% names(wide))) wide$`Non-Proliferating` <- NA_real_
  
  wide <- wide %>% dplyr::filter(!is.na(Proliferating), !is.na(`Non-Proliferating`))
  n_pairs <- nrow(wide)
  if (n_pairs < 1) {
    return(dplyr::tibble(
      n_pairs = 0, median_Prolif = NA_real_, median_NonProlif = NA_real_,
      median_diff_pp = NA_real_, wilcox_stat = NA_real_, p_value = NA_real_,
      effect_rbc = NA_real_
    ))
  }
  
  w <- suppressWarnings(stats::wilcox.test(
    x = wide$Proliferating, y = wide$`Non-Proliferating`,
    paired = TRUE, exact = FALSE, conf.int = FALSE
  ))
  V <- unname(w$statistic)
  r_rb <- if (!is.null(V)) 1 - (2 * as.numeric(V)) / (n_pairs * (n_pairs + 1)) else NA_real_
  
  dplyr::tibble(
    n_pairs           = n_pairs,
    median_Prolif     = stats::median(wide$Proliferating, na.rm = TRUE),
    median_NonProlif  = stats::median(wide$`Non-Proliferating`, na.rm = TRUE),
    median_diff_pp    = stats::median(wide$Proliferating - wide$`Non-Proliferating`, na.rm = TRUE),
    wilcox_stat       = as.numeric(V),
    p_value           = w$p.value,
    effect_rbc        = r_rb
  )
}

## ---- Run paired stats across all combos -----------------------------------
stats_all <- cytokine_an %>%
  group_by(Readout, Group, Condition, Compartment) %>%
  group_modify(~ paired_stats(.x)) %>%
  ungroup() %>%
  group_by(Readout) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write_csv(stats_all, file.path(res_root, "_ALL_readouts_summary.csv"))

## ---- Build pid_diffs (per-PID paired differences) -------------------------
##   This is what the per-readout loop below needs. Build once.
pid_diffs <- cytokine_an %>%
  mutate(
    Prolif_Status = stringr::str_squish(Prolif_Status),
    Prolif_Status = stringr::str_replace_all(Prolif_Status, "[\u2010-\u2015\u2212]", "-"),
    Prolif_Status = dplyr::recode(Prolif_Status,
                                  "Non Proliferating" = "Non-Proliferating",
                                  .default = Prolif_Status)
  ) %>%
  group_by(Readout, Group, Condition, Compartment, PID, Prolif_Status) %>%
  summarise(Value_for_test = stats::median(Value_for_test, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = Prolif_Status, values_from = Value_for_test,
              values_fill = NA_real_) %>%
  filter(!is.na(Proliferating), !is.na(`Non-Proliferating`)) %>%
  mutate(diff_P_minus_NP = Proliferating - `Non-Proliferating`)

## ---- Per-readout subfolders -----------------------------------------------
readouts <- cytokine_an %>% distinct(Readout) %>% pull(Readout) %>% sort()

purrr::walk(readouts, function(rd) {
  safe_rd <- sanitize_for_path(rd)
  subdir  <- file.path(res_root, safe_rd)
  dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
  
  stats_rd <- stats_all %>% dplyr::filter(Readout == rd)
  write_csv(stats_rd, file.path(subdir, paste0(safe_rd, "_paired_stats.csv")))
  
  diffs_rd <- pid_diffs %>% dplyr::filter(Readout == rd)
  write_csv(diffs_rd, file.path(subdir, paste0(safe_rd, "_paired_differences.csv")))
})

## ---- Per-readout slope plots (Prolif vs NonProlif) ------------------------
cytokine_an <- cytokine_an %>%
  mutate(Prolif_Status = factor(Prolif_Status,
                                levels = c("Non-Proliferating","Proliferating")))

purrr::walk(readouts, function(rd) {
  safe_rd <- sanitize_for_path(rd)
  subdir  <- file.path(res_root, safe_rd)
  
  data_rd <- cytokine_an %>% filter(Readout == rd)
  combos  <- data_rd %>% distinct(Group, Condition, Compartment)
  
  purrr::pwalk(combos, function(Group, Condition, Compartment) {
    df <- data_rd %>%
      filter(Group == !!Group, Condition == !!Condition, Compartment == !!Compartment) %>%
      group_by(PID, Prolif_Status) %>%
      summarise(Value_for_test = stats::median(Value_for_test, na.rm = TRUE),
                .groups = "drop")
    
    wide <- df %>%
      tidyr::pivot_wider(names_from = Prolif_Status, values_from = Value_for_test,
                         values_fill = NA_real_)
    if (!("Proliferating"     %in% names(wide))) wide$Proliferating <- NA_real_
    if (!("Non-Proliferating" %in% names(wide))) wide$`Non-Proliferating` <- NA_real_
    wide <- wide %>% filter(!is.na(Proliferating), !is.na(`Non-Proliferating`))
    n_pairs <- nrow(wide)
    if (n_pairs < 1) return(invisible(NULL))
    
    stat_row <- stats_all %>%
      filter(Readout == !!rd, Group == !!Group,
             Condition == !!Condition, Compartment == !!Compartment) %>%
      slice(1)
    if (nrow(stat_row) == 0) {
      w <- suppressWarnings(stats::wilcox.test(wide$Proliferating, wide$`Non-Proliferating`,
                                               paired = TRUE, exact = FALSE))
      median_diff <- stats::median(wide$Proliferating - wide$`Non-Proliferating`, na.rm = TRUE)
      pval <- w$p.value; padj <- NA_real_
    } else {
      median_diff <- stat_row$median_diff_pp
      pval        <- stat_row$p_value
      padj        <- stat_row$p_adj
    }
    
    subtitle_txt <- glue("n_pairs = {n_pairs} | median Δ(P − NP) = {round(median_diff, 3)} | p = {signif(pval, 3)}{ifelse(is.na(padj),'', paste0(' | FDR = ', signif(padj,3)))}")
    
    w_in <- max(4.5, min(7, 4.5 + 0.05 * n_pairs))
    h_in <- 4.5
    
    p <- ggplot(df, aes(x = Prolif_Status, y = Value_for_test, group = PID)) +
      geom_line(alpha = 0.35) +
      geom_point(aes(color = Prolif_Status), size = 2) +
      stat_summary(fun = stats::median, geom = "crossbar", width = 0.5, fatten = 0, alpha = 0.7) +
      scale_x_discrete(drop = FALSE) +
      labs(title = glue("{rd} — {Group} — {Condition} — {Compartment}"),
           subtitle = subtitle_txt,
           x = NULL, y = "% of Parent (batch-normalized)") +
      anti_clip_theme()
    
    fname <- file.path(subdir, paste0(
      sanitize_for_path(glue("{rd}__{Group}__{Condition}__{Compartment}")), ".png"))
    ggsave(fname, plot = p, width = w_in, height = h_in, dpi = 300, bg = "white")
  })
})

## ============================================================================
## PROLIFERATING-ONLY COMPARISONS (Comp1–Comp4)
## ============================================================================
comp_dir <- file.path(proj_root, "Results", "Flow_ProlifOnly_GroupComparisons", "StatTables_CSV")
dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

prolif <- cytokine_long %>%
  filter(!is.na(PID_Type), PID_Type != "Batch Control") %>%
  mutate(Condition = str_squish(Condition)) %>%
  filter(Prolif_Status == "Proliferating",
         Metric %in% c("Freq. of Parent", "Freq. of Grandparent")) %>%
  mutate(
    PathCore    = str_replace(GatePath, "\\s*\\|\\s*Freq\\..*$", "") %>%
      str_replace_all(",", " ") %>% str_squish(),
    LineageGate = map_chr(str_split(PathCore, "/"), function(segs) {
      segs <- segs[segs != ""]
      idx  <- which(segs == "Proliferating")
      if (length(idx) && idx[1] > 1) segs[idx[1] - 1] else NA_character_
    }),
    LineageGate = case_when(
      LineageGate %in% c("gd TCR","gdTCR") ~ "gd TCR",
      TRUE                                  ~ LineageGate
    ),
    Time_Point     = str_squish(Time_Point),
    Time_Point_std = recode(Time_Point,
                            "12 Wks" = "12wks", "44 Wks" = "44wks",
                            "Entry"  = "Entry", "Dx"     = "Dx",
                            .default = Time_Point),
    Compartment    = coalesce(Compartment, LineageGate),
    Readout = if (length(T_cols) > 0) {
      pmap_chr(across(all_of(T_cols)), function(...) {
        toks <- c(...); toks <- toks[!is.na(toks) & toks != ""]
        if (length(toks) == 0) NA_character_ else str_squish(paste(toks, collapse = " / "))
      })
    } else NA_character_,
    Value_for_test = pmin(pmax(Value_RawNormalized, 0), 100),
    DenominatorLabel = case_when(
      Metric == "Freq. of Parent"      ~ paste0("Freq. of Parent: ",      DenominatorGate),
      Metric == "Freq. of Grandparent" ~ paste0("Freq. of Grandparent: ", DenominatorGate),
      TRUE                              ~ Metric
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
  rename(Maternal_IGRA = `Maternal IGRA`, Infant_IGRA = `Infant IGRA`) %>%
  select(PID, PID_Type, Group, Time_Point_std, Condition,
         Metric, DenominatorGate, DenominatorLabel,
         ParentGate, GrandparentGate, LineageGate, Compartment,
         Readout, Maternal_IGRA, Infant_IGRA, Value_for_test)

## ---- Effect-size helpers --------------------------------------------------
cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  diffs <- outer(x, y, "-")
  (sum(diffs > 0) - sum(diffs < 0)) / (length(x) * length(y))
}

# Helper: collapse to one row per PID before unpaired tests (median across replicates)
collapse_per_pid <- function(df, grp_col) {
  df %>%
    dplyr::filter(!is.na(.data[[grp_col]])) %>%
    dplyr::group_by(PID, .data[[grp_col]]) %>%
    dplyr::summarise(Value_for_test = stats::median(Value_for_test, na.rm = TRUE),
                     .groups = "drop")
}

## ---- Comp1: Infants paired 12wks → 44wks within stim (BCG, DosR) ---------
comp1 <- prolif %>%
  filter(PID_Type == "Infant", Condition %in% c("BCG","DosR"),
         Time_Point_std %in% c("12wks","44wks")) %>%
  group_by(Readout, Condition, Metric, DenominatorGate, DenominatorLabel, LineageGate) %>%
  group_modify(~{
    wide <- .x %>%
      group_by(PID, Time_Point_std) %>%
      summarise(Value = stats::median(Value_for_test, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Time_Point_std, values_from = Value) %>%
      filter(is.finite(`12wks`), is.finite(`44wks`))
    n_pairs <- nrow(wide)
    if (n_pairs < 1) {
      tibble(n_pairs = 0, median_12w = NA_real_, median_44w = NA_real_,
             median_diff_pp = NA_real_, p_value = NA_real_, effect_rbc = NA_real_)
    } else {
      w <- suppressWarnings(wilcox.test(wide$`12wks`, wide$`44wks`, paired = TRUE, exact = FALSE))
      V <- unname(w$statistic)
      tibble(
        n_pairs        = n_pairs,
        median_12w     = stats::median(wide$`12wks`, na.rm = TRUE),
        median_44w     = stats::median(wide$`44wks`, na.rm = TRUE),
        median_diff_pp = stats::median(wide$`44wks` - wide$`12wks`, na.rm = TRUE),
        p_value        = w$p.value,
        effect_rbc     = if (is.null(V)) NA_real_ else 1 - (2*as.numeric(V))/(n_pairs*(n_pairs+1))
      )
    }
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(Condition, Metric, DenominatorGate, LineageGate, desc(abs(median_diff_pp)))

write_csv(comp1, file.path(comp_dir, "Comp1_Infants_Paired_12w_vs_44w_byStim.csv"))

## ---- Comp2: Infants by Maternal IGRA (Pos vs Neg) by TP × Condition ------
comp2 <- prolif %>%
  filter(PID_Type == "Infant",
         Time_Point_std %in% c("12wks","44wks"),
         !is.na(Maternal_IGRA), Maternal_IGRA %in% c("Pos","Neg")) %>%
  group_by(Readout, Time_Point_std, Condition, Metric, DenominatorGate, DenominatorLabel, LineageGate) %>%
  group_modify(~{
    dat  <- collapse_per_pid(.x, "Maternal_IGRA")
    gpos <- dat %>% filter(Maternal_IGRA == "Pos") %>% pull(Value_for_test)
    gneg <- dat %>% filter(Maternal_IGRA == "Neg") %>% pull(Value_for_test)
    if (!length(gpos) || !length(gneg)) {
      tibble(n_pos = length(gpos), n_neg = length(gneg),
             median_pos = NA_real_, median_neg = NA_real_,
             median_diff_pp = NA_real_, p_value = NA_real_, cliffs_delta = NA_real_)
    } else {
      w <- suppressWarnings(wilcox.test(gpos, gneg, paired = FALSE, exact = FALSE))
      tibble(
        n_pos = length(gpos), n_neg = length(gneg),
        median_pos = stats::median(gpos, na.rm = TRUE),
        median_neg = stats::median(gneg, na.rm = TRUE),
        median_diff_pp = stats::median(gpos, na.rm = TRUE) - stats::median(gneg, na.rm = TRUE),
        p_value = w$p.value,
        cliffs_delta = cliffs_delta(gpos, gneg)
      )
    }
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(Time_Point_std, Condition, Metric, DenominatorGate, LineageGate, desc(abs(median_diff_pp)))

write_csv(comp2, file.path(comp_dir, "Comp2_Infants_MaternalIGRA_Pos_vs_Neg_byTP_byStim.csv"))

## ---- Comp3: Infants by Infant IGRA (Pos vs Neg) --------------------------
comp3 <- prolif %>%
  filter(PID_Type == "Infant",
         Time_Point_std %in% c("12wks","44wks"),
         !is.na(Infant_IGRA), Infant_IGRA %in% c("Pos","Neg")) %>%
  group_by(Readout, Time_Point_std, Condition, Metric, DenominatorGate, DenominatorLabel, LineageGate) %>%
  group_modify(~{
    dat  <- collapse_per_pid(.x, "Infant_IGRA")
    gpos <- dat %>% filter(Infant_IGRA == "Pos") %>% pull(Value_for_test)
    gneg <- dat %>% filter(Infant_IGRA == "Neg") %>% pull(Value_for_test)
    if (!length(gpos) || !length(gneg)) {
      tibble(n_pos = length(gpos), n_neg = length(gneg),
             median_pos = NA_real_, median_neg = NA_real_,
             median_diff_pp = NA_real_, p_value = NA_real_, cliffs_delta = NA_real_)
    } else {
      w <- suppressWarnings(wilcox.test(gpos, gneg, paired = FALSE, exact = FALSE))
      tibble(
        n_pos = length(gpos), n_neg = length(gneg),
        median_pos = stats::median(gpos, na.rm = TRUE),
        median_neg = stats::median(gneg, na.rm = TRUE),
        median_diff_pp = stats::median(gpos, na.rm = TRUE) - stats::median(gneg, na.rm = TRUE),
        p_value = w$p.value,
        cliffs_delta = cliffs_delta(gpos, gneg)
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

## ---- Comp4: Mothers KW (Pos vs Neg vs ATB) -------------------------------
comp4_overall <- prolif %>%
  filter(PID_Type == "Mother",
         !is.na(Maternal_IGRA), Maternal_IGRA %in% c("Pos","Neg","ATB")) %>%
  group_by(Readout, Condition, Metric, DenominatorGate, DenominatorLabel, LineageGate) %>%
  group_modify(~{
    dat <- collapse_per_pid(.x, "Maternal_IGRA") %>% drop_na()
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
         !is.na(Maternal_IGRA), Maternal_IGRA %in% c("Pos","Neg","ATB")) %>%
  group_by(Readout, Condition, Metric, DenominatorGate, DenominatorLabel, LineageGate) %>%
  group_modify(~{
    dat <- collapse_per_pid(.x, "Maternal_IGRA") %>% drop_na()
    if (length(unique(dat$Maternal_IGRA)) < 2) {
      tibble(comp = character(), p = numeric())
    } else {
      pw  <- suppressWarnings(pairwise.wilcox.test(
        dat$Value_for_test, dat$Maternal_IGRA, p.adjust.method = "BH", exact = FALSE))
      mat <- pw$p.value
      if (is.null(mat)) { tibble(comp = character(), p = numeric())
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

write_csv(comp4_overall,  file.path(comp_dir, "Comp4_Mothers_KruskalWallis_byStim.csv"))
write_csv(comp4_pairwise, file.path(comp_dir, "Comp4_Mothers_PairwiseWilcoxon_byStim.csv"))

## ============================================================================
## PLOTTING SUITE — Proliferating Only
## ============================================================================
plots_root <- file.path(proj_root, "Results", "Flow_ProlifOnly_GroupComparisons", "Plots")
dirs <- list(
  comp1_slopes     = file.path(plots_root, "Infants_12wk_vs_44wk_Paired"),
  comp2_violin     = file.path(plots_root, "Infants_byMaternalIGRA_PosVsNeg"),
  comp3_violin     = file.path(plots_root, "Infants_byInfantIGRA_PosVsNeg"),
  comp4_box        = file.path(plots_root, "Mothers_byMaternalIGRA_PosNegATB"),
  summary_heatmaps = file.path(plots_root, "Summary_Heatmaps_AllComparisons")
)
invisible(lapply(dirs, function(d) dir.create(d, recursive = TRUE, showWarnings = FALSE)))

## ---- Comp1 slope plots ----------------------------------------------------
plot_comp1_slope <- function(df, condition, readout, lineage, outdir) {
  dat <- df %>%
    filter(PID_Type == "Infant", Condition == condition,
           Time_Point_std %in% c("12wks","44wks"),
           Readout == readout, LineageGate == lineage) %>%
    group_by(PID, Time_Point_std) %>%
    summarise(Value = stats::median(Value_for_test, na.rm = TRUE),
              DenominatorLabel = dplyr::first(DenominatorLabel), .groups = "drop") %>%
    pivot_wider(names_from = Time_Point_std, values_from = Value) %>%
    filter(is.finite(`12wks`), is.finite(`44wks`)) %>%
    mutate(diff = `44wks` - `12wks`)
  n_pairs <- nrow(dat); if (n_pairs < 1) return(invisible(NULL))
  
  w <- suppressWarnings(stats::wilcox.test(dat$`12wks`, dat$`44wks`, paired = TRUE, exact = FALSE))
  med12   <- stats::median(dat$`12wks`, na.rm = TRUE)
  med44   <- stats::median(dat$`44wks`, na.rm = TRUE)
  meddiff <- stats::median(dat$diff, na.rm = TRUE)
  pval    <- w$p.value
  
  denom_lbl <- unique(dat$DenominatorLabel)[1]
  title_txt <- glue("{mk_readout_label(lineage, readout, denom_lbl)} — {condition}")
  subt_txt  <- glue("n={n_pairs} | med 12w={round(med12,2)}, 44w={round(med44,2)} | Δ={round(meddiff,2)} | p={signif(pval,3)} {sig_stars(pval)}")
  
  long <- dat |>
    select(PID, `12wks`, `44wks`) |>
    pivot_longer(cols = c(`12wks`,`44wks`), names_to = "Time", values_to = "Value") |>
    mutate(Time = factor(Time, levels = c("12wks","44wks")))
  
  wraps <- wrap_title_sub(title_txt, subt_txt, 90, 110, 4.8, 0.28)
  w_in  <- max(5.2, min(8.6, 5.2 + 0.05*n_pairs))
  
  p <- ggplot(long, aes(x = Time, y = Value, group = PID)) +
    geom_line(alpha = 0.38, linewidth = 0.55) +
    geom_point(aes(color = Time), size = 2) +
    stat_summary(fun = stats::median, geom = "crossbar", width = 0.5, fatten = 0, alpha = 0.92) +
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    scale_x_discrete(drop = FALSE) +
    labs(title = wraps$title, subtitle = wraps$subtitle,
         x = NULL, y = "% of Parent (batch-normalized)") +
    anti_clip_theme()
  
  fname <- file.path(outdir, paste0(
    sanitize_for_path(glue("{lineage}__{readout}__{condition}__slope")), ".png"))
  ggsave(fname, p, width = w_in, height = wraps$height, units = "in", dpi = 300, bg = "white")
}

comp1_top <- comp1 %>%
  filter(Condition %in% c("BCG","DosR")) %>%
  mutate(Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel)) %>%
  arrange(Condition, desc(abs(median_diff_pp))) %>%
  group_by(Condition) %>% slice_head(n = 25) %>% ungroup()

invisible(pmap(list(comp1_top$Condition, comp1_top$Readout, comp1_top$LineageGate),
               ~ plot_comp1_slope(prolif, ..1, ..2, ..3, dirs$comp1_slopes)))

## ---- Comp2/3 violin plots -------------------------------------------------
plot_violin_pn <- function(df, tp, stim, readout, lineage, igra_col, outdir, label) {
  dat <- df %>%
    filter(PID_Type == "Infant", Time_Point_std == tp, Condition == stim,
           Readout == readout, LineageGate == lineage,
           !is.na(.data[[igra_col]]), .data[[igra_col]] %in% c("Pos","Neg"))
  if (nrow(dat) < 2 || length(unique(dat[[igra_col]])) < 2) return(invisible(NULL))
  
  dat_pid <- dat %>%
    group_by(PID, .data[[igra_col]]) %>%
    summarise(Value_for_test = stats::median(Value_for_test, na.rm = TRUE),
              DenominatorLabel = dplyr::first(DenominatorLabel),
              .groups = "drop")
  
  g_pos <- dat_pid$Value_for_test[dat_pid[[igra_col]] == "Pos"]
  g_neg <- dat_pid$Value_for_test[dat_pid[[igra_col]] == "Neg"]
  g_pos <- g_pos[is.finite(g_pos)]; g_neg <- g_neg[is.finite(g_neg)]
  if (length(g_pos) < 1 || length(g_neg) < 1) return(invisible(NULL))
  
  w <- suppressWarnings(stats::wilcox.test(g_pos, g_neg, paired = FALSE, exact = FALSE))
  pval <- w$p.value
  
  denom_lbl <- unique(dat_pid$DenominatorLabel)[1]
  title_txt <- glue("{mk_readout_label(lineage, readout, denom_lbl)} — {stim} @ {tp}")
  subt_txt  <- glue("Pos n={length(g_pos)}, Neg n={length(g_neg)} | p={signif(pval,3)} {sig_stars(pval)}")
  wraps <- wrap_title_sub(title_txt, subt_txt, 90, 110, 4.8, 0.28)
  
  p <- ggplot(dat_pid, aes(x = .data[[igra_col]], y = Value_for_test, fill = .data[[igra_col]])) +
    geom_violin(width = 0.9, alpha = 0.25, color = NA, trim = TRUE) +
    geom_boxplot(width = 0.46, outlier.shape = NA, alpha = 0.95) +
    geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.7) +
    scale_y_continuous(limits = c(0,100), expand = expansion(mult = c(0.02, 0.06))) +
    coord_cartesian(clip = "off") +
    labs(title = wraps$title, subtitle = wraps$subtitle,
         x = NULL, y = "% of Parent (batch-normalized)") +
    anti_clip_theme()
  
  fname <- file.path(outdir, paste0(
    sanitize_for_path(glue("{lineage}__{readout}__{stim}__{tp}__{label}_violin")), ".png"))
  ggsave(fname, p, width = 6.4, height = wraps$height, units = "in", dpi = 300, bg = "white")
}

comp2_top <- comp2 %>%
  mutate(Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel)) %>%
  arrange(Time_Point_std, Condition, desc(abs(median_diff_pp))) %>%
  group_by(Time_Point_std, Condition) %>% slice_head(n = 25) %>% ungroup()

invisible(pmap(list(comp2_top$Time_Point_std, comp2_top$Condition,
                    comp2_top$Readout, comp2_top$LineageGate),
               ~ plot_violin_pn(prolif, ..1, ..2, ..3, ..4,
                                "Maternal_IGRA", dirs$comp2_violin, "MaternalIGRA")))

comp3_top <- comp3 %>%
  mutate(Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel)) %>%
  arrange(Time_Point_std, Condition, desc(abs(median_diff_pp))) %>%
  group_by(Time_Point_std, Condition) %>% slice_head(n = 25) %>% ungroup()

invisible(pmap(list(comp3_top$Time_Point_std, comp3_top$Condition,
                    comp3_top$Readout, comp3_top$LineageGate),
               ~ plot_violin_pn(prolif, ..1, ..2, ..3, ..4,
                                "Infant_IGRA", dirs$comp3_violin, "InfantIGRA")))

## ---- Comp4 box plots ------------------------------------------------------
plot_comp4_box <- function(df, stim, readout, lineage, outdir) {
  dat <- df %>%
    filter(PID_Type == "Mother", Condition == stim,
           Readout == readout, LineageGate == lineage,
           !is.na(Maternal_IGRA), Maternal_IGRA %in% c("Pos","Neg","ATB"))
  if (nrow(dat) < 3 || length(unique(dat$Maternal_IGRA)) < 2) return(invisible(NULL))
  
  dat_pid <- dat %>%
    group_by(PID, Maternal_IGRA) %>%
    summarise(Value_for_test = stats::median(Value_for_test, na.rm = TRUE),
              DenominatorLabel = dplyr::first(DenominatorLabel),
              .groups = "drop")
  
  kw   <- suppressWarnings(stats::kruskal.test(Value_for_test ~ Maternal_IGRA, data = dat_pid))
  p_kw <- kw$p.value
  
  lvls <- c("Pos","Neg","ATB")
  cnt  <- dat_pid %>% count(Maternal_IGRA) %>%
    right_join(tibble(Maternal_IGRA = lvls), by = "Maternal_IGRA") %>%
    mutate(n = coalesce(n, 0L))
  axis_labels <- setNames(paste0(cnt$Maternal_IGRA, " (n=", cnt$n, ")"), cnt$Maternal_IGRA)
  
  denom_lbl <- unique(dat_pid$DenominatorLabel)[1]
  title_txt <- glue("{mk_readout_label(lineage, readout, denom_lbl)} — {stim}")
  subt_txt  <- glue("Kruskal–Wallis p={signif(p_kw,3)} {sig_stars(p_kw)}")
  wraps <- wrap_title_sub(title_txt, subt_txt, 90, 110, 4.8, 0.28)
  
  p <- ggplot(dat_pid, aes(x = Maternal_IGRA, y = Value_for_test, fill = Maternal_IGRA)) +
    geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.95) +
    geom_jitter(width = 0.13, height = 0, size = 1.6, alpha = 0.7) +
    scale_y_continuous(limits = c(0,100), expand = expansion(mult = c(0.02, 0.06))) +
    scale_x_discrete(limits = lvls, labels = axis_labels, drop = FALSE) +
    coord_cartesian(clip = "off") +
    labs(title = wraps$title, subtitle = wraps$subtitle,
         x = NULL, y = "% of Parent (batch-normalized)") +
    anti_clip_theme()
  
  fname <- file.path(outdir, paste0(
    sanitize_for_path(glue("{lineage}__{readout}__{stim}__Mothers_KW_box")), ".png"))
  ggsave(fname, p, width = 6.8, height = wraps$height, units = "in", dpi = 300, bg = "white")
}

if (nrow(comp4_overall) > 0) {
  comp4_top <- comp4_overall %>%
    filter(!is.na(kw_p)) %>%
    mutate(Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel)) %>%
    arrange(Condition, kw_p) %>%
    group_by(Condition) %>% slice_head(n = 25) %>% ungroup()
  
  invisible(pmap(list(comp4_top$Condition, comp4_top$Readout, comp4_top$LineageGate),
                 ~ plot_comp4_box(prolif, ..1, ..2, ..3, dirs$comp4_box)))
}

## ---- Summary heatmaps -----------------------------------------------------
select_sig_plus_top <- function(df, group_vars, effect_col, p_col, top_n = 30, alpha = 0.05) {
  df  <- df %>% mutate(.is_sig = !!rlang::sym(p_col) < alpha)
  sig <- df %>% filter(.is_sig)
  top <- df %>%
    group_by(across(all_of(group_vars))) %>%
    arrange(desc(abs(!!rlang::sym(effect_col))), .by_group = TRUE) %>%
    slice_head(n = top_n) %>% ungroup()
  bind_rows(sig, top) %>% distinct()
}

compute_heatmap_height <- function(df, facet_vars, row_label_col,
                                   base = 5.5, per_row = 0.18, min_h = 6) {
  if (nrow(df) == 0) return(min_h)
  if (length(facet_vars)) {
    fc <- df %>% group_by(across(all_of(facet_vars))) %>%
      summarise(n_rows = n_distinct(!!rlang::sym(row_label_col)), .groups = "drop")
    max_rows <- if (nrow(fc)) max(fc$n_rows) else n_distinct(df[[row_label_col]])
  } else max_rows <- n_distinct(df[[row_label_col]])
  max(min_h, base + per_row * max_rows)
}

# Comp1 heatmap
if (nrow(comp1) > 0) {
  comp1_hm <- comp1 %>%
    mutate(Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel),
           star = sig_stars(p_value)) %>%
    select_sig_plus_top(group_vars = c("Condition"), effect_col = "median_diff_pp",
                        p_col = "p_value", top_n = 30, alpha = 0.05) %>%
    group_by(Condition) %>% arrange(desc(abs(median_diff_pp)), .by_group = TRUE) %>% ungroup()
  
  if (nrow(comp1_hm) > 0) {
    h_in <- compute_heatmap_height(comp1_hm, "Condition", "Readout_label")
    p1 <- ggplot(comp1_hm,
                 aes(x = "Effect", y = reorder(Readout_label, median_diff_pp), fill = median_diff_pp)) +
      geom_tile(color = "grey70", linewidth = 0.2) +
      geom_text(aes(label = ifelse(star == "", "", star)), size = 3.6, fontface = "bold") +
      scale_fill_gradient2(low = "#1b9e77", mid = "grey92", high = "#d95f02",
                           midpoint = 0, name = "Median Δ(44−12)") +
      facet_wrap(~Condition, scales = "free_y", ncol = 1) +
      labs(title = "Comp1: 12wks → 44wks paired change (Infants, Proliferating)",
           subtitle = "Tiles = median Δ(44−12). Stars = raw p.",
           x = NULL, y = NULL) +
      anti_clip_theme(12) +
      theme(strip.text = element_text(face = "bold"),
            axis.text.y = element_text(size = 8.5),
            legend.position = "right") +
      guides(fill = guide_colorbar(barwidth = 1.5, barheight = 8))
    ggsave(file.path(dirs$summary_heatmaps, "Comp1_Heatmap_medianDiff_p.png"),
           p1, width = 8.0, height = h_in, dpi = 300, bg = "white")
  }
}

# Comp2/3 heatmaps
plot_pn_heatmap <- function(comp_df, title_label, outfile) {
  if (nrow(comp_df) == 0) return(invisible(NULL))
  hm <- comp_df %>%
    mutate(Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel),
           star = sig_stars(p_value)) %>%
    select_sig_plus_top(group_vars = c("Time_Point_std","Condition"),
                        effect_col = "median_diff_pp", p_col = "p_value",
                        top_n = 30, alpha = 0.05) %>%
    group_by(Time_Point_std, Condition) %>%
    arrange(desc(abs(median_diff_pp)), .by_group = TRUE) %>% ungroup() %>%
    mutate(Condition = factor(Condition, levels = c("BCG","DosR","E6C10","GAG","MED")))
  if (nrow(hm) == 0) return(invisible(NULL))
  h_in <- compute_heatmap_height(hm, "Time_Point_std", "Readout_label")
  p <- ggplot(hm, aes(x = Condition, y = reorder(Readout_label, median_diff_pp), fill = median_diff_pp)) +
    geom_tile(color = "grey70", linewidth = 0.2) +
    geom_text(aes(label = ifelse(star == "", "", star)), size = 3.6, fontface = "bold") +
    scale_x_discrete(drop = TRUE) +
    scale_fill_gradient2(low = "#4575b4", mid = "grey92", high = "#d73027",
                         midpoint = 0, name = "Median (Pos−Neg)") +
    facet_grid(Time_Point_std ~ ., scales = "free_y") +
    labs(title = title_label, subtitle = "Tiles = median (Pos−Neg). Stars = raw p.",
         x = NULL, y = NULL) +
    anti_clip_theme(12) +
    theme(axis.text.y = element_text(size = 8.5),
          strip.text.y = element_text(face = "bold", angle = 0),
          legend.position = "right") +
    guides(fill = guide_colorbar(barwidth = 1.5, barheight = 8))
  ggsave(outfile, p, width = 8.4, height = h_in, dpi = 300, bg = "white")
}

plot_pn_heatmap(comp2,
                "Comp2: Maternal IGRA (Pos vs Neg), Infants (Proliferating)",
                file.path(dirs$summary_heatmaps, "Comp2_Heatmap_medianDiff_p.png"))

plot_pn_heatmap(comp3,
                "Comp3: Infant IGRA (Pos vs Neg), Infants (Proliferating)",
                file.path(dirs$summary_heatmaps, "Comp3_Heatmap_medianDiff_p.png"))

# Comp4 heatmap
if (exists("comp4_overall") && nrow(comp4_overall) > 0) {
  hm <- comp4_overall %>%
    mutate(Readout_label = mk_readout_label(LineageGate, Readout, DenominatorLabel),
           star = sig_stars(kw_p)) %>%
    select_sig_plus_top(group_vars = c("Condition"), effect_col = "kw_p",
                        p_col = "kw_p", top_n = 30, alpha = 0.05) %>%
    mutate(Condition = factor(Condition, levels = c("BCG","DosR","GAG","E6C10","MED"))) %>%
    group_by(Condition) %>% arrange(kw_p, .by_group = TRUE) %>% ungroup()
  if (nrow(hm) > 0) {
    h_in <- compute_heatmap_height(hm, character(0), "Readout_label")
    p4 <- ggplot(hm, aes(x = Condition, y = reorder(Readout_label, -kw_p), fill = -log10(kw_p))) +
      geom_tile(color = "grey70", linewidth = 0.2) +
      geom_text(aes(label = ifelse(star == "", "", star)), size = 3.6, fontface = "bold") +
      scale_x_discrete(drop = TRUE) +
      scale_fill_viridis_c(name = expression(-log[10]~p), option = "plasma", na.value = "grey90") +
      labs(title = "Comp4: Mothers (Pos vs Neg vs ATB), Kruskal–Wallis",
           subtitle = "Tiles = −log10(p); stars = raw p.",
           x = NULL, y = NULL) +
      anti_clip_theme(12) +
      theme(axis.text.y = element_text(size = 8.5), legend.position = "right") +
      guides(fill = guide_colorbar(barwidth = 1.5, barheight = 8))
    ggsave(file.path(dirs$summary_heatmaps, "Comp4_Heatmap_KW_p.png"),
           p4, width = 8.6, height = h_in, dpi = 300, bg = "white")
  }
}

message("Script 2 complete. Results in: ", file.path(proj_root, "Results"))