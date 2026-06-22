===============================================================================
 BCG Parse — IGRA+ vs IGRA- DGE :  README / file guide
 Output root: DGE/IGRA_pos_vs_neg_7sets/
 Script:      BCG_PARSE_DGE_IGRA_pseudobulk_MAST.R
===============================================================================

-------------------------------------------------------------------------------
 1. WHAT THIS ANALYSIS ASKS
-------------------------------------------------------------------------------
For each cell-type set, at each timepoint (week12 and week44, analysed
SEPARATELY), which genes differ between IGRA-positive and IGRA-negative infants
-- and specifically, which of those differences are part of the BCG RESPONSE
rather than a baseline difference already present without stimulation.

IGRA status is a DONOR-level property (4 IGRA- donors, 6 IGRA+ donors), so the
biological unit of replication is the donor, not the cell. That single fact
drives every methods choice below.

Direction convention everywhere: POSITIVE fold-change = HIGHER in IGRA+.

-------------------------------------------------------------------------------
 2. CELL-TYPE SETS (12)
-------------------------------------------------------------------------------
Some sets pool more than one annotation label (the "+" units). Pooling is done
by label, not by re-clustering.

  CD4_Naive_T      = CD4 Naive T
  CD4_TCM          = CD4 TCM
  CD8_Naive_T      = CD8 Naive T
  CD8_TCM_TEM      = CD8 TCM + CD8 TEM
  Treg             = Treg
  gdT              = gdT (Vd1) + gdT (cytotoxic) + gdT (Vg9Vd2)
  CD14_Mono_Infl   = CD14 Monocyte + CD14 Monocyte (inflammatory)
  MAIT             = MAIT
  NK               = NK + NK CD56bright + NK (activated)
  B                = B naive + B intermediate + B memory
  Plasmablast      = Plasmablast
  DC               = pDC + Activated DC

Excluded: doublet / low-quality clusters and CD4 T (activated/stress), a
QC-flagged population.

-------------------------------------------------------------------------------
 3. FOLDER MAP
-------------------------------------------------------------------------------
  Pseudobulk_DESeq2/      PRIMARY result. Donor-level DESeq2.
  BCG_response_specific/  The BCG-vs-baseline filtering applied to the above.
  MAST_sensitivity/       SECONDARY cell-level check + pseudoreplication flags.
  Plots/                  Heatmaps and violins of the BCG-specific genes.

===============================================================================
 4. Pseudobulk_DESeq2/   --  the primary, reportable result
===============================================================================
WHAT PSEUDOBULK MEANS
  For each donor, within one stratum (e.g. CD4_TCM cells that are BCG-treated at
  week12), the RAW UMI COUNTS of all that donor's cells are SUMMED into a single
  profile -- collapsing thousands of cells back into one mini bulk-RNA-seq
  sample per donor. The IGRA comparison is then 6 IGRA+ profiles vs 4 IGRA-
  profiles, which honours the true sample size. (Cell-level tests instead treat
  every cell as independent and so massively overstate N -- see section 6.)

HOW THE TEST WORKS
  DESeq2 fits, per gene, the model  ~ IGRA_status  and tests Positive vs
  Negative. This is run TWICE per (set x timepoint): once in BCG-treated cells
  and once in media-treated cells -- the two runs feed section 5.

FILES
  pseudobulk_sample_table.csv
      One row per donor-pseudobulk: set, sample_id, stim, timepoint, n_cells
      (how many cells were summed), IGRA_status. Donors with n_cells below the
      threshold are dropped before testing.
  pseudobulk_checkpoint.qs2
      The pseudobulk count matrix + sample table (re-loadable; not for sharing).
  PB_<set>_<timepoint>_<stim>.csv
      DESeq2 results for one stratum. Columns:
        gene            gene symbol
        baseMean        mean normalised expression across the donors
        log2FoldChange  IGRA+ vs IGRA-  (positive = up in IGRA+)
        lfcSE           standard error of the log2FC
        stat            Wald statistic
        pvalue          raw p-value
        padj            BH-adjusted p-value  <-- use this for significance
        set, timepoint, stim
        n_pos, n_neg    number of donors per IGRA arm that contributed
                        (CHECK THESE: a result from 2 vs 2 donors is fragile)
  ALL__pseudobulk_DESeq2_LONG.csv
      All PB_* rows concatenated. (Large; gitignored -- regenerate by re-running.)

===============================================================================
 5. BCG_response_specific/   --  isolating the BCG-driven IGRA signal
===============================================================================
WHY
  A gene that differs by IGRA in UNSTIMULATED (media) cells is a baseline
  difference, not something BCG did. We only want genes whose IGRA difference is
  specific to (or changed by) BCG stimulation.

HOW IT COMPARES
  Per (set x timepoint), the BCG DESeq2 result and the media DESeq2 result are
  joined gene-by-gene, and each gene is classified:
        sig_bcg      significant by IGRA in BCG    (padj<0.05 AND |log2FC|>=0.585)
        sig_media    significant by IGRA in media  (same thresholds)
        opp_sign     BCG and media effects point in opposite directions
        bcg_specific = sig_bcg AND ( NOT sig_media  OR  opp_sign )
  So a "BCG-response-specific" gene is one that moves with IGRA under BCG but
  not under media, OR flips direction between the two.

  NOTE (important caveat): this is a difference-of-significance rule -- it can
  call a gene "specific" when the BCG and media effects are actually similar in
  size and only differ in which side of p=0.05 they fell. The stricter test of
  the same idea is a stim x IGRA INTERACTION model; if a Section 2b interaction
  table is present, prefer it for anything you intend to claim.

FILES
  SPEC_<set>_<timepoint>.csv
      Every gene with both BCG and media stats, sorted with bcg_specific first.
      Columns: gene, lfc_bcg, padj_bcg, lfc_med, padj_med, sig_bcg, sig_media,
               opp_sign, bcg_specific, set, timepoint.
      (If a set had no media result, lfc_med/padj_med are NA and specificity vs
       media is unverified for that set -- noted in the run log.)
  ALL__BCG_response_specific_LONG.csv
      The bcg_specific == TRUE genes across all sets/timepoints, concatenated.

===============================================================================
 6. MAST_sensitivity/   --  secondary check, NOT the primary result
===============================================================================
WHAT IT IS
  MAST is the cell-level test (every cell = one data point), run on BCG-treated
  cells only, IGRA+ vs IGRA-, per (set x timepoint). It is included as a
  SENSITIVITY analysis -- a second method to see whether the pseudobulk picture
  holds up -- not as the headline result.

WHY SECONDARY
  Because MAST treats each cell as independent, it effectively claims an N in
  the thousands for a comparison that biologically has N=10 donors. This
  "pseudoreplication" produces long, over-confident hit lists, many driven by
  one or two donors. Two columns expose this directly (below).

FILES
  MAST_<set>_<timepoint>_BCG.csv
      Standard FindMarkers(MAST) output plus per-donor diagnostics. Columns:
        gene
        p_val               raw p-value
        avg_log2FC          IGRA+ vs IGRA-  (positive = up in IGRA+)
        pct.1               fraction of IGRA+ cells expressing the gene
        pct.2               fraction of IGRA- cells expressing the gene
        p_val_adj           adjusted p-value
        n_donors_pos_detect # of IGRA+ DONORS expressing it  (computed for hits)
        n_donors_neg_detect # of IGRA- DONORS expressing it
        one_donor_flag      TRUE if min(pos,neg) donors <= 1
                            --> a "hit" carried by a single donor; do not trust
        set, timepoint
  OVERLAP_pseudobulk_vs_MAST.csv
      Per (set x timepoint): n_pseudobulk_sig, n_MAST_sig, n_shared.
      A large gap (e.g. MAST 300 / pseudobulk 15) is the inflation made visible.

HOW TO READ THE TWO TOGETHER
  Trust genes that PSEUDOBULK calls and MAST corroborates. MAST-only hits are
  exploratory: check one_donor_flag and the n_donors columns before believing
  any of them. MAST earns its keep mainly for rare sets where pseudobulk runs
  out of donors.

===============================================================================
 7. Plots/
===============================================================================
  Heatmap_<set>_<timepoint>.png
      Z-scored mean expression of the BCG-specific genes across the FOUR
      IGRA x stim groups (IGRA+ BCG, IGRA- BCG, IGRA+ media, IGRA- media).
      A real BCG-specific gene should separate by IGRA in the BCG columns but
      NOT in the media columns -- the heatmap lets you eyeball exactly that.
  Vln_<set>_<timepoint>.png
      Violin plots of the top BCG-specific genes in BCG cells only, IGRA+ vs
      IGRA-, with a Wilcoxon annotation.

===============================================================================
 8. THRESHOLDS USED
===============================================================================
  Significance (pseudobulk & specificity):  padj < 0.05  AND  |log2FC| >= 0.585 (1.5x)
  MAST significance:                         p_val_adj < 0.05
  Min cells per donor-pseudobulk:            10   (else that donor is dropped)
  Min donors per IGRA arm to run DESeq2:     2    (3+ is much safer; watch n_pos/n_neg)
  Min cells per IGRA arm for MAST:           30

===============================================================================
 9. ONE-LINE SUMMARY OF THE LOGIC
===============================================================================
  Pseudobulk DESeq2 (donor-level)  ->  what you can defensibly claim.
  BCG-vs-media filter              ->  keep only BCG-driven IGRA differences.
  MAST + per-donor flags           ->  cross-check; expose one-donor artefacts.
===============================================================================
