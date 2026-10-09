# Step 08: supplementary data tables.
#   Data S8-2  all GO enrichment results, unfiltered, with fold enrichment (FE)  (05)
#   Data S8-1  per-orthogroup high-light response: log2FC, ANOVA p, trajectory,
#              Arabidopsis loci, HL4d category                                 (03, 04)
#   Data S0-2  orthogroup membership for all eight species                     (02)
#
# Input : output/GO_*.csv, output/anova_*.Rdata, output/limma_*_SL_4D_vs_WL_3h.csv,
#         output/transition_pattern_pergene.csv, output/Orth_At.Rdata, output/Orth_ath.Rdata
# Output: output/SupplementaryData20_GO_enrichment_all.csv, output/SupplementaryData19_per_orthogroup_response.csv,
#         output/SupplementaryData18_orthogroups.csv
# Run from this directory:  Rscript 08_suptables.R
suppressPackageStartupMessages({ library(tidyverse) })
OUT <- 'output'

# ---------------- Data S8-2: all GO enrichment results ----------------
files <- list(
  list(f = 'GO_SL4D_gt_WL3h_bothSpecies',  comp = 'HL4d_vs_LL',   coh = 'Both_up'),
  list(f = 'GO_SL3h_gt_WL3h_bothSpecies',  comp = 'HL3h_vs_LL',   coh = 'Both_up'),
  list(f = 'GO_SL4D_gt_SL3h_bothSpecies',  comp = 'HL4d_vs_HL3h', coh = 'Both_up'),
  list(f = 'GO_SL4DWL3h_padj005_Sp_gt_Wa', comp = 'HL4d_vs_LL',   coh = 'Sp_biased'),
  list(f = 'GO_SL4DWL3h_padj005_Wa_gt_Sp', comp = 'HL4d_vs_LL',   coh = 'Wa_biased'),
  list(f = 'GO_SLWL3h_padj005_Sp_gt_Wa',   comp = 'HL3h_vs_LL',   coh = 'Sp_biased'),
  list(f = 'GO_SLWL3h_padj005_Wa_gt_Sp',   comp = 'HL3h_vs_LL',   coh = 'Wa_biased'),
  list(f = 'GO_SL4D3h_padj005_Sp_gt_Wa',   comp = 'HL4d_vs_HL3h', coh = 'Sp_biased'),
  list(f = 'GO_SL4D3h_padj005_Wa_gt_Sp',   comp = 'HL4d_vs_HL3h', coh = 'Wa_biased'),
  list(f = 'GO_pattern_Sp_UF', comp = 'trajectory', coh = 'Sp_UF'),
  list(f = 'GO_pattern_Sp_FU', comp = 'trajectory', coh = 'Sp_FU'),
  list(f = 'GO_pattern_Sp_DF', comp = 'trajectory', coh = 'Sp_DF'),
  list(f = 'GO_pattern_Wa_FU', comp = 'trajectory', coh = 'Wa_FU'),
  list(f = 'GO_pattern_Wa_FD', comp = 'trajectory', coh = 'Wa_FD'),
  list(f = 'GO_pattern_Wa_DF', comp = 'trajectory', coh = 'Wa_DF'))

S1 <- purrr::map_dfr(files, function(x) {
  p <- file.path(OUT, paste0(x$f, '.csv'))
  if (!file.exists(p)) { warning('missing: ', p); return(NULL) }
  d <- read.csv(p)
  d %>% mutate(
      Fw  = suppressWarnings(as.numeric(sub('^< *', '', Fisherweight01))),
      Sig = as.numeric(Significant), Exp = as.numeric(Expected),
      Ann = as.numeric(Annotated),
      FE  = ifelse(Exp > 0, round(Sig / Exp, 2), NA_real_)) %>%   # FE undefined when Expected rounds to 0
    transmute(comparison = x$comp, cohort = x$coh,
              GO.ID, Term, Ontology,
              Fisher_weight01 = Fw, Annotated = Ann,
              Significant = Sig, Expected = Exp, FE,
              SYMBOL = if ('SYMBOL' %in% colnames(d)) SYMBOL else NA_character_)
})
write.csv(S1, file.path(OUT, 'SupplementaryData20_GO_enrichment_all.csv'), row.names = FALSE, na = '')
cat('Data S8-2 rows:', nrow(S1), '\n')

# The GO terms quoted in the main text (L398-431) can be looked up here, e.g.
print(S1 %>% filter(comparison == 'HL4d_vs_LL', cohort == 'Both_up',
                    GO.ID %in% c('GO:0016117', 'GO:0009813', 'GO:0010287', 'GO:0009535')) %>%
        select(GO.ID, Term, Fisher_weight01, FE))

# ---------------- Data S8-1: per-orthogroup response ----------------
load(file.path(OUT, 'anova_SL_4D_vs_WL_3h.Rdata')); S4 <- out_tbl
load(file.path(OUT, 'anova_SL_3h_vs_WL_3h.Rdata')); S3 <- out_tbl
load(file.path(OUT, 'Orth_At.Rdata'))
pp  <- read.csv(file.path(OUT, 'transition_pattern_pergene.csv'), stringsAsFactors = FALSE)
sp4 <- read.csv(file.path(OUT, 'limma_Sp_SL_4D_vs_WL_3h.csv'))
wa4 <- read.csv(file.path(OUT, 'limma_Wa_SL_4D_vs_WL_3h.csv'))

both4 <- intersect(sp4$OrthoID[sp4$adj.P.Val < 0.05 & sp4$logFC > 0],
                   wa4$OrthoID[wa4$adj.P.Val < 0.05 & wa4$logFC > 0])
at_map <- Orth_At %>%
  mutate(At_gene = sub('\\..*$', '', At_Tran)) %>%
  group_by(OrthoID) %>%
  summarise(At_genes = paste(unique(At_gene), collapse = ';'), .groups = 'drop')

S2 <- S4 %>%
  transmute(OrthoID, Sp_diff_HL4d = Sp_diff_log, Wa_diff_HL4d = Wa_diff_log, p_adj_HL4d = p_adj) %>%
  left_join(S3 %>% transmute(OrthoID, Sp_diff_HL3h = Sp_diff_log, Wa_diff_HL3h = Wa_diff_log, p_adj_HL3h = p_adj),
            by = 'OrthoID') %>%
  left_join(pp %>% transmute(OrthoID, Sp_pattern = Sp_pat, Wa_pattern = Wa_pat), by = 'OrthoID') %>%
  left_join(at_map, by = 'OrthoID') %>%
  mutate(Category_HL4d = case_when(
    OrthoID %in% both4                        ~ 'Both_up',
    is.na(p_adj_HL4d) | p_adj_HL4d >= 0.05    ~ 'NS',
    Sp_diff_HL4d > Wa_diff_HL4d               ~ 'Sp_biased',
    TRUE                                      ~ 'Wa_biased'))
write.csv(S2, file.path(OUT, 'SupplementaryData19_per_orthogroup_response.csv'), row.names = FALSE, na = '')
cat('Data S8-1 rows:', nrow(S2), '\n')
print(S2 %>% count(Category_HL4d))            # -> Both_up 204, Sp_biased 1792, Wa_biased 435

# ---------------- Data S0-2: orthogroup membership ----------------
load(file.path(OUT, 'Orth_ath.Rdata'))
ren <- c(
  Araport11_pep_20250214 = 'Arabidopsis_thaliana',
  Colocasia_esculenta    = 'Colocasia_esculenta',
  Sp_withorga            = 'Spirodela_polyrhiza',
  La_longest.isoform     = 'Landoltia_punctata',
  Laeq_longest.isoform   = 'Lemna_aequinoctialis',
  Lgib                   = 'Lemna_gibba',
  Wo_longest.isoform     = 'Wolffiella_hyalina',
  Wa_withorga            = 'Wolffia_australiana')
d <- Orth_ath
for (old in names(ren)) if (old %in% colnames(d)) colnames(d)[colnames(d) == old] <- ren[[old]]
spcols <- unname(ren)
n_of <- function(x) { x[is.na(x)] <- ''; ifelse(x == '', 0L, lengths(strsplit(x, ', '))) }
out <- d %>%
  select(OrthoID, Orthogroup, all_of(spcols)) %>%
  mutate(across(all_of(spcols), ~ ifelse(is.na(.x), '', .x)))
counts <- as.matrix(sapply(spcols, function(c) n_of(out[[c]]))); colnames(counts) <- spcols
cnt_df <- as.data.frame(counts); colnames(cnt_df) <- paste0('n_', spcols)
out <- bind_cols(out, cnt_df) %>%
  mutate(n_total = rowSums(counts), n_species = rowSums(counts > 0))
write.csv(out, file.path(OUT, 'SupplementaryData18_orthogroups.csv'), row.names = FALSE, na = '')
cat('Data S0-2 rows:', nrow(out), '\n')
cat('single-copy in the 7 non-Arabidopsis taxa:',
    sum(rowSums(counts[, setdiff(spcols, 'Arabidopsis_thaliana')] == 1) == 7), '\n')   # -> 2,819
