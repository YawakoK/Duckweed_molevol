# Step 03: counts -> TPM -> orthogroup-level expression -> batch correction ->
#          within-species limma contrasts -> species x condition ANOVA.
#
# Input : data/featureCounts_Sp.txt, data/featureCounts_Wa.txt   (01_map_and_count.sh)
#         output/Orth_Sp.Rdata, output/Orth_Wa.Rdata              (02)
# Output: output/TPM_orthogroup.Rdata                 orthogroup x sample TPM (sum over paralogs)
#         output/logTPM_orthogroup_batch_corrected.Rdata  log2(TPM+1), limma::removeBatchEffect
#         output/limma_<Sp|Wa>_<A>_vs_<B>.csv/.Rdata   within-species contrasts (3 per species)
#         output/anova_<A>_vs_<B>.csv/.Rdata           per-OG species x condition interaction,
#                                                      BH-adjusted, with per-species mean log2FC
#
# Conditions: WL_3h = LL (low light), SL_3h = HL_3h, SL_4D = HL_4d.
# Batch: replicates 1-3 Azenta (NovaSeq 6000), 4-6 BGI (DNBSEQ-G400).
# Run from this directory:  Rscript 03_tpm_batch_limma_anova.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(tibble)
})
OUT <- 'output'; dir.create(OUT, showWarnings = FALSE)

load(file.path(OUT, 'Orth_Sp.Rdata'))
load(file.path(OUT, 'Orth_Wa.Rdata'))

# ------------------------------------------------------------------
# featureCounts table -> TPM -> orthogroup sums, for one species
# ------------------------------------------------------------------
sample_cols <- function(sp) c(
  paste0(sp, '_SL_3h_', 1:3), paste0(sp, '_SL_4D_', 1:3), paste0(sp, '_WL_3h_', 1:3),
  paste0(sp, '_SL_3h_', 4:6), paste0(sp, '_SL_4D_', 4:6), paste0(sp, '_WL_3h_', 4:6))

read_tpm <- function(sp) {
  d <- read.csv(sprintf('data/featureCounts_%s.txt', sp), sep = '\t', header = FALSE)
  d <- d[c(-1, -2), c(-2, -3, -4, -5)]            # drop header lines and Chr/Start/End/Strand
  col <- c(paste0(sp, '_Tran'), paste0(sp, '_Length'), sample_cols(sp))
  colnames(d) <- col
  TPM <- NULL
  for (i in 3:20) {
    tmp <- (as.numeric(d[, i]) / as.numeric(d[[paste0(sp, '_Length')]])) * 1000
    tmp <- (tmp / sum(tmp)) * 1e6
    TPM <- cbind(TPM, tmp)
  }
  d <- cbind(d[, 1:2], TPM)
  colnames(d) <- col
  d <- cbind(setNames(list(sapply(d[[paste0(sp, '_Tran')]],
                                  function(x) strsplit(x, '_T0')[[1]][1])),
                      paste0(sp, '_Gene')), d)
  d
}

cat('[step 1] Wa counts -> TPM -> orthogroups\n')
Wa <- read_tpm('Wa')
Wa_Orthmerged_sum <- merge(Orth_Wa, Wa) %>%
  group_by(OrthoID) %>%
  summarise(across(matches('^Wa_(SL|WL)_'), \(x) sum(x, na.rm = TRUE))) %>%
  as.data.frame()

cat('[step 2] Sp counts -> TPM -> orthogroups\n')
Sp <- read_tpm('Sp')
Sp_Orthmerged_sum <- merge(Orth_Sp, Sp) %>%
  group_by(OrthoID) %>%
  summarise(across(matches('^Sp_(SL|WL)_'), \(x) sum(x, na.rm = TRUE))) %>%
  as.data.frame()

cat('[step 3] merge Sp + Wa\n')
TPM_orthogroup <- merge(Sp_Orthmerged_sum, Wa_Orthmerged_sum)
save(TPM_orthogroup, file = file.path(OUT, 'TPM_orthogroup.Rdata'))
cat('  orthogroups expressed in both species:', nrow(TPM_orthogroup), '\n')  # -> 10,916 (main text L384)

# ------------------------------------------------------------------
# log2(TPM + 1) and batch correction
# ------------------------------------------------------------------
cat('[step 4] batch correction (limma::removeBatchEffect)\n')
expr_mat <- TPM_orthogroup %>% column_to_rownames('OrthoID') %>% as.matrix()
log_expr_mat <- log2(expr_mat + 1)

sample_info <- data.frame(
  sample    = colnames(log_expr_mat),
  species   = str_extract(colnames(log_expr_mat), '^Sp|Wa'),
  treatment = str_extract(colnames(log_expr_mat), 'SL_3h|SL_4D|WL_3h'),
  replicate = as.integer(str_extract(colnames(log_expr_mat), '[0-9]+$')))
sample_info$batch <- ifelse(sample_info$replicate %in% 1:3, 'batch1', 'batch2')
sample_info$group <- paste(sample_info$species, sample_info$treatment, sep = '_')
group  <- factor(sample_info$group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

log_expr_corrected <- removeBatchEffect(log_expr_mat, batch = sample_info$batch, design = design)
logTPM_orthogroup_batch_corrected <- as.data.frame(log_expr_corrected) %>% rownames_to_column('OrthoID')
save(logTPM_orthogroup_batch_corrected, file = file.path(OUT, 'logTPM_orthogroup_batch_corrected.Rdata'))
write.csv(logTPM_orthogroup_batch_corrected, file.path(OUT, 'logTPM_orthogroup_batch_corrected.csv'),
          row.names = FALSE)

# ------------------------------------------------------------------
# within-species limma contrasts
# ------------------------------------------------------------------
cat('[step 5] within-species limma\n')
run_within_species <- function(logTPM_df, species, condA, condB) {
  cond_regex <- paste(condA, condB, sep = '|')
  cols <- grep(paste0('^', species, '_(', cond_regex, ')_'), names(logTPM_df), value = TRUE)
  stopifnot(length(cols) > 0)
  expr_mat <- logTPM_df %>% dplyr::select(OrthoID, all_of(cols)) %>%
    column_to_rownames('OrthoID') %>% as.matrix()
  sample_tbl <- data.frame(sample = colnames(expr_mat)) %>%
    mutate(Condition = sub(paste0('^', species, '_(', cond_regex, ')_\\d+$'), '\\1', sample),
           Condition = factor(Condition, levels = c(condA, condB)))
  design <- model.matrix(~ 0 + Condition, data = sample_tbl)
  colnames(design) <- levels(sample_tbl$Condition)
  fit  <- limma::lmFit(expr_mat, design)
  fit2 <- limma::contrasts.fit(fit, limma::makeContrasts(contrasts = paste0(condA, '-', condB), levels = design))
  fit2 <- limma::eBayes(fit2)
  res  <- limma::topTable(fit2, number = Inf, sort.by = 'P') %>% rownames_to_column('OrthoID')
  stub <- file.path(OUT, sprintf('limma_%s_%s_vs_%s', species, condA, condB))
  save(res, file = paste0(stub, '.Rdata'))
  readr::write_csv(res, paste0(stub, '.csv'))
  invisible(res)
}
for (sp in c('Sp', 'Wa'))
  for (pair in list(c('SL_4D', 'SL_3h'), c('SL_3h', 'WL_3h'), c('SL_4D', 'WL_3h')))
    run_within_species(logTPM_orthogroup_batch_corrected, sp, pair[1], pair[2])

# ------------------------------------------------------------------
# species x condition ANOVA per orthogroup
# ------------------------------------------------------------------
cat('[step 6] species x condition ANOVA\n')
run_pairwise_anova <- function(logTPM_df, condA, condB) {
  cond_regex <- paste(condA, condB, sep = '|')
  cols <- grep(paste0('^(Sp|Wa)_(', cond_regex, ')_'), names(logTPM_df), value = TRUE)
  expr_mat <- logTPM_df %>% dplyr::select(OrthoID, all_of(cols)) %>%
    column_to_rownames('OrthoID') %>% as.matrix()
  sample_tbl <- data.frame(sample = colnames(expr_mat)) %>%
    mutate(Species   = str_extract(sample, '^Sp|Wa'),
           Condition = sub(paste0('^(Sp|Wa)_(', cond_regex, ')_\\d+$'), '\\2', sample),
           Species   = factor(Species, levels = c('Sp', 'Wa')),
           Condition = factor(Condition, levels = c(condA, condB)))
  anova_results <- apply(expr_mat, 1, function(y) {
    df  <- data.frame(y = y, Species = sample_tbl$Species, Condition = sample_tbl$Condition)
    fit <- tryCatch(aov(y ~ Species * Condition, data = df), error = function(e) NULL)
    if (is.null(fit)) return(c(statistic = NA, p.value = NA))
    s <- summary(fit)[[1]]
    if ('Species:Condition' %in% rownames(s))
      return(c(statistic = s['Species:Condition', 'F value'], p.value = s['Species:Condition', 'Pr(>F)']))
    c(statistic = NA, p.value = NA)
  })
  anova_results <- t(anova_results)
  data.frame(OrthoID = rownames(expr_mat),
             statistic = anova_results[, 'statistic'],
             p.value   = anova_results[, 'p.value']) %>%
    mutate(p_adj = p.adjust(p.value, method = 'BH'))
}

build_summary <- function(logTPM_df, condA, condB) {
  cond_regex <- paste(condA, condB, sep = '|')
  cols <- grep(paste0('^(Sp|Wa)_(', cond_regex, ')_'), names(logTPM_df), value = TRUE)
  long <- logTPM_df %>% dplyr::select(OrthoID, all_of(cols)) %>%
    pivot_longer(-OrthoID, names_to = 'sample', values_to = 'expr') %>%
    mutate(Species   = str_extract(sample, '^Sp|Wa'),
           Condition = sub(paste0('^(Sp|Wa)_(', cond_regex, ')_\\d+$'), '\\2', sample))
  diffs <- long %>%
    group_by(OrthoID, Species, Condition) %>%
    summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Condition, values_from = mean_expr) %>%
    mutate(diff_log = .data[[condA]] - .data[[condB]]) %>%
    dplyr::select(OrthoID, Species, diff_log) %>%
    pivot_wider(names_from = Species, values_from = diff_log, names_glue = '{Species}_diff_log')
  out_tbl <- run_pairwise_anova(logTPM_df, condA, condB) %>% dplyr::left_join(diffs, by = 'OrthoID')
  stub <- file.path(OUT, sprintf('anova_%s_vs_%s', condA, condB))
  save(out_tbl, file = paste0(stub, '.Rdata'))
  readr::write_csv(out_tbl, paste0(stub, '.csv'))
  invisible(out_tbl)
}
SL4D3h_summary   <- build_summary(logTPM_orthogroup_batch_corrected, 'SL_4D', 'SL_3h')
SLWL3h_summary   <- build_summary(logTPM_orthogroup_batch_corrected, 'SL_3h', 'WL_3h')
SL4DWL3h_summary <- build_summary(logTPM_orthogroup_batch_corrected, 'SL_4D', 'WL_3h')
cat('[done]\n')
