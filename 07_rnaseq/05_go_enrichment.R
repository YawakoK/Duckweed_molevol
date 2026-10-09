# Step 05: GO enrichment (topGO, weight01 Fisher) for the response cohorts (Data S8-2).
#
# GO terms are assigned to an orthogroup through its Arabidopsis members
# (TAIR ATH_GO_GOSLIM gene -> GO).  Universe = all orthogroups containing an
# Arabidopsis gene (cohort-based runs) or that universe intersected with the
# trajectory-classified orthogroups (trajectory runs).  Fold enrichment (FE) =
# Significant / Expected is added downstream (08).  Terms with weight01 Fisher
# p < 0.05, >= 10 annotated and >= 3 query orthogroups are reported as enriched
# (Methods); nothing is filtered here.
#
# Input : output/Orth_ath.Rdata, output/Orth_At.Rdata                (02)
#         output/anova_*.Rdata, output/limma_*.csv                   (03)
#         output/transition_pattern_pergene.csv                      (04)
#         data/ATH_GO_GOSLIM_gene2GO.txt   TAIR gene -> GO (2 columns, space separated)
#         data/Araport11_symbols.tsv       transcript -> gene symbol (from Araport11 headers)
# Output: output/GO_<cohort>.csv  (18 files, see `runs` below)
# Run from this directory:  Rscript 05_go_enrichment.R      (~10-20 min)
suppressPackageStartupMessages({
  library(tidyverse)
  library(GO.db)
  library(topGO)
})
OUT <- 'output'

# ---- data ----
load(file.path(OUT, 'Orth_ath.Rdata'))
load(file.path(OUT, 'Orth_At.Rdata'))
At_GO <- read.table('data/ATH_GO_GOSLIM_gene2GO.txt', header = FALSE, sep = ' ') %>% unique()
load(file.path(OUT, 'anova_SL_4D_vs_SL_3h.Rdata')); SL4D3h_summary   <- out_tbl
load(file.path(OUT, 'anova_SL_3h_vs_WL_3h.Rdata')); SLWL3h_summary   <- out_tbl
load(file.path(OUT, 'anova_SL_4D_vs_WL_3h.Rdata')); SL4DWL3h_summary <- out_tbl

# ---- orthogroup -> GO ----
At_GO_df <- unique(At_GO[!is.na(At_GO$V1) & !is.na(At_GO$V2), ])
Orth_At_gene <- Orth_At %>%
  dplyr::mutate(At_gene = sub('\\..*$', '', At_Tran)) %>%
  dplyr::distinct(OrthoID, At_gene)
orth_go_tbl <- Orth_At_gene %>%
  dplyr::inner_join(At_GO_df, by = c('At_gene' = 'V1')) %>%
  dplyr::transmute(OrthoID, GO = V2) %>%
  dplyr::distinct()
geneGO <- lapply(split(orth_go_tbl$GO, orth_go_tbl$OrthoID), unique)
geneNames_all <- unique(Orth_At_gene$OrthoID)
geneGO_full <- setNames(rep(list(character(0)), length(geneNames_all)), geneNames_all)
geneGO_full[names(geneGO)] <- geneGO
geneGO <- geneGO_full

GO_terms <- Term(ls(GOTERM)) %>% as.data.frame()
GO_terms <- data.frame('GO.ID' = rownames(GO_terms), term = GO_terms[, 1])

geneNames <- names(geneGO)
all_genes <- geneGO

# ---- annotation helpers (orthogroup id, TAIR loci and symbols per term) ----
build_tair_symbol_map <- function(tsv_path) {
  stopifnot(file.exists(tsv_path))
  d <- read.delim(tsv_path, header = FALSE, col.names = c('gene_id', 'SYMBOL'),
                  quote = '', stringsAsFactors = FALSE)
  d$SYMBOL <- trimws(d$SYMBOL); d$SYMBOL[d$SYMBOL == ''] <- NA_character_
  tibble::tibble(At_gene = sub('\\..*$', '', d$gene_id), SYMBOL = d$SYMBOL) %>%
    dplyr::group_by(At_gene) %>%
    dplyr::summarize(SYMBOL = dplyr::first(SYMBOL[!is.na(SYMBOL)]), .groups = 'drop')
}
collapse_unique <- function(x) {
  x <- unique(x); x <- x[!is.na(x) & x != '']
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = ';')
}
add_term_details_topgo <- function(go_table, go_data, focus_genes, dir_map,
                                   orth_to_og, orth_to_tair, tair_to_symbol) {
  gene_sets <- topGO::genesInTerm(go_data, go_table$GO.ID)
  og_col <- character(nrow(go_table)); tair_col <- character(nrow(go_table))
  sym_col <- character(nrow(go_table))
  use_dir <- !is.null(dir_map)
  if (use_dir) {
    sp_counts <- integer(nrow(go_table)); wa_counts <- integer(nrow(go_table))
    eq_counts <- integer(nrow(go_table))
  }
  for (i in seq_len(nrow(go_table))) {
    genes_in_term <- gene_sets[[i]]
    if (length(genes_in_term) > 0) {
      genes_in_term <- unique(genes_in_term)
      if (!is.null(focus_genes)) genes_in_term <- intersect(genes_in_term, focus_genes)
    }
    if (length(genes_in_term) > 0) {
      og_col[i] <- collapse_unique(orth_to_og[genes_in_term])
      tair_vec <- unique(unlist(orth_to_tair[genes_in_term], use.names = FALSE))
      tair_col[i] <- collapse_unique(tair_vec)
      sym_col[i] <- collapse_unique(tair_to_symbol[tair_vec])
      if (use_dir) {
        dirs <- dir_map[genes_in_term]
        sp_counts[i] <- sum(dirs == 'Sp_gt_Wa', na.rm = TRUE)
        wa_counts[i] <- sum(dirs == 'Wa_gt_Sp', na.rm = TRUE)
        eq_counts[i] <- sum(dirs == 'Sp_eq_Wa', na.rm = TRUE)
      }
    } else {
      og_col[i] <- NA_character_; tair_col[i] <- NA_character_; sym_col[i] <- NA_character_
      if (use_dir) { sp_counts[i] <- 0L; wa_counts[i] <- 0L; eq_counts[i] <- 0L }
    }
  }
  if (use_dir) {
    go_table$Sp_gt_Wa <- sp_counts; go_table$Wa_gt_Sp <- wa_counts; go_table$Sp_eq_Wa <- eq_counts
  }
  go_table$OG <- og_col; go_table$tairID <- tair_col; go_table$SYMBOL <- sym_col
  go_table
}

orth_map   <- Orth_ath %>% dplyr::select(OrthoID, Orthogroup) %>% dplyr::distinct()
orth_to_og <- setNames(orth_map$Orthogroup, orth_map$OrthoID)
orth_to_tair <- Orth_At_gene %>% dplyr::distinct(OrthoID, At_gene)
orth_to_tair <- lapply(split(orth_to_tair$At_gene, orth_to_tair$OrthoID), unique)
tair_symbol_tbl <- build_tair_symbol_map('data/Araport11_symbols.tsv')
tair_to_symbol  <- setNames(tair_symbol_tbl$SYMBOL, tair_symbol_tbl$At_gene)

# ---- topGO wrapper ----
makeGOtable <- function(A, B, focus_genes = NULL, dir_map = NULL) {
  one_ont <- function(ont) {
    GOdata <- new('topGOdata', ontology = ont, allGenes = A, annot = annFUN.gene2GO, gene2GO = B)
    r_cl <- runTest(GOdata, algorithm = 'classic',  statistic = 'fisher')
    r_w  <- runTest(GOdata, algorithm = 'weight01', statistic = 'fisher')
    len  <- score(r_cl) %>% length()
    tb <- GenTable(GOdata, Fisher = r_cl, orderBy = 'Fisherweight01', Fisherweight01 = r_w, topNodes = len)
    tb <- cbind(tb, ontology = ont)
    tb <- cbind(tb, Fisherweight01FDR = p.adjust(tb$Fisherweight01, 'BH'))
    add_term_details_topgo(tb, GOdata, focus_genes, dir_map, orth_to_og, orth_to_tair, tair_to_symbol)
  }
  result.table <- rbind(one_ont('MF'), one_ont('CC'), one_ont('BP'))
  result.table <- dplyr::select(result.table, -Term)
  result.table <- merge(result.table, GO_terms, by = 'GO.ID')
  result.table <- result.table %>%
    dplyr::rename(Term = term, Ontology = ontology) %>%
    dplyr::select(GO.ID, Term, Ontology, Fisherweight01, Annotated, Significant, Expected, Fisher,
                  dplyr::any_of(c('Sp_gt_Wa', 'Wa_gt_Sp', 'Sp_eq_Wa', 'OG', 'tairID', 'SYMBOL')))
  result.table %>%
    dplyr::mutate(Fisherweight01_num = suppressWarnings(as.numeric(sub('^<\\s*', '', Fisherweight01)))) %>%
    dplyr::arrange(Fisherweight01_num) %>%
    dplyr::select(-Fisherweight01_num)
}
make_dir_map <- function(pair_summary, focus_ortho) {
  dir_tbl <- pair_summary %>%
    dplyr::mutate(direction = dplyr::case_when(
      Sp_diff_log > Wa_diff_log ~ 'Sp_gt_Wa',
      Sp_diff_log < Wa_diff_log ~ 'Wa_gt_Sp',
      TRUE ~ 'Sp_eq_Wa')) %>%
    dplyr::select(OrthoID, direction) %>%
    dplyr::filter(OrthoID %in% focus_ortho)
  setNames(dir_tbl$direction, dir_tbl$OrthoID)
}
run_go <- function(focus_ortho, name, dir_map = NULL) {
  Focus_Genes <- intersect(unique(focus_ortho), geneNames)
  geneList <- factor(as.integer(geneNames %in% Focus_Genes))
  names(geneList) <- geneNames
  res <- makeGOtable(geneList, all_genes, focus_genes = Focus_Genes, dir_map = dir_map)
  out_csv <- file.path(OUT, paste0('GO_', name, '.csv'))
  readr::write_csv(res, out_csv)
  cat('  wrote', out_csv, '(', nrow(res), 'GO terms,', length(Focus_Genes), 'focus OGs )\n')
  res
}

# =====================================================================
# A. species x condition cohorts (ANOVA interaction, BH p < 0.05)
#    HL4d_vs_LL = SL4DWL3h, HL3h_vs_LL = SLWL3h, HL4d_vs_HL3h = SL4D3h
# =====================================================================
for (cmp in list(list(s = SL4DWL3h_summary, n = 'SL4DWL3h'),
                 list(s = SL4D3h_summary,   n = 'SL4D3h'),
                 list(s = SLWL3h_summary,   n = 'SLWL3h'))) {
  cat('[GO]', cmp$n, '...\n')
  tmp <- cmp$s
  Focus <- tmp$OrthoID[tmp$p_adj < 0.05]
  run_go(Focus, paste0(cmp$n, '_padj005'), make_dir_map(tmp, Focus))
  run_go(tmp$OrthoID[tmp$p_adj < 0.05 & tmp$Sp_diff_log > tmp$Wa_diff_log], paste0(cmp$n, '_padj005_Sp_gt_Wa'))
  run_go(tmp$OrthoID[tmp$p_adj < 0.05 & tmp$Sp_diff_log < tmp$Wa_diff_log], paste0(cmp$n, '_padj005_Wa_gt_Sp'))
}

# =====================================================================
# B. induced in both species (within-species limma adj.P < 0.05 & logFC > 0)
# =====================================================================
lim <- function(sp, a, b) read.csv(file.path(OUT, sprintf('limma_%s_%s_vs_%s.csv', sp, a, b)))
cat('[GO] both-species induced ...\n')
for (pr in list(c('SL_4D', 'SL_3h', 'SL4D_gt_SL3h'), c('SL_4D', 'WL_3h', 'SL4D_gt_WL3h'), c('SL_3h', 'WL_3h', 'SL3h_gt_WL3h'))) {
  s <- lim('Sp', pr[1], pr[2]); w <- lim('Wa', pr[1], pr[2])
  Focus <- intersect(s$OrthoID[s$adj.P.Val < 0.05 & s$logFC > 0],
                     w$OrthoID[w$adj.P.Val < 0.05 & w$logFC > 0])      # HL4d vs LL: 204 (main text L397)
  run_go(Focus, paste0(pr[3], '_bothSpecies'))
}

# =====================================================================
# C. trajectory cohorts (Fig. 5A patterns).  Universe = trajectory-classified
#    orthogroups that carry a GO annotation.
# =====================================================================
pp <- read.csv(file.path(OUT, 'transition_pattern_pergene.csv'), stringsAsFactors = FALSE)
univ <- intersect(unique(pp$OrthoID), names(geneGO))
cat(sprintf('[pattern-GO] pattern OGs=%d ; annotated background universe=%d\n', nrow(pp), length(univ)))
geneNames <- univ
all_genes <- geneGO[univ]
sets <- list(
  Sp_UF = pp$OrthoID[pp$Sp_pat == 'UF'],
  Sp_FU = pp$OrthoID[pp$Sp_pat == 'FU'],
  Sp_DF = pp$OrthoID[pp$Sp_pat == 'DF'],
  Wa_FU = pp$OrthoID[pp$Wa_pat == 'FU'],
  Wa_FD = pp$OrthoID[pp$Wa_pat == 'FD'],
  Wa_DF = pp$OrthoID[pp$Wa_pat == 'DF'])
for (nm in names(sets)) {
  cat(sprintf('[GO] pattern set %s (n=%d) ...\n', nm, length(unique(sets[[nm]]))))
  run_go(sets[[nm]], paste0('pattern_', nm))
}
cat('[done]\n')
