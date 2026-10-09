# Step 07: Supplementary Figs. 4-7 - expression of representative orthogroups.
#   S4  induced by high light in both species (carotenoid / photoprotective)
#   S5  S. polyrhiza-biased (NPQ, anthocyanin/flavonoid)
#   S6  W. australiana-biased (mitochondrial respiratory chain, ribosome)
#   S7  ROS-scavenging enzymes
#
# Panels are addressed by Arabidopsis locus and resolved to the orthogroup through
# Orth_At, never by a hard-coded OrthoID (OrthoFinder renumbers orthogroups every run).
#
# Input : output/logTPM_orthogroup_batch_corrected.Rdata, output/anova_SL_4D_vs_WL_3h.Rdata (03),
#         output/Orth_At.Rdata (02)
# Output: output/SupplementaryFig4_shared_HL_induced_genes.pdf ... output/SupplementaryFig7_Wa_biased_HL_response_genes.pdf
# Run from this directory:  Rscript 07_supfigs.R
suppressPackageStartupMessages({ library(tidyverse); library(ggtext)
                                 library(patchwork); library(ggh4x) })
PT  <- 1 / ggplot2::.pt
OUT <- 'output'

load(file.path(OUT, 'logTPM_orthogroup_batch_corrected.Rdata')); logTPM <- logTPM_orthogroup_batch_corrected
load(file.path(OUT, 'anova_SL_4D_vs_WL_3h.Rdata')); SL4DWL3h_summary <- out_tbl
load(file.path(OUT, 'Orth_At.Rdata')); AtMap <- Orth_At

oid_of <- function(at) {
  o <- AtMap$OrthoID[startsWith(AtMap$At_Tran, paste0(at, '.'))][1]
  if (is.na(o) || !(o %in% logTPM$OrthoID)) NA_character_ else o
}

gene_pair <- function(at, gene_name, wa_ylim = NULL) {
  oid <- oid_of(at)
  if (is.na(oid)) { warning('no orthogroup for ', gene_name, ' (', at, ')'); return(NULL) }
  dat <- logTPM %>% filter(OrthoID == oid) %>%
    pivot_longer(cols = matches('^(Sp|Wa)_(SL_3h|SL_4D|WL_3h)_\\d+$'),
                 names_to = c('Species', 'Condition', 'Rep'),
                 names_pattern = '(Sp|Wa)_(SL_3h|SL_4D|WL_3h)_(\\d+)', values_to = 'Expr') %>%
    mutate(Species = factor(Species, levels = c('Sp', 'Wa')),
           Cond = factor(recode(Condition, WL_3h = 'LL', SL_3h = '3h', SL_4D = '4d'), levels = c('LL', '3h', '4d')))
  md   <- dat %>% group_by(Species, Cond) %>% summarise(M = mean(Expr, na.rm = TRUE), .groups = 'drop')
  padj <- SL4DWL3h_summary$p_adj[SL4DWL3h_summary$OrthoID == oid][1]
  plab <- if (is.na(padj)) 'NS' else paste0('*p* = ', formatC(padj, format = 'e', digits = 1))
  p <- ggplot(dat, aes(Cond, Expr, color = Species)) +
    geom_line(data = md, aes(y = M, group = Species), linewidth = 1.5 * PT) +
    geom_point(data = md, aes(y = M), size = 1.8) +
    geom_jitter(width = 0.12, height = 0, size = 0.9, alpha = 0.5) +
    facet_wrap(~Species, ncol = 2, scales = 'free_y') +
    scale_color_manual(values = c(Sp = '#8375a6', Wa = '#809c4c')) +
    labs(title = paste0('**', gene_name, '**', ' (', plab, ')'), x = NULL, y = 'log2(TPM + 1)') +
    theme_classic(base_size = 9) +
    theme(legend.position = 'none', strip.background = element_blank(),
          strip.text = element_text(face = 'bold', size = 8),
          axis.text.x = element_text(size = 6.5), axis.text.y = element_text(size = 7),
          axis.title.y = element_text(size = 8),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 1.5 * PT),
          axis.line = element_blank(), panel.spacing = grid::unit(0.35, 'lines'),
          plot.title = ggtext::element_markdown(size = 9), plot.margin = margin(4, 4, 4, 4))
  if (!is.null(wa_ylim))
    p <- p + ggh4x::facetted_pos_scales(y = list(Species == 'Wa' ~ scale_y_continuous(limits = wa_ylim)))
  p
}

build <- function(plots, file, nrow_) {
  plots <- plots[!sapply(plots, is.null)]
  fig <- wrap_plots(plots, ncol = 3) + plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(face = 'bold', size = 11, hjust = 0, vjust = 1), plot.tag.position = c(0, 1))
  ggsave(file.path(OUT, file), fig, width = 6.5, height = 1.9 * nrow_ + 0.3, units = 'in')
  cat('wrote', file, '\n')
}

# ---- S4: induced in both species ----
build(list(
  gene_pair('AT3G22840', 'ELIP1'), gene_pair('AT5G17230', 'PSY1'),
  gene_pair('AT3G04870', 'ZDS1'),  gene_pair('AT5G06130', 'ORLIKE'),
  gene_pair('AT3G10230', 'LCYB'),  gene_pair('AT4G25700', 'BCH1')),
  'SupplementaryFig4_shared_HL_induced_genes.pdf', 2)

# ---- S5: S. polyrhiza-biased ----
build(list(
  gene_pair('AT1G44575', 'PSBS/NPQ4'),
  gene_pair('AT4G09820', 'TT8', wa_ylim = c(0, 1)),
  gene_pair('AT5G13930', 'CHS'),   gene_pair('AT1G56500', 'SOQ1'),
  gene_pair('AT1G18400', 'FLAP1'), gene_pair('AT1G61720', 'BAN')),
  'SupplementaryFig5_Sp_biased_HL_response_genes.pdf', 2)

# ---- S6: W. australiana-biased ----
# AT5G25450 is an unnamed Complex III subunit (GO:0045275), labelled by locus.
build(list(
  gene_pair('AT5G67590', 'FRO1'),
  gene_pair('AT2G46540', 'NDUFB10'),
  gene_pair('AT5G25450', 'AT5G25450 (CIII)'),
  gene_pair('AT2G33040', 'ATP3'),
  gene_pair('ATMG00410', 'ATP6'),
  gene_pair('AT1G43170', 'RPL3A')),
  'SupplementaryFig7_Wa_biased_HL_response_genes.pdf', 2)

# ---- S7: ROS scavenging ----
build(list(
  gene_pair('AT4G35090', 'CAT2'), gene_pair('AT2G31570', 'GPX2'),
  gene_pair('AT4G32320', 'APX6')),
  'SupplementaryFig6_ROS_scavenging_genes.pdf', 1)
