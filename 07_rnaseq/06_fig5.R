# Step 06: Fig. 5 (A: trajectory frequencies; B, C: Sp vs Wa response scatter;
#          D-I: six representative orthogroups).
#
# Input : output/logTPM_orthogroup_batch_corrected.Rdata, output/anova_*.Rdata,
#         output/limma_*_vs_WL_3h.csv (03); output/transition_pattern_pergene.csv (04);
#         output/Orth_At.Rdata (02) - used only to verify the six panel orthogroups.
# Output: output/Fig5.pdf
# Run from this directory:  Rscript 06_fig5.R
suppressPackageStartupMessages({ library(tidyverse); library(ggtext)
                                 library(patchwork); library(stringr) })
PT  <- 1 / ggplot2::.pt
OUT <- 'output'

# ---- data ----
load(file.path(OUT, 'logTPM_orthogroup_batch_corrected.Rdata')); logTPM <- logTPM_orthogroup_batch_corrected
load(file.path(OUT, 'anova_SL_3h_vs_WL_3h.Rdata')); SLWL3h_summary   <- out_tbl
load(file.path(OUT, 'anova_SL_4D_vs_WL_3h.Rdata')); SL4DWL3h_summary <- out_tbl
load(file.path(OUT, 'Orth_At.Rdata'))
sp3 <- read.csv(file.path(OUT, 'limma_Sp_SL_3h_vs_WL_3h.csv'))
wa3 <- read.csv(file.path(OUT, 'limma_Wa_SL_3h_vs_WL_3h.csv'))
sp4 <- read.csv(file.path(OUT, 'limma_Sp_SL_4D_vs_WL_3h.csv'))
wa4 <- read.csv(file.path(OUT, 'limma_Wa_SL_4D_vs_WL_3h.csv'))
pp  <- read.csv(file.path(OUT, 'transition_pattern_pergene.csv'), stringsAsFactors = FALSE)

# ---- A : trajectory frequencies with trajectory icons (FF excluded) ----
lev <- c('UU', 'UF', 'UD', 'FU', 'FD', 'DU', 'DF', 'DD')
N   <- nrow(pp)
long <- bind_rows(tibble(pattern = pp$Sp_pat, Species = 'Sp'),
                  tibble(pattern = pp$Wa_pat, Species = 'Wa')) %>%
  filter(pattern %in% lev) %>% count(Species, pattern, name = 'n') %>%
  mutate(pattern = factor(pattern, levels = lev),
         Species = factor(Species, levels = c('Sp', 'Wa')))
ymax <- max(long$n); ymid <- -0.25 * ymax; amp <- 0.11 * ymax
dlt <- c(U = 1, F = 0, D = -1)
sch <- tibble(pattern = lev, i = seq_along(lev),
              t1 = substr(lev, 1, 1), t2 = substr(lev, 2, 2)) %>%
  mutate(y1 = 0, y2 = dlt[t1], y3 = dlt[t1] + dlt[t2]) %>%
  pivot_longer(c(y1, y2, y3), names_to = 'tp', values_to = 'traj') %>%
  mutate(x = i + recode(tp, y1 = -0.24, y2 = 0, y3 = 0.24),
         y = ymid + (traj / 2) * amp)
dg <- position_dodge(width = 0.7)
A <- ggplot(long, aes(pattern, n, fill = Species)) +
  geom_col(position = dg, width = 0.62) +
  geom_text(aes(label = n), position = dg, vjust = -0.35, size = 7 * PT, color = 'gray20') +
  geom_hline(yintercept = ymid, linetype = 'dotted', linewidth = 0.5 * PT, color = 'gray70') +
  geom_line(data = sch, aes(x = x, y = y, group = pattern), inherit.aes = FALSE,
            linewidth = 2 * PT, color = 'gray15',
            arrow = arrow(length = grid::unit(2.4, 'pt'), type = 'closed')) +
  geom_point(data = sch, aes(x = x, y = y), inherit.aes = FALSE, size = 0.9, color = 'gray15') +
  scale_fill_manual(values = c(Sp = '#8375a6', Wa = '#809c4c'),
                    labels = c(Sp = '*Spirodela polyrhiza*', Wa = '*Wolffia australiana*')) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  scale_y_continuous(breaks = c(0, 200, 400, 600),
                     limits = c(ymid - amp * 1.55, ymax * 1.16), expand = expansion(add = 0)) +
  labs(x = NULL, y = paste0('Orthogroups (of ', format(N, big.mark = ','), ')'), fill = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.background = element_rect(fill = scales::alpha('white', 0.85), color = NA),
        legend.key.size = grid::unit(3.5, 'mm'),
        legend.text = ggtext::element_markdown(size = 8),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10, family = 'mono', face = 'bold', margin = margin(t = 3)),
        axis.ticks.x = element_blank(), axis.ticks.length.x = grid::unit(0, 'pt'),
        axis.title.y = element_text(size = 8.5),
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 1.0 * PT),
        axis.line = element_blank(), plot.margin = margin(6, 10, 8, 14))

# ---- B / C : Sp vs Wa response scatter (NS / both up / species-biased) ----
both_up_3h <- intersect(sp3$OrthoID[sp3$adj.P.Val < 0.05 & sp3$logFC > 0],
                        wa3$OrthoID[wa3$adj.P.Val < 0.05 & wa3$logFC > 0])
both_up_4d <- intersect(sp4$OrthoID[sp4$adj.P.Val < 0.05 & sp4$logFC > 0],
                        wa4$OrthoID[wa4$adj.P.Val < 0.05 & wa4$logFC > 0])
cat_palette <- c('NS' = 'gray75', 'Both up' = '#E69F00', 'Sp-biased' = '#8375a6', 'Wa-biased' = '#809c4c')

make_scatter <- function(sdf, both_set, time_label) {
  cats <- sdf %>%
    mutate(Category = case_when(
      OrthoID %in% both_set              ~ 'Both up',
      is.na(p_adj) | p_adj >= 0.05       ~ 'NS',
      Sp_diff_log > Wa_diff_log          ~ 'Sp-biased',
      TRUE                               ~ 'Wa-biased'),
      Category = factor(Category, levels = c('NS', 'Both up', 'Sp-biased', 'Wa-biased')))
  both_n <- sum(cats$Category == 'Both up'); sp_n <- sum(cats$Category == 'Sp-biased'); wa_n <- sum(cats$Category == 'Wa-biased')
  cat(sprintf('%s vs LL: both up %d, Sp-biased %d, Wa-biased %d\n', time_label, both_n, sp_n, wa_n))
  #   -> HL3h: 39 / 167 / 29 (main text L407-408); HL4d: 204 / 1792 / 435 (L397, L412, L420)
  ggplot(cats, aes(Sp_diff_log, Wa_diff_log, color = Category)) +
    geom_point(data = filter(cats, Category == 'NS'), size = 0.45, alpha = 0.35) +
    geom_point(data = filter(cats, Category != 'NS'), size = 0.85) +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = 'gray40', linewidth = 0.75 * PT) +
    geom_hline(yintercept = 0, linetype = 'dotted', color = 'gray60', linewidth = 0.5 * PT) +
    geom_vline(xintercept = 0, linetype = 'dotted', color = 'gray60', linewidth = 0.5 * PT) +
    scale_color_manual(values = cat_palette, breaks = c('Both up', 'Sp-biased', 'Wa-biased'),
                       labels = c(sprintf('Both up (%d)', both_n), sprintf('Sp-biased (%d)', sp_n),
                                  sprintf('Wa-biased (%d)', wa_n))) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
    labs(x = paste0('*S. polyrhiza* ', time_label, ' - LL<br>log2 fold change (TPM + 1)'),
         y = paste0('*W. australiana* ', time_label, ' - LL<br>log2 fold change (TPM + 1)'),
         color = NULL, title = paste0(time_label, ' vs LL response')) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = 'bold', size = 10),
          axis.title.x = ggtext::element_markdown(size = 8), axis.title.y = ggtext::element_markdown(size = 8),
          axis.text = element_text(size = 7.5),
          legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.text = element_text(size = 7),
          legend.background = element_rect(fill = scales::alpha('white', 0.7), color = NA),
          legend.key.size = grid::unit(3, 'mm'),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 1.5 * PT),
          axis.line = element_blank(), plot.margin = margin(6, 6, 6, 6))
}
B <- make_scatter(SLWL3h_summary,   both_up_3h, 'HL3h')
C <- make_scatter(SL4DWL3h_summary, both_up_4d, 'HL4d')

# ---- D-I : representative orthogroups, Sp | Wa facets ----
# Panels are addressed by Arabidopsis locus and resolved to the orthogroup through
# Orth_At (OrthoIDs are renumbered whenever OrthoFinder is re-run).
oid_of <- function(at) {
  o <- Orth_At$OrthoID[startsWith(Orth_At$At_Tran, paste0(at, '.'))][1]
  if (is.na(o) || !(o %in% logTPM$OrthoID)) stop('no expressed orthogroup for ', at)
  o
}
gene_pair <- function(at, gene_title) {
  oid <- oid_of(at)
  dat <- logTPM %>% filter(OrthoID == oid) %>%
    pivot_longer(cols = matches('^(Sp|Wa)_(SL_3h|SL_4D|WL_3h)_\\d+$'),
                 names_to = c('Species', 'Condition', 'Rep'),
                 names_pattern = '(Sp|Wa)_(SL_3h|SL_4D|WL_3h)_(\\d+)', values_to = 'Expr') %>%
    mutate(Species = factor(Species, levels = c('Sp', 'Wa')),
           Cond = factor(recode(Condition, WL_3h = 'LL', SL_3h = '3h', SL_4D = '4d'), levels = c('LL', '3h', '4d')))
  md   <- dat %>% group_by(Species, Cond) %>% summarise(M = mean(Expr), .groups = 'drop')
  padj <- SL4DWL3h_summary$p_adj[SL4DWL3h_summary$OrthoID == oid][1]
  plab <- if (is.na(padj)) 'NS' else paste0('*p* = ', formatC(padj, format = 'e', digits = 1))
  cat(sprintf('  %-7s %s -> %s  p_adj(HL4d vs LL) = %s\n', gene_title, at, oid, format(padj, digits = 2)))
  ggplot(dat, aes(Cond, Expr, color = Species)) +
    geom_line(data = md, aes(y = M, group = Species), linewidth = 1.5 * PT) +
    geom_point(data = md, aes(y = M), size = 1.8) +
    geom_jitter(width = 0.12, height = 0, size = 0.9, alpha = 0.5) +
    facet_wrap(~Species, ncol = 2, scales = 'free_y') +
    scale_color_manual(values = c(Sp = '#8375a6', Wa = '#809c4c')) +
    labs(title = paste0('**', gene_title, '**', ' (', plab, ')'), x = NULL, y = 'log2(TPM + 1)') +
    theme_classic(base_size = 9) +
    theme(legend.position = 'none', strip.background = element_blank(),
          strip.text = element_text(face = 'bold', size = 8),
          axis.text.x = element_text(size = 6.5), axis.text.y = element_text(size = 7),
          axis.title.y = element_text(size = 8),
          panel.border = element_rect(color = 'black', fill = NA, linewidth = 1.5 * PT),
          axis.line = element_blank(), panel.spacing = grid::unit(0.35, 'lines'),
          plot.title = ggtext::element_markdown(size = 9), plot.margin = margin(6, 5, 6, 5))
}
cat('panel orthogroups:\n')
D  <- gene_pair('AT3G22840', 'ELIP1')     # induced in both
E  <- gene_pair('AT3G15840', 'PIFI')      # Sp-biased
Fp <- gene_pair('AT5G42800', 'DFR')       # Sp-biased
G  <- gene_pair('AT4G32320', 'APX6')      # Sp-biased
H  <- gene_pair('AT1G22840', 'CYTC-1')    # Wa-biased
I  <- gene_pair('AT1G07770', 'RPS15A')    # Wa-biased

# ---- compose ----
scatter_row <- B + C + plot_layout(ncol = 2, widths = c(1, 1))
gene_row1   <- D + E + Fp + plot_layout(ncol = 3, widths = c(1, 1, 1))
gene_row2   <- G + H + I  + plot_layout(ncol = 3, widths = c(1, 1, 1))
final <- A / scatter_row / patchwork::free(gene_row1) / patchwork::free(gene_row2) +
  plot_layout(heights = c(0.96, 1.00, 0.72, 0.72)) +
  plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'))) &
  theme(plot.tag = element_text(face = 'bold', size = 11, hjust = 0, vjust = 1), plot.tag.position = c(0, 1))
ggsave(file.path(OUT, 'Fig5.pdf'), final, width = 6.5, height = 8, units = 'in')
cat('wrote output/Fig5.pdf\n')
