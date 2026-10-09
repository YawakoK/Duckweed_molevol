# Fig. 1C: relative molecular evolutionary rate per species.
#
# Input : data/fig1_relative_rates.csv  (06_summarise_absrel.py)
#         one row per single-copy orthogroup (2,819), each column a species'
#         root-to-tip branch length divided by that of Colocasia esculenta in
#         the same gene tree.
# Output: output/Fig1C_relative_rate.pdf/.png, output/Fig1C_medians_letters.csv
#
# Statistics (Methods): Friedman test (species as treatments, gene trees as
# blocks) followed by post-hoc pairwise Wilcoxon signed-rank tests with Holm
# correction.  Compact letters are assigned from the pairwise matrix and
# ordered by decreasing median rate; rows follow the tip order of Fig. 1B.
# Run from this directory:  Rscript 08_fig1c_relative_rate.R
suppressPackageStartupMessages({ library(tidyverse) })

IN  <- Sys.getenv('IN',  'data')     # set IN=output to use freshly regenerated tables
OUT <- 'output'
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

wide <- read.csv(file.path(IN, 'fig1_relative_rates.csv')) %>% drop_na()
sp   <- setdiff(colnames(wide), 'OG')
cat('gene trees:', nrow(wide), '  species:', paste(sp, collapse = ' '), '\n')

# ---- Friedman test (blocks = gene trees) ----
m  <- as.matrix(wide[, sp])
ft <- friedman.test(m)
cat(sprintf('Friedman: chi2 = %.1f, df = %d, p = %.3g\n',
            ft$statistic, ft$parameter, ft$p.value))

# ---- post-hoc pairwise Wilcoxon signed-rank, Holm ----
pairs <- combn(sp, 2, simplify = FALSE)
pv <- sapply(pairs, function(p)
  wilcox.test(m[, p[1]], m[, p[2]], paired = TRUE)$p.value)
pv <- p.adjust(pv, method = 'holm')
names(pv) <- sapply(pairs, function(p) paste(p, collapse = '-'))
cat('pairwise comparisons:', length(pv),
    '  all p < 0.05:', all(pv < 0.05), '\n')          # -> main text L170 "all pairs p < 0.05"

# compact letters, assigned greedily: species sharing a letter are NOT
# significantly different from each other under the Holm-adjusted tests.
diff_pair <- function(a, b) {
  k1 <- paste(a, b, sep = '-'); k2 <- paste(b, a, sep = '-')
  pp <- if (k1 %in% names(pv)) pv[[k1]] else pv[[k2]]
  pp < 0.05
}
ord    <- names(sort(sapply(sp, function(x) median(m[, x])), decreasing = TRUE))
groups <- list()
for (s_ in ord) {
  placed <- FALSE
  for (g in seq_along(groups)) {
    if (all(!sapply(groups[[g]], function(o) diff_pair(s_, o)))) {
      groups[[g]] <- c(groups[[g]], s_); placed <- TRUE; break
    }
  }
  if (!placed) groups[[length(groups) + 1]] <- s_
}
cld <- setNames(rep('', length(sp)), sp)
for (g in seq_along(groups))
  for (s_ in groups[[g]]) cld[s_] <- paste0(cld[s_], letters[g])

med <- wide %>% pivot_longer(all_of(sp), names_to = 'Species', values_to = 'ratio') %>%
  group_by(Species) %>% summarise(med = median(ratio), .groups = 'drop') %>%
  arrange(desc(med))
print(as.data.frame(med), row.names = FALSE)   # -> main text L172-173: Sp 0.250, Wo 1.260

# Rows follow the tip order of the tree in Fig. 1B so the two panels line up;
# only the compact letters are ordered by decreasing median rate.
tree_order <- c('Sp', 'La', 'Lgib', 'Laeq', 'Wa', 'Wo')   # top -> bottom in Fig. 1B
stopifnot(setequal(tree_order, sp))
lvl <- rev(tree_order)                                    # ggplot y goes bottom-up
long <- wide %>%
  pivot_longer(all_of(sp), names_to = 'Species', values_to = 'ratio') %>%
  mutate(Species = factor(Species, levels = lvl))

lab <- tibble(Species = factor(names(cld), levels = lvl),
              letter  = unname(cld))

p <- ggplot(long, aes(y = Species, x = ratio)) +
  geom_violin(fill = '#2C7FB8', color = 'gray30', linewidth = 0.25,
              trim = TRUE, scale = 'width', adjust = 1.5) +
  geom_boxplot(width = 0.12, fill = 'white', outlier.shape = NA,
               linewidth = 0.2) +
  geom_text(data = lab, aes(x = 2.95, y = Species, label = letter),
            inherit.aes = FALSE, hjust = 1, size = 6 / ggplot2::.pt) +
  coord_cartesian(xlim = c(0, 3)) +
  labs(x = 'Relative Branch length', y = NULL) +
  scale_x_continuous(breaks = c(0, 1, 2, 3)) +
  theme_bw(base_size = 6) +
  theme(panel.grid.minor = element_blank(),
        panel.border  = element_rect(linewidth = 0.4),
        panel.grid.major.x = element_line(linewidth = 0.2),
        panel.grid.major.y = element_blank(),
        axis.text.y   = element_blank(),   # species are identified by the tree in (B)
        axis.ticks.y  = element_blank(),
        axis.text.x   = element_text(size = 6, color = 'black'),
        axis.title.x  = element_text(size = 6.5, face = 'bold'),
        axis.ticks.x  = element_line(linewidth = 0.25),
        plot.margin   = margin(1, 1, 1, 1))

# Final panel size for Fig. 1C as laid out with the tree in Fig. 1B.
ggsave(file.path(OUT, 'Fig1C_relative_rate.pdf'), p, width = 3.22, height = 4.39, units = 'cm')
ggsave(file.path(OUT, 'Fig1C_relative_rate.png'), p, width = 3.22, height = 4.39, units = 'cm', dpi = 600)

write.csv(tibble(Species = names(cld), letter = unname(cld)) %>%
            left_join(med, by = 'Species') %>% arrange(desc(med)),
          file.path(OUT, 'Fig1C_medians_letters.csv'), row.names = FALSE)
cat('wrote', file.path(OUT, 'Fig1C_relative_rate.pdf/.png'), '\n')
