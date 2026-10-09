# Step 04: high-light response trajectory of every orthogroup (Fig. 5A, Data S8-1).
#
# Each of the two consecutive transitions (LL -> HL_3h, HL_3h -> HL_4d) is scored
# within species from the limma contrasts as U (up), F (flat) or D (down):
#   BH-adjusted p < 0.05 & |log2FC| >= 0.5  -> U or D, otherwise F
# giving nine possible trajectories (UU, UF, ..., DD) per species.
#
# Input : output/limma_<Sp|Wa>_SL_3h_vs_WL_3h.Rdata, output/limma_<Sp|Wa>_SL_4D_vs_SL_3h.Rdata (03)
# Output: output/transition_pattern_pergene.csv   per-orthogroup log2FC, p, pattern (both species)
#         output/transition_pattern_9x9.csv       Sp x Wa contingency table
# Run from this directory:  Rscript 04_transition_patterns.R
suppressPackageStartupMessages(library(tidyverse))
OUT <- 'output'

ld <- function(f) { e <- new.env(); load(f, envir = e); get('res', envir = e) }
Sp1 <- ld(file.path(OUT, 'limma_Sp_SL_3h_vs_WL_3h.Rdata'))
Sp2 <- ld(file.path(OUT, 'limma_Sp_SL_4D_vs_SL_3h.Rdata'))
Wa1 <- ld(file.path(OUT, 'limma_Wa_SL_3h_vs_WL_3h.Rdata'))
Wa2 <- ld(file.path(OUT, 'limma_Wa_SL_4D_vs_SL_3h.Rdata'))

PADJ <- 0.05; FC <- 0.5
cls <- function(lfc, padj) { padj[is.na(padj)] <- 1
  out <- rep('F', length(lfc)); out[padj < PADJ & lfc >=  FC] <- 'U'
  out[padj < PADJ & lfc <= -FC] <- 'D'; out }
mk <- function(d1, d2) {
  inner_join(d1 %>% transmute(OrthoID, l1 = logFC, p1 = adj.P.Val),
             d2 %>% transmute(OrthoID, l2 = logFC, p2 = adj.P.Val), by = 'OrthoID') %>%
    mutate(t1 = cls(l1, p1), t2 = cls(l2, p2), pat = paste0(t1, t2)) }
Sp <- mk(Sp1, Sp2); Wa <- mk(Wa1, Wa2)

J <- inner_join(Sp %>% transmute(OrthoID, Sp_l1 = l1, Sp_p1 = p1, Sp_l2 = l2, Sp_p2 = p2, Sp_pat = pat),
                Wa %>% transmute(OrthoID, Wa_l1 = l1, Wa_p1 = p1, Wa_l2 = l2, Wa_p2 = p2, Wa_pat = pat),
                by = 'OrthoID')
lev <- c('UU', 'UF', 'UD', 'FU', 'FF', 'FD', 'DU', 'DF', 'DD')
J <- J %>% mutate(Sp_pat = factor(Sp_pat, levels = lev),
                  Wa_pat = factor(Wa_pat, levels = lev),
                  shared = as.character(Sp_pat) == as.character(Wa_pat))
N <- nrow(J); ct <- table(Sp = J$Sp_pat, Wa = J$Wa_pat)
shr <- sum(diag(ct)); ff <- ct['FF', 'FF']

cat(sprintf('Rule: adj.P<%.2f & |logFC|>=%.1f -> U/D else F.  t1=LL->HL3h, t2=HL3h->HL4d\n', PADJ, FC))
cat(sprintf('Universe (expressed/testable both species, both transitions): N = %d\n\n', N))   # -> 10,916
cat(sprintf('Identical 9-pattern Sp==Wa : %d / %d = %.1f%%\n', shr, N, 100 * shr / N))
cat(sprintf('  FF/FF: %d (%.1f%%)  |  Identical EX FF/FF : %d / %d = %.1f%%\n\n',
            ff, 100 * ff / N, shr - ff, N - ff, 100 * (shr - ff) / (N - ff)))                 # -> ~3.8 % (main text L391)
cat('=== 9x9 contingency (row=Sp, col=Wa) ===\n'); print(addmargins(ct))
mg <- tibble(pattern = lev,
             Sp_n = as.integer(table(J$Sp_pat)[lev]),
             Wa_n = as.integer(table(J$Wa_pat)[lev]),
             shared_n = sapply(lev, function(p) sum(J$Sp_pat == p & J$Wa_pat == p))) %>%
      mutate(Sp_pct = round(100 * Sp_n / N, 1), Wa_pct = round(100 * Wa_n / N, 1),
             share_within_SpPat = ifelse(Sp_n > 0, round(100 * shared_n / Sp_n, 1), NA))
cat('\n=== Per-species marginal pattern usage ===\n'); print(as.data.frame(mg), row.names = FALSE)
#   -> main text L389-391: Sp UF 225; Wa FU 655, FD 238, DF 245; FF 95.8 % / 88.2 % (Fig. 5 legend)

write.csv(as.data.frame.matrix(addmargins(ct)), file.path(OUT, 'transition_pattern_9x9.csv'))
write.csv(J %>% mutate(across(c(Sp_pat, Wa_pat), as.character)),
          file.path(OUT, 'transition_pattern_pergene.csv'), row.names = FALSE)
cat('\nwrote output/transition_pattern_9x9.csv and output/transition_pattern_pergene.csv\n')
