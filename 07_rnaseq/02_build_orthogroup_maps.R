# Step 02: orthogroup <-> gene maps from the OrthoFinder table.
#
# Input : data/Orthogroups.tsv          (02_phylogenomics/01, eight species)
#         data/Sp_gene_positions.txt, data/Wa_gene_positions.txt (chr start end gene)
# Output: output/Orth_ath.Rdata  full table with a stable OrthoID (Orth00001.. in file order)
#         output/Orth_At.Rdata   OrthoID x Arabidopsis transcript
#         output/Orth_Sp.Rdata   OrthoID x S. polyrhiza gene (+ position, copy number)
#         output/Orth_Wa.Rdata   OrthoID x W. australiana gene (+ position, copy number)
# Run from this directory:  Rscript 02_build_orthogroup_maps.R
suppressPackageStartupMessages({ library(tidyverse) })
OUT <- 'output'; dir.create(OUT, showWarnings = FALSE)

Orth_ath <- read.csv('data/Orthogroups.tsv', sep = '\t', header = TRUE)
cat('Orthogroups.tsv: rows=', nrow(Orth_ath), ' cols=', ncol(Orth_ath), '\n', sep = '')
stopifnot('Araport11_pep_20250214' %in% colnames(Orth_ath))
stopifnot('Sp_withorga' %in% colnames(Orth_ath))
stopifnot('Wa_withorga' %in% colnames(Orth_ath))

# OrthoID: numbered in the row order of Orthogroups.tsv (OrthoFinder renumbers OGs
# every run, so OrthoIDs are only meaningful together with this table).
Orth_ath$OrthoID <- paste0('Orth', formatC(1:length(Orth_ath$Orthogroup), width = 5, flag = '0'))

Orth_At <- data.frame(OrthoID = Orth_ath$OrthoID,
                      At_Tran = Orth_ath$Araport11_pep_20250214) %>%
  separate_rows(At_Tran, sep = ', ') %>%
  filter(At_Tran != '' & At_Tran != '*')
tmp <- as.data.frame(table(Orth_At$OrthoID)); colnames(tmp) <- c('OrthoID', 'At_CN')
Orth_At <- merge(Orth_At, tmp)

Orth_Sp <- data.frame(OrthoID = Orth_ath$OrthoID, Sp_Tran = Orth_ath$Sp_withorga) %>%
  separate_rows(Sp_Tran, sep = ', ') %>%
  filter(Sp_Tran != '' & Sp_Tran != '*')
tmp <- as.data.frame(table(Orth_Sp$OrthoID)); colnames(tmp) <- c('OrthoID', 'Sp_CN')
Orth_Sp <- merge(Orth_Sp, tmp)
Orth_Sp <- cbind(Sp_Gene = sapply(seq_len(nrow(Orth_Sp)),
                                  function(x) strsplit(Orth_Sp$Sp_Tran[x], '_T0')[[1]][1]),
                 Orth_Sp)

Orth_Wa <- data.frame(OrthoID = Orth_ath$OrthoID, Wa_Tran = Orth_ath$Wa_withorga) %>%
  separate_rows(Wa_Tran, sep = ', ') %>%
  filter(Wa_Tran != '' & Wa_Tran != '*')
tmp <- as.data.frame(table(Orth_Wa$OrthoID)); colnames(tmp) <- c('OrthoID', 'Wa_CN')
Orth_Wa <- merge(Orth_Wa, tmp)
Orth_Wa <- cbind(Wa_Gene = sapply(seq_len(nrow(Orth_Wa)),
                                  function(x) strsplit(Orth_Wa$Wa_Tran[x], '_T0')[[1]][1]),
                 Orth_Wa)

# gene positions (chr, start, end, gene) from the combined GFF3s
Sp_pos <- read.csv('data/Sp_gene_positions.txt', sep = ' ', header = FALSE)
colnames(Sp_pos) <- c('Sp_chr', 'Sp_start', 'Sp_end', 'Sp_Gene')
Sp_pos$Sp_Gene <- sapply(seq_len(nrow(Sp_pos)),
                         function(x) strsplit(Sp_pos$Sp_Gene[x], '_T0')[[1]][1])
Orth_Sp <- merge(Sp_pos, Orth_Sp)

Wa_pos <- read.csv('data/Wa_gene_positions.txt', sep = ' ', header = FALSE)
colnames(Wa_pos) <- c('Wa_chr', 'Wa_start', 'Wa_end', 'Wa_Gene')
Wa_pos$Wa_Gene <- sapply(seq_len(nrow(Wa_pos)),
                         function(x) strsplit(Wa_pos$Wa_Gene[x], '_T0')[[1]][1])
Orth_Wa <- merge(Wa_pos, Orth_Wa)

cat('Orth_ath rows:', nrow(Orth_ath),
    ' | Orth_At rows:', nrow(Orth_At),
    ' | Orth_Sp rows:', nrow(Orth_Sp),
    ' | Orth_Wa rows:', nrow(Orth_Wa), '\n')

save(Orth_ath, file = file.path(OUT, 'Orth_ath.Rdata'))
save(Orth_At,  file = file.path(OUT, 'Orth_At.Rdata'))
save(Orth_Sp,  file = file.path(OUT, 'Orth_Sp.Rdata'))
save(Orth_Wa,  file = file.path(OUT, 'Orth_Wa.Rdata'))
cat('wrote output/Orth_{ath,At,Sp,Wa}.Rdata\n')
