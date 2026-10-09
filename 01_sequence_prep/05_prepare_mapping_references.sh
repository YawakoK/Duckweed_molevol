#!/bin/bash
# Combined nuclear + organellar references for S. polyrhiza (Sp) and W. australiana (Wa):
# used for read mapping (07_rnaseq/01) and as the Sp / Wa protein and CDS sets for
# OrthoFinder and the codon analyses (02_phylogenomics).
#
# Sources
#   Sp nuclear    Sp_polyrhiza_9509-REF-OXFORD-3.0 genome + CSHL2022v1 gene models (Lemna Genome Hub)
#   Sp plastid    NC_015891 ; Sp mitochondrion NC_017840            (NCBI GenBank + GFF)
#   Wa nuclear    Wo_australiana_8730-REF-CSHL-1.0 genome + CSHL2022v1 gene models (Lemna Genome Hub)
#   Wa plastid    NC_015899                                          (NCBI)
#   Wa mitochondrion: the chrM contig of the Wa assembly, gene models called with MFannot
#                 (https://github.com/BFL-lab/Mfannot, default parameters) and converted with
#                 agat_convert_mfannot2gff.pl
#
# Steps (AGAT 1.0.0; gffread; seqkit)
#   1. standardise each organellar GFF with agat_convert_sp_gxf2gxf.pl and prefix feature IDs
#      (Sp_chrC_/Sp_chrM_/Wa_chrC_/Wa_chrM_) so that they cannot collide with nuclear IDs
#   2. extract organellar protein (-p) and CDS with agat_sp_extract_sequences.pl
#   3. concatenate genome FASTAs and GFF3s; remove tRNA/rRNA features from the combined GFF3
#      (-> *_withorga_SP.gff3, the annotation given to featureCounts)
#   4. extract CDS from the combined annotation (agat_sp_extract_sequences.pl -t cds) and
#      protein sets = distributed primary proteins + organellar proteins
#   5. gene positions for 07_rnaseq/data: awk '$3=="gene"{print $1,$4,$5,$9}' | sed 's/ID=//'
#
# The exact commands as run (2025-05) are reproduced below for Sp; Wa is identical
# with the Wa files, plus the MFannot conversion.
set -euo pipefail

# ---- 1-2. organellar gene models ----
for g in Sp_Chroloplast_NC_015891 Sp_Mitochondria_NC_017840 Wa_Chroloplast_NC_015899; do
  agat_convert_sp_gxf2gxf.pl -g $g.gff -o $g.agated.gff
done
sed -i "s/=nbis/=Sp_chrM_nbis/g" Sp_Mitochondria_NC_017840.agated.gff
sed -i "s/=nbis/=Sp_chrC_nbis/g" Sp_Chroloplast_NC_015891.agated.gff
sed -i "s/=nbis/=Wa_chrC_nbis/g" Wa_Chroloplast_NC_015899.agated.gff
for g in Sp_Chroloplast_NC_015891 Sp_Mitochondria_NC_017840 Wa_Chroloplast_NC_015899; do
  agat_sp_extract_sequences.pl -g $g.agated.gff -f $g.fasta -o $g.agated.protein.faa -p
  agat_sp_extract_sequences.pl -g $g.agated.gff -f $g.fasta -o $g.agated.cds.fa
  awk '{print $1}' $g.agated.protein.faa > $g.agated.protein.renamed.faa
done

# Wa mitochondrion: MFannot -> GFF3 -> protein / CDS (tRNA and rRNA dropped)
seqkit grep -n -p "chrM" Wo_australiana_8730-REF-CSHL-1.0.fasta > Wo_chrM.fasta
agat_convert_mfannot2gff.pl -m mfannot_Wo_chrM.fasta.new -o Wa_mfannot.gff3
grep -v tRNA Wa_mfannot.gff3 | grep -v rRNA > Wa_mfannot.mRNA.gff3
agat_convert_sp_gxf2gxf.pl -g Wa_mfannot.mRNA.gff3 -o Wa_Mitochondria_mfannot.agated.gff
agat_sp_extract_sequences.pl -g Wa_Mitochondria_mfannot.agated.gff -f Wo_chrM.fasta -t exon -p \
    -o Wa_Mitochondria_mfannot.protein.faa
agat_sp_extract_sequences.pl -g Wa_Mitochondria_mfannot.agated.gff -f Wo_chrM.fasta -t exon \
    -o Wa_Mitochondria_mfannot.cds_.fa
sed -E 's/^>.*transcript=([^ ]+).*/>\1/' Wa_Mitochondria_mfannot.protein.faa > Wa_Mitochondria_mfannot.protein.renamed.fasta
sed -E 's/^>.*transcript=([^ ]+).*/>\1/' Wa_Mitochondria_mfannot.cds_.fa    > Wa_Mitochondria_mfannot.cds.fa

# ---- 3. combined references ----
# Sp: chromosome names in the gff3 (chr1..chr20) are harmonised to the genome FASTA (Chr01..Chr20)
sed -E 's/\bchr10\b/Chr10/g;s/\bchr11\b/Chr11/g;s/\bchr12\b/Chr12/g;s/\bchr13\b/Chr13/g;s/\bchr14\b/Chr14/g;s/\bchr15\b/Chr15/g;s/\bchr16\b/Chr16/g;s/\bchr17\b/Chr17/g;s/\bchr18\b/Chr18/g;s/\bchr19\b/Chr19/g;s/\bchr20\b/Chr20/g;s/\bchr1\b/Chr01/g;s/\bchr2\b/Chr02/g;s/\bchr3\b/Chr03/g;s/\bchr4\b/Chr04/g;s/\bchr5\b/Chr05/g;s/\bchr6\b/Chr06/g;s/\bchr7\b/Chr07/g;s/\bchr8\b/Chr08/g;s/\bchr9\b/Chr09/g;' \
    Sp_polyrhiza_9509-REF-OXFORD-3.0_CSHL2022v1.genes.gff3 > Sp_polyrhiza_9509-REF-OXFORD-3.0_CSHL2022v1.genes.renamed.gff3
cat Sp9509_oxford_v3.fasta Sp_Mitochondria_NC_017840.fasta Sp_Chroloplast_NC_015891.fasta > Sp9509_oxford_v3_withorga.fasta
cat Sp_polyrhiza_9509-REF-OXFORD-3.0_CSHL2022v1.genes.renamed.gff3 Sp_Mitochondria_NC_017840.agated.gff \
    Sp_Chroloplast_NC_015891.agated.gff > Sp_polyrhiza_9509-REF-OXFORD-3.0_CSHL2022v1.withorga.gff3
grep -v trna Sp_polyrhiza_9509-REF-OXFORD-3.0_CSHL2022v1.withorga.gff3 | grep -v rrna \
    > Sp_polyrhiza_9509-REF-OXFORD-3.0_CSHL2022v1.withorga_SP.gff3

# Wa: the assembly already carries chrC and chrM contigs; they are replaced by the
# NCBI plastid record and the MFannot-annotated chrM
seqkit fx2tab -n Wo_australiana_8730-REF-CSHL-1.0.fasta | grep -v chrC | grep -v chrM > noCnoM.tsv
seqkit grep -f noCnoM.tsv Wo_australiana_8730-REF-CSHL-1.0.fasta > Wo_australiana_8730-REF-CSHL-1.0_noCnoM.fasta
cat Wo_australiana_8730-REF-CSHL-1.0_noCnoM.fasta Wo_chrM.fasta Wa_Chroloplast_NC_015899.fasta \
    > Wo_australiana_8730-REF-CSHL-1.0_withorga.fasta
cat Wo_australiana_8730-REF-CSHL-1.0_CSHL2022v1.genes.gff3 Wa_Mitochondria_mfannot.agated.gff \
    Wa_Chroloplast_NC_015899.agated.gff > Wo_australiana_8730-REF-CSHL-1.0_CSHL2022v1.withorga.gff3
grep -v trna Wo_australiana_8730-REF-CSHL-1.0_CSHL2022v1.withorga.gff3 | grep -v rrna \
    > Wo_australiana_8730-REF-CSHL-1.0_CSHL2022v1.withorga_SP.gff3

# ---- 4. protein and CDS sets for OrthoFinder / codon analyses ----
cat Sp_polyrhiza_9509-REF-OXFORD-3.0_CSHL2022v1.genes.proteins.primary.fasta \
    Sp_Mitochondria_NC_017840.agated.protein.renamed.faa \
    Sp_Chroloplast_NC_015891.agated.protein.renamed.faa > Sp_withorga.fasta
cat Wo_australiana_8730-REF-CSHL-1.0_CSHL2022v1.genes.proteins.primary.fasta \
    Wa_Mitochondria_mfannot.protein.renamed.fasta \
    Wa_Chroloplast_NC_015899.agated.protein.renamed.faa > Wa_withorga.fasta

agat_sp_extract_sequences.pl -g Sp_polyrhiza_9509-REF-OXFORD-3.0_CSHL2022v1.withorga_SP.gff3 \
    -f Sp9509_oxford_v3.fasta -t cds -o Sp9509_oxford_v3.agted.cds.fa
agat_sp_extract_sequences.pl -g Wo_australiana_8730-REF-CSHL-1.0_CSHL2022v1.withorga_SP.gff3 \
    -f Wo_australiana_8730-REF-CSHL-1.0_noCnoM.fasta -t cds -o Wo_australiana_8730-REF-CSHL-1.0_noCnoM.agted.cds.fa
cat Sp9509_oxford_v3.agted.cds.fa Sp_Mitochondria_NC_017840.agated.cds.fa Sp_Chroloplast_NC_015891.agated.cds.fa \
    | awk '{print $1}' > Sp_withorga.cds.renamed.fa
cat Wo_australiana_8730-REF-CSHL-1.0_noCnoM.agted.cds.fa Wa_Mitochondria_mfannot.cds.fa Wa_Chroloplast_NC_015899.agated.cds.fa \
    | awk '{print $1}' > Wa_withorga.cds.renamed.fa

# ---- 5. gene positions (07_rnaseq/data/{Sp,Wa}_gene_positions.txt) ----
awk '$3=="gene"{print $1,$4,$5,$9}' Sp_polyrhiza_9509-REF-OXFORD-3.0_CSHL2022v1.withorga_SP.gff3 | sed 's/ID=//g' > Sp_gene_positions.txt
awk '$3=="gene"{print $1,$4,$5,$9}' Wo_australiana_8730-REF-CSHL-1.0_CSHL2022v1.withorga_SP.gff3 | sed 's/ID=//g' > Wa_gene_positions.txt
