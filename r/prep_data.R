setwd('~/lab_repos/pi3k')
library(data.table)
library(ggplot2)
library(cowplot)

maf <- fread('data/input.pik3cabd.maf.txt.gz')
maf <- maf[Hugo_Symbol %in% c('PIK3CA','PIK3CD'),]
maf <- maf[FILTER %in% c('PASS','.'),]

required.fields <- c('Tumor_Sample_Barcode','Hugo_Symbol','Chromosome','Start_Position','End_Position',
'Reference_Allele','Tumor_Seq_Allele2','Variant_Type','Variant_Classification','tm','Amino_Acid_Change',
'oncotree_organtype','Center')

maf <- maf[,(required.fields),with=F]
maf <- maf[order(Tumor_Sample_Barcode,Hugo_Symbol,Start_Position),]
write.tsv(maf,'data/data_mutations.txt')
