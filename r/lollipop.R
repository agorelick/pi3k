setwd('~/lab_repos/pi3k')

library(data.table)

d <- fread('data/input.pik3cabd.maf.txt.gz')
d <- d[Hugo_Symbol=='PIK3CD']
setnames(d,'Tumor_Seq_Allele2','Variant_Allele')
d <- d[Variant_Classification=='Missense_Mutation' & Variant_Type=='SNP',]

required_fields <- c('Chromosome','Start_Position','End_Position','Reference_Allele','Variant_Allele')
dd <- d[,(required_fields),with=F]
rpkg::write.tsv(dd,'data/pik3cd_lollipop_data.txt')


