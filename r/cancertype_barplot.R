## Barplot showing prevalence of PIK3CD mutations across cancer types from IMPACT data
## also shade based on known hotspot mutations
##
## here we are using IMPACT samples from the 52K MAF (impact only because they are annotated for MSI)
##
## 2020-05-25


library(rpkg); library(data.table)
library(ggplot2); library(ggrepel);
library(magrittr)

setwd('~/lab_repos/pi3k')


## we want to restrict ourselves to MSKIMPACT samples
## which were among the 52K data freeze for hotspots, so load the sample IDs from the 52K MAF
## then load the most up-to-date MSKIMPACT data and subset for these PATIENTS

supermaf_samples <- fread('data/data_samples.txt')
valid_samples <- grep('P-00',supermaf_samples$Tumor_Sample_Barcode,value=T)

## subset clinical data for valid patients
info <- fread('data/data_clinical_sample_20200119.txt',skip=4,
              select=c('SAMPLE_ID','PATIENT_ID','CANCER_TYPE','CANCER_TYPE_DETAILED','MSI_SCORE','GENE_PANEL'))
info <- info[PATIENT_ID %in% valid_samples,]

## exclude MSI-High samples
info$MSI_SCORE <- as.numeric(info$MSI_SCORE)
info <- info[info$MSI_SCORE < 10 | is.na(info$MSI_SCORE),]
info[is.na(CANCER_TYPE),CANCER_TYPE:='']

## collapse rarer tumortypes into Other
tbl <- table.freq(info$CANCER_TYPE[info$CANCER_TYPE %nin% c('','Cancer of Unknown Primary')])
top30 <- tbl$value[1:30]
info$CANCER_TYPE_LABEL <- info$CANCER_TYPE
info[CANCER_TYPE %nin% top30,CANCER_TYPE_LABEL:='Other']
tbl <- table.freq(info$CANCER_TYPE_LABEL)

## load mutation data, subset for the retained samples
d <- fread('data/data_mutation_extended_20200119.txt',
           select=c('Tumor_Sample_Barcode','Hugo_Symbol','Variant_Classification','Variant_Type',
                    'HGVSp_Short'))
d <- d[Tumor_Sample_Barcode %in% info$SAMPLE_ID,]
d <- d[Hugo_Symbol=='PIK3CD',]
#tbl <- table.freq(d$Variant_Classification)
d <- d[Hugo_Symbol=='PIK3CD' & Variant_Classification=='Missense_Mutation' & Variant_Type=='SNP']
d$Amino_Acid_Position <- as.integer(substr(d$HGVSp_Short,4,(nchar(d$HGVSp_Short)-1)))
d <- merge(d, info, by.x='Tumor_Sample_Barcode', by.y='SAMPLE_ID', all.x=T)

## annotate PIK3CD hotspots
hs <- fread('data/output.PIK3CD.txt')
hotspots <- sort(unique(as.integer(hs$Amino_Acid_Position[hs$qvalue < 0.1])))
d$hotspot <- d$Amino_Acid_Position %in% hotspots

summarize_cancertype <- function(d) {
    samples_mutated <- length(unique(d$Tumor_Sample_Barcode))
    samples_hotspot <- length(unique(d$Tumor_Sample_Barcode[d$hotspot==T]))
    list(samples_mutated=samples_mutated,samples_hotspot=samples_hotspot)
}
muts <- d[,summarize_cancertype(.SD),by=CANCER_TYPE_LABEL]

dat <- merge(tbl, muts, by.x='value', by.y='CANCER_TYPE_LABEL', all.x=T)
dat[is.na(dat)] <- 0
names(dat)[1:2] <- c('type','samples_overall')
write.tsv(dat,'data/cancertype_barplot_data.txt')


dat$samples_mutated_no_hotspot <- dat$samples_mutated - dat$samples_hotspot
dat$prop_samples_hotspot <- dat$samples_hotspot / dat$samples_overall
dat$prop_samples_mutated_no_hotspot <- dat$samples_mutated_no_hotspot / dat$samples_overall
dat$prop_mutated <- dat$samples_mutated/dat$samples_overall
dat <- dat[order(dat$prop_mutated,decreasing=T),]

m <- melt(dat[,c(1,6,7),with=F])
m$sample <- 'Non-hotspot'
m$sample[m$variable=='prop_samples_hotspot'] <- 'Hotspot'
m$sample <- factor(m$sample, levels=rev(c('Non-hotspot','Hotspot')))
orders <- c(dat$type[dat$type!='Other'],'Other')
m$type <- factor(m$type,levels=orders)
mycols <- c('grey80','black'); names(mycols) <- c('Non-hotspot','Hotspot')

p <- ggplot(m,aes(x=type,y=value)) +
    geom_bar(stat='identity',aes(fill=sample)) +
    theme_std(base_size=16) + 
    scale_y_continuous(expand=c(0,0),labels=scales::percent) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
          legend.position='right') +
    labs(x=NULL,y='Samples',subtitle='52K MAF subset for MSK-IMPACT samples') +
    scale_fill_manual(values=mycols)
p

ggsave('figures/barplot_cancertypes_props.pdf',width=10,height=7)





