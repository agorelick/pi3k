rm(list=ls())

setwd('~/lab_repos/pi3k')
library(here)
library(ggplot2)
library(rpkg)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get a summary table of the aligned residues based on 3 methods
# Needleman-Wunsch, BLASTp, COBALT (aligning PIK3CA, PIK3CD, PIK3CB)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## load NW pairwise alignment for PI3K
d <- fread('processed_data/global_pairwise_alignment_nw.txt')
d <- d[family=='PI3K']
d[,pik3ca_tm:=paste0(al1,pos1)]
d[,pik3cd_tm:=paste0(al2,pos2)]
d[is.na(pos1),pik3ca_tm:='-']
d[is.na(pos2),pik3cd_tm:='-']

res <- d[,c('pos2','pik3cd_tm','pos1','pik3ca_tm'),with=F]
names(res) <- c('PIK3CD_pos','PIK3CD_aa','PIK3CA_pos_NW','PIK3CA_NW')
res <- res[PIK3CD_aa != '-']

## load the alignment data from BLASTp
al <- read.table('original_data/cobalt/blastp_PIK3CD_PIK3CA_BLOSUM80.fa',sep='\n')
pd <- paste(al$V1[2:19],collapse='')
pa <- paste(al$V1[21:38],collapse='')
al <- data.table(pik3cd=strsplit(pd,'')[[1]], pik3ca=strsplit(pa,'')[[1]])
al <- cbind(aligned_pos=1:nrow(al), al)
al$pik3cd_pos <- NA
al$pik3cd_pos[al$pik3cd!='-'] <- 1:sum(al$pik3cd!='-') + 19
al$pik3ca_pos <- NA
al$pik3ca_pos[al$pik3ca!='-'] <- 1:sum(al$pik3ca!='-') + 19
al$pik3cd_tm <- paste0(al$pik3cd,al$pik3cd_pos)
al$pik3ca_tm <- paste0(al$pik3ca,al$pik3ca_pos)
al[grepl('NA',pik3cd_tm),pik3cd_tm:='-']
al[grepl('NA',pik3ca_tm),pik3ca_tm:='-']
al$i <- 1:nrow(al)

result <- merge(res, al[,c('pik3cd_tm','pik3ca_tm'),with=F], by.x='PIK3CD_aa', by.y='pik3cd_tm', all.x=T)
result <- result[order(PIK3CD_pos),]
setnames(result,'pik3ca_tm','PIK3CA_BLAST')

## add the paralogous PIK3CA and PIK3CB residues based on the COBALT alignments
cob <- fread('processed_data/cobalt_alignment_PI3K_annotated.tsv')
cob[,PIK3CA_tm:=paste0(PIK3CA,PIK3CA_pos)]
cob[,PIK3CB_tm:=paste0(PIK3CB,PIK3CB_pos)]
cob[,PIK3CD_tm:=paste0(PIK3CD,PIK3CD_pos)]
cob[grepl('NA',PIK3CA_tm),PIK3CA_tm:='-']
cob[grepl('NA',PIK3CB_tm),PIK3CB_tm:='-']
cob[grepl('NA',PIK3CD_tm),PIK3CD_tm:='-']
cob <- cob[,c('PIK3CD_tm','PIK3CA_tm','PIK3CB_tm'),with=F]
result <- merge(result, cob, by.x='PIK3CD_aa', by.y='PIK3CD_tm', all.x=T)
result <- result[order(PIK3CD_pos),]
setnames(result,'PIK3CA_tm','PIK3CA_COBALT')
setnames(result,'PIK3CB_tm','PIK3CB_COBALT')
result <- result[,c(2,1,4:7,3),with=F]
result$aligned <- F
result[PIK3CA_NW==PIK3CA_BLAST,aligned:=T]


# ~~~~~~~~
# run PROS
# ~~~~~~~~

## append with hotspots data for aligned positions
## exclude positions with no mutations in either pik3ca/pik3cd
hotspots_a <- fread(here('original_data/mutational_data/output.PIK3CA.txt'),select=c('Hugo_Symbol','qvalue','Mutation_Count','Total_Samples','probability','Analysis_Type','tm','Amino_Acid_Position'))
hotspots_d <- fread(here('original_data/mutational_data/output.PIK3CD.txt'),select=c('Hugo_Symbol','qvalue','Mutation_Count','Total_Samples','probability','Analysis_Type','tm','Amino_Acid_Position'))

collapsed_residue <- function(hs) {
    sample_count <- hs$Total_Samples[hs$Analysis_Type=='pancan']
    mutation_count <- hs$Mutation_Count[hs$Analysis_Type=='pancan']
    binomial_parameter <- hs$probability[hs$Analysis_Type=='pancan']
    pancan <- any(hs$qvalue < 0.1 & hs$Analysis_Type=='pancan')
    cancertype <- any(hs$qvalue < 0.1 & hs$Analysis_Type!='pancan')
    q_pancan <- hs$qvalue[hs$Analysis_Type=='pancan']
    q_cancertype <- min(hs$qvalue[hs$Analysis_Type!='pancan'])
    if(pancan==T) {
        highest <- 'Hotspot (pan-cancer)'
    } else if(cancertype) {
        highest <- 'Hotspot (cancer-specific)'
    } else {
        highest <- 'Mutated'
    }
    list(class=highest, mutation_count=mutation_count, sample_count=sample_count, binomial_parameter=binomial_parameter, q_pancan=q_pancan, q_cancertype=q_cancertype)
}
hotspots_a <- hotspots_a[,collapsed_residue(.SD),by=Amino_Acid_Position]
hotspots_a$Amino_Acid_Position <- as.integer(hotspots_a$Amino_Acid_Position)
hotspots_a <- hotspots_a[!is.na(Amino_Acid_Position)]
hotspots_d <- hotspots_d[,collapsed_residue(.SD),by=Amino_Acid_Position]
hotspots_d$Amino_Acid_Position <- as.integer(hotspots_d$Amino_Acid_Position)
hotspots_d <- hotspots_d[!is.na(Amino_Acid_Position)]


## merge in hotspots data for PIK3CA and PIK3CD
d <- merge(result, hotspots_d, by.x='PIK3CD_pos', by.y='Amino_Acid_Position',all.x=T)
d <- merge(d, hotspots_a, by.x='PIK3CA_pos_NW', by.y='Amino_Acid_Position',all.x=T)
names(d) <- gsub('[.]x','_PIK3CD',names(d))
names(d) <- gsub('[.]y','_PIK3CA',names(d))
d[is.na(mutation_count_PIK3CA),mutation_count_PIK3CA:=0]
d[is.na(mutation_count_PIK3CD),mutation_count_PIK3CD:=0]


## annotate each aligned position/test with the binomial parameter for A OR B
d$binomial_parameter <- d$binomial_parameter_PIK3CA + d$binomial_parameter_PIK3CD
d$aligned_mutations <- d$mutation_count_PIK3CA + d$mutation_count_PIK3CD
d$aligned_samples <- d$sample_count_PIK3CD ## this is the same as for PIK3CA
d$i <- 1:nrow(d)
d$PROS_test <- F
d[aligned==T & mutation_count_PIK3CA > 0 & mutation_count_PIK3CD > 0, PROS_test:=T]
d[PROS_test==F,binomial_parameter:=NA]
d[PROS_test==F,aligned_mutations:=NA]
d[PROS_test==F,aligned_samples:=NA]
setnames(d,'class_PIK3CA','PIK3CA_NW_hotspot')
setnames(d,'class_PIK3CD','PIK3CD_hotspot')

## test every paralogous position
test <- function(d) {
    p <- tryCatch({
        binom.test(d$aligned_mutations,d$aligned_samples,d$binomial_parameter,alternative='greater')$p.value
    }, error=function(e) {
        as.numeric(NA)
    })
    d$aligned_pval <- p
    d
}
dd <- d[,test(.SD),by=i]

## get q-value
dd$aligned_qval <- p.adjust(dd$aligned_pval,method='BH')
dd$aligned_hotspot <- dd$aligned_qval < 0.1

## format the output table
setnames(dd,'mutation_count_PIK3CD','PIK3CD_mutated')
setnames(dd,'mutation_count_PIK3CA','PIK3CA_NW_mutated')
setnames(dd,'q_pancan_PIK3CA','PIK3CA_NW_q_pancan')
setnames(dd,'q_cancertype_PIK3CA','PIK3CA_NW_q_cancertype')
setnames(dd,'q_pancan_PIK3CD','PIK3CD_q_pancan')
setnames(dd,'q_cancertype_PIK3CD','PIK3CD_q_cancertype')

out <- dd[,c('PIK3CD_pos','PIK3CD_aa','PIK3CA_NW','PIK3CA_BLAST','PIK3CA_COBALT','PIK3CB_COBALT',
             'aligned','aligned_hotspot','aligned_pval','aligned_qval','binomial_parameter',
             'aligned_mutations','aligned_samples','PROS_test',
             'PIK3CD_mutated',
             'PIK3CD_hotspot',
             'PIK3CD_q_pancan',
             'PIK3CD_q_cancertype',
             'PIK3CA_NW_mutated',
             'PIK3CA_NW_hotspot',
             'PIK3CA_NW_q_pancan',
             'PIK3CA_NW_q_cancertype'),with=F]
out[PROS_test==T,]
out <- out[order(PIK3CD_pos),]
write.tsv(out,here('processed_data/aligned_residue_summary.txt'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

negative_alleles <- c('L806','S367','R894','G971','R338','D911')

#dd <- fread(here('processed_data/aligned_residue_enrichment_blosum80_blastp_intersect.txt'))
dd <- fread(here('processed_data/aligned_residue_summary.txt'))
dd <- dd[PROS_test==T,]
dd$nlog10q <- -log10(dd$aligned_qval)

positions_to_plot <- dd[grepl('Hotspot',PIK3CD_hotspot) | grepl('Hotspot',PIK3CA_NW_hotspot) | aligned_hotspot==T | PIK3CD_aa %in% negative_alleles,(PIK3CD_pos)]
nrow(dd) - length(positions_to_plot)

pd <- dd[PIK3CD_pos %in% positions_to_plot,]
pd <- pd[order(nlog10q,aligned_mutations,decreasing=F),]
pd$aligned_pos <- 1:nrow(pd)

cols <- c('#BFBFBF','#6599ff','#0433ff')
names(cols) <- c('Mutated','Hotspot (cancer-specific)','Hotspot (pan-cancer)')

p_signif <- ggplot(pd,aes(x=aligned_pos,y=nlog10q)) +
    geom_bar(stat='identity',fill='black') + 
    ang::theme_ang(base_size=12) +
    scale_y_continuous(expand=c(0,0),position='right') + 
    coord_flip() +
    labs(x=NULL,y='-log10(Q-value)') +
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())

#### N mutations in PIK3CA
yinfo <- break_axis(pd$PIK3CA_NW_mutated,lowerticksize=100,upperticksize=500,maxlower=500, ratio_lower_to_upper=0.666)
pd$y <- yinfo$newy
p_n_pik3ca <- ggplot(pd,aes(x=aligned_pos,y=y)) +
    geom_bar(stat='identity',fill='black') + 
    ang::theme_ang(base_size=12) +
    scale_y_continuous(expand=c(0,0),position='right',limits=yinfo$limits,breaks=yinfo$breaks,
                       labels=yinfo$labels) + 
    coord_flip() +
    labs(x=NULL,y='PIK3CA\nmutated') +
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
          axis.text.x=element_text(angle=-90,hjust=1,vjust=-0.5))

#### N mutations in PIK3CD
p_n_pik3cd <- ggplot(pd,aes(x=aligned_pos,y=PIK3CD_mutated)) +
    geom_bar(stat='identity',fill='black') + 
    ang::theme_ang(base_size=12) +
    scale_y_continuous(expand=c(0,0),position='right') +
    coord_flip() +
    labs(x=NULL,y='PIK3CD\nmutated') +
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
          axis.text.x=element_text(angle=-90, hjust=1))

## hotspot annotations in PIK3CD
p_pik3cd <- ggplot(pd,aes(x=aligned_pos,y=1)) +
    geom_tile(aes(fill=PIK3CD_hotspot),color='white') + 
    scale_fill_manual(values=cols,name=NULL) + 
    geom_text(aes(label=PIK3CD_aa)) +
    ang::theme_ang(base_size=12) +
    coord_flip() +
    theme(
          axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),axis.line.y=element_blank(),
          legend.position='none') +
    labs(x=NULL,y='PIK3CD') +
    scale_y_continuous(expand=c(0,0),position='right') 

## hotspot annotations in PIK3CA
p_pik3ca <- ggplot(pd,aes(x=aligned_pos,y=1)) +
    geom_tile(aes(fill=PIK3CA_NW_hotspot),color='white') + 
    scale_fill_manual(values=cols,name=NULL) + 
    geom_text(aes(label=PIK3CA_NW)) +
    ang::theme_ang(base_size=12) +
    coord_flip() +
    theme(
          axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),axis.line.y=element_blank(),
          legend.position='none') +
    labs(x=NULL,y='PIK3CA') +
    scale_y_continuous(expand=c(0,0),position='right') 

library(cowplot)
p <- plot_grid(p_pik3cd,p_n_pik3cd,p_pik3ca,p_n_pik3ca,p_signif,
               ncol=6,rel_widths=c(1,1,1,1,3),align='h')
ggsave(here('figures/aligned_residue_enrichment_conservative_v4.pdf'),width=8,height=10)




