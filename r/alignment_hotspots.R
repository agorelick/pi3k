setwd('~/lab_repos/pi3k')
library(here)
library(ggplot2)
library(rpkg)



## load the alignment data from COBALT
al <- read.table('data/cobalt_alignment.fa',sep='\n')
pd <- paste(al$V1[2:20],collapse='')
pa <- paste(al$V1[22:40],collapse='')
al <- data.table(pik3cd=strsplit(pd,'')[[1]], pik3ca=strsplit(pa,'')[[1]])
al <- cbind(aligned_pos=1:nrow(al), al)
al$pik3cd_pos <- NA
al$pik3cd_pos[al$pik3cd!='-'] <- 1:sum(al$pik3cd!='-')
al$pik3ca_pos <- NA
al$pik3ca_pos[al$pik3ca!='-'] <- 1:sum(al$pik3ca!='-')
al$pik3cd_tm <- paste0('PIK3CD ',al$pik3cd_pos)
al$pik3ca_tm <- paste0('PIK3CA ',al$pik3ca_pos)
al[grepl('NA',pik3cd_tm),pik3cd_tm:='']
al[grepl('NA',pik3ca_tm),pik3ca_tm:='']
al$i <- 1:nrow(al)


## append with hotspots data for aligned positions
## exclude positions with no mutations in either pik3ca/pik3cd
hotspots_a <- fread(here('data/output.PIK3CA.txt'),select=c('Hugo_Symbol','qvalue','Mutation_Count','Total_Samples','probability','Analysis_Type','tm','Amino_Acid_Position'))
hotspots_d <- fread(here('data/output.PIK3CD.txt'),select=c('Hugo_Symbol','qvalue','Mutation_Count','Total_Samples','probability','Analysis_Type','tm','Amino_Acid_Position'))
d <- merge(al, hotspots_d, by.x='pik3cd_tm', by.y='tm',all.x=T)
d <- merge(d, hotspots_a, by.x='pik3ca_tm', by.y='tm',all.x=T,allow.cartesian=T)
d <- d[pik3ca_tm!='' & pik3cd_tm!='',] ## subset for residues present in both proteins (no gaps)
names(d) <- gsub('[.]x','_pik3cd',names(d))
names(d) <- gsub('[.]y','_pik3ca',names(d))
d <- d[!is.na(Mutation_Count_pik3ca) & !is.na(Mutation_Count_pik3cd)]## exclude positions with no mutations in either protein
d$pik3cd_tm <- paste0('PIK3CD ',d$pik3cd,d$pik3cd_pos)
d$pik3ca_tm <- paste0('PIK3CA ',d$pik3ca,d$pik3ca_pos)
pik3ca_hotspots <- unique(d$pik3ca_tm[d$qvalue_pik3ca < 0.1])
pik3cd_hotspots <- unique(d$pik3cd_tm[d$qvalue_pik3cd < 0.1])
d$pik3ca_hotspot <- d$pik3ca_tm %in% pik3ca_hotspots
d$pik3cd_hotspot <- d$pik3cd_tm %in% pik3cd_hotspots
d <- d[Analysis_Type_pik3cd=='pancan' & Analysis_Type_pik3cd==Analysis_Type_pik3ca,] ## only compare A/D within the same analysis type
d <- d[order(aligned_pos),]


## annotate each aligned position/test with the binomial parameter for A OR B
#d$parameter <- d$probability_pik3ca + d$probability_pik3cd - d$probability_pik3ca*d$probability_pik3cd
d$parameter <- d$probability_pik3ca + d$probability_pik3cd #
d$aligned_mutations <- d$Mutation_Count_pik3ca + d$Mutation_Count_pik3cd
d$aligned_samples <- d$Total_Samples_pik3cd ## this is the same as Total_Samples_pik3ca
d$i <- 1:nrow(d)


## test every position within the same analysis-type
test <- function(d) {
    tst <- binom.test(d$aligned_mutations,d$aligned_samples,d$parameter,alternative='greater')
    d$aligned_pval <- tst$p.value
    d
}
dd <- d[,test(.SD),by=i]


## get q-value within analysis-type
correct_within_type <- function(d) {
    d$aligned_qval <- p.adjust(d$aligned_pval,method='BH')
    d
}
dd <- dd[,correct_within_type(.SD),by=Analysis_Type_pik3cd]
dd <- dd[order(aligned_pval,decreasing=F),]
dd$aligned_hotspot <- dd$aligned_qval < 0.1


## format the output table
dd <- dd[,c('aligned_pos','pik3ca_tm','pik3cd_tm','pik3ca_hotspot','pik3cd_hotspot','aligned_hotspot',
           'aligned_qval','Mutation_Count_pik3cd','Mutation_Count_pik3ca','aligned_mutations',
           'aligned_samples','parameter','Analysis_Type_pik3ca'),with=F]
names(dd) <- c('aligned_pos','pik3ca_residue','pik3cd_residue','pik3ca_ever_hotspot','pik3cd_ever_hotspot',
               'aligned_hotspot','aligned_qvalue',
               'pik3cd_mutated','pik3ca_mutated','aligned_samples_mutated',
               'total_samples','probability','aligned_hotspot_analysis_type')
dd <- dd[order(aligned_pos),]
write.tsv(dd,here('data/aligned_residue_enrichment.txt'))


## plot results
#experimental_alleles <- paste('PIK3CD',c('R38','R88','G124','N334','A414','C416','E525','R636',
#                                          'L806','A835','R894','G971','E1021','L1023','E1025'))

dd <- fread(here('data/aligned_residue_enrichment.txt'))
dd$nlog10q <- -log10(dd$aligned_qvalue)
dd <- dd[order(nlog10q,aligned_samples_mutated,decreasing=T),]

positions_to_plot <- dd$aligned_pos[dd$pik3ca_ever_hotspot==T | 
                                    dd$pik3cd_ever_hotspot==T |
                                    dd$aligned_hotspot==T] 
                                    #dd$pik3cd_residue %in% experimental_alleles]
pd <- dd[aligned_pos %in% positions_to_plot,]

pd$aligned_pos <- factor(pd$aligned_pos,levels=rev(pd$aligned_pos))
pd$pik3ca_residue <- gsub('PIK3CA ','',pd$pik3ca_residue)
pd$pik3cd_residue <- gsub('PIK3CD ','',pd$pik3cd_residue)

pd$pik3cd_hotspot <- 'Not significant'
pd[pik3cd_ever_hotspot==T,pik3cd_hotspot:='Hotspot']
pd$pik3ca_hotspot <- 'Not significant'
pd[pik3ca_ever_hotspot==T,pik3ca_hotspot:='Hotspot']
pd$aligned_hotspot <- as.character(pd$aligned_hotspot)
pd[aligned_hotspot=='TRUE',aligned_hotspot:='Hotspot']
pd[aligned_hotspot=='FALSE',aligned_hotspot:='Not significant']


cols <- c('#BFBFBF','#ff4040')
names(cols) <- c('Not significant','Hotspot')


p_signif <- ggplot(pd,aes(x=aligned_pos,y=nlog10q)) +
    geom_bar(stat='identity',fill='black',color='white') +
    theme_std(base_size=16) +
    scale_y_continuous(expand=c(0,0),position='right') + 
    coord_flip() +
    labs(x=NULL,y='-log10(Q-value)') +
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())



#### N mutations in PIK3CA
yinfo <- break_axis(pd$pik3ca_mutated,lowerticksize=100,upperticksize=500,maxlower=500,
                    #minupper=500,
                    ratio_lower_to_upper=0.666)
pd$y <- yinfo$newy
p_n_pik3ca <- ggplot(pd,aes(x=aligned_pos,y=y)) +
    geom_bar(stat='identity',fill='black',color='white') +
    theme_std(base_size=16) +
    scale_y_continuous(expand=c(0,0),position='right',limits=yinfo$limits,breaks=yinfo$breaks,
                       labels=yinfo$labels) + 
    coord_flip() +
    labs(x=NULL,y='PIK3CA\nmutated') +
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
          axis.text.x=element_text(angle=-90,hjust=1,vjust=-0.5))



#### N mutations in PIK3CD
#yinfo <- break_axis(pd$pik3cd_mutated,lowerticksize=50,upperticksize=500,maxlower=250,
#                    ratio_lower_to_upper=0.666)
#pd$y <- yinfo$newy
p_n_pik3cd <- ggplot(pd,aes(x=aligned_pos,y=pik3cd_mutated)) +
    geom_bar(stat='identity',fill='black',color='white') +
    theme_std(base_size=16) +
    scale_y_continuous(expand=c(0,0),position='right') +
    coord_flip() +
    labs(x=NULL,y='PIK3CD\nmutated') +
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
          axis.text.x=element_text(angle=-90,hjust=1,vjust=-0.5))

p_pik3cd <- ggplot(pd,aes(x=aligned_pos,y=1)) +
    geom_tile(aes(fill=pik3cd_hotspot),color='white') + 
    scale_fill_manual(values=cols,name=NULL) + 
    geom_text(aes(label=pik3cd_residue)) +
    theme_std(base_size=16) +
    coord_flip() +
    theme(
          axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),axis.line.y=element_blank(),
          legend.position='none') +
    labs(x=NULL,y='PIK3CD') +
    scale_y_continuous(expand=c(0,0),position='right') 


p_pik3ca <- ggplot(pd,aes(x=aligned_pos,y=1)) +
    geom_tile(aes(fill=pik3ca_hotspot),color='white') + 
    scale_fill_manual(values=cols,name=NULL) + 
    geom_text(aes(label=pik3ca_residue)) +
    theme_std(base_size=16) +
    coord_flip() +
    theme(
          axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),axis.line.y=element_blank(),
          legend.position='none') +
    labs(x=NULL,y='PIK3CA') +
    scale_y_continuous(expand=c(0,0),position='right') 


p_aligned <- ggplot(pd,aes(x=aligned_pos,y=1)) +
    geom_tile(aes(fill=aligned_hotspot),color='white') + 
    scale_fill_manual(values=cols,name=NULL) + 
    #geom_text(aes(label=aligned_residue)) +
    theme_std(base_size=16) +
    coord_flip() +
    theme(
          axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),axis.line.y=element_blank(),
          legend.position='none') +
    labs(x=NULL,y='Aligned residue') +
    scale_y_continuous(expand=c(0,0),position='right') 


library(cowplot)
p <- plot_grid(p_pik3cd,p_n_pik3cd,p_pik3ca,p_n_pik3ca,p_aligned,p_signif,
               ncol=6,rel_widths=c(1,1,1,1,1,3),align='h')
p
ggsave(here('figures/aligned_residue_enrichment_conservative_v2.pdf'),width=8,height=10)











## If I want to include data for non-pancan results...

## append with hotspots data for aligned positions
hotspots_a <- fread(here('data/output.PIK3CA.txt'),select=c('Hugo_Symbol','qvalue','Mutation_Count','Total_Samples','probability','Analysis_Type','tm','Amino_Acid_Position'))
hotspots_d <- fread(here('data/output.PIK3CD.txt'),select=c('Hugo_Symbol','qvalue','Mutation_Count','Total_Samples','probability','Analysis_Type','tm','Amino_Acid_Position'))
d <- merge(al, hotspots_d, by.x='pik3cd_tm', by.y='tm',all.x=T)
d <- merge(d, hotspots_a, by.x='pik3ca_tm', by.y='tm',all.x=T,allow.cartesian=T)
d <- d[pik3ca_tm!='' & pik3cd_tm!='',] ## subset for residues present in both proteins (no gaps)
names(d) <- gsub('[.]x','_pik3cd',names(d))
names(d) <- gsub('[.]y','_pik3ca',names(d))
d <- d[!is.na(Mutation_Count_pik3ca) & !is.na(Mutation_Count_pik3cd)]## exclude positions with no mutations in either protein
d$pik3cd_tm <- paste0('PIK3CD ',d$pik3cd,d$pik3cd_pos)
d$pik3ca_tm <- paste0('PIK3CA ',d$pik3ca,d$pik3ca_pos)
pik3ca_hotspots <- unique(d$pik3ca_tm[d$qvalue_pik3ca < 0.1])
pik3cd_hotspots <- unique(d$pik3cd_tm[d$qvalue_pik3cd < 0.1])
d$pik3ca_hotspot <- d$pik3ca_tm %in% pik3ca_hotspots
d$pik3cd_hotspot <- d$pik3cd_tm %in% pik3cd_hotspots


## get best q-value for pik3ca, pik3cd, and aligned residues
get_best_q <- function(dd) {
    best_pik3ca_q <- min(dd$qvalue_pik3ca,na.rm=T)
    best_pik3ca_types <- paste(sortunique(dd$Analysis_Type_pik3ca[dd$qvalue_pik3ca==best_pik3ca_q]),collapse=',')
    best_pik3cd_q <- min(dd$qvalue_pik3cd,na.rm=T)
    best_pik3cd_types <- paste(sortunique(dd$Analysis_Type_pik3cd[dd$qvalue_pik3cd==best_pik3cd_q]),collapse=',')
    dd$best_pik3ca_q <- best_pik3ca_q
    dd$best_pik3ca_types <- best_pik3ca_types
    dd$best_pik3cd_q <- best_pik3cd_q
    dd$best_pik3cd_types <- best_pik3cd_types
    dd
}
d <- d[,get_best_q(.SD),by=i]
d[best_pik3ca_q >= 0.1,best_pik3ca_types:='Not significant']
d[best_pik3cd_q >= 0.1,best_pik3cd_types:='Not significant']
d <- d[Analysis_Type_pik3cd==Analysis_Type_pik3ca,] ## only compare A/D within the same analysis type


## annotate each aligned position/test with the binomial parameter for A OR B
d$parameter <- d$probability_pik3ca + d$probability_pik3cd - d$probability_pik3ca*d$probability_pik3cd
d$aligned_mutations <- d$Mutation_Count_pik3ca + d$Mutation_Count_pik3cd
d$aligned_samples <- d$Total_Samples_pik3cd ## this is the same as Total_Samples_pik3ca
d$i <- 1:nrow(d)


## test every position within the same analysis-type
test <- function(d) {
    tst <- binom.test(d$aligned_mutations,d$aligned_samples,d$parameter,alternative='greater')
    d$aligned_pval <- tst$p.value
    d
}
dd <- d[,test(.SD),by=i]


## get q-value within analysis-type
correct_within_type <- function(d) {
    d$aligned_qval <- p.adjust(d$aligned_pval,method='BH')
    d
}
dd <- dd[,correct_within_type(.SD),by=Analysis_Type_pik3cd]
dd <- dd[order(aligned_pval,decreasing=F),]
dd$aligned_hotspot <- dd$aligned_qval < 0.1


get_best_aligned_q <- function(dd) {
    best_aligned_q <- min(dd$aligned_qval,na.rm=T)
    best_aligned_types <- paste(sortunique(dd$Analysis_Type_pik3ca[dd$aligned_qval==best_aligned_q]),collapse=',')
    dd$best_aligned_q <- best_aligned_q
    dd$best_aligned_types <- best_aligned_types
    dd
}
dd <- dd[,get_best_aligned_q(.SD),by=i]
dd[best_aligned_q >= 0.1,best_aligned_types:='Not significant']
dd$aligned_hotspot <- dd$best_aligned_q < 0.1


## format the output table
dd <- dd[,c('i','pik3ca_tm','pik3cd_tm','pik3ca_hotspot','pik3cd_hotspot','aligned_hotspot',
            'best_pik3cd_q','best_pik3ca_q','best_aligned_q',
            'best_pik3cd_types','best_pik3ca_types','best_aligned_types',
           'Mutation_Count_pik3cd','Mutation_Count_pik3ca','aligned_mutations',
           'aligned_samples','parameter'),with=F]


names(dd) <- c('aligned_pos','pik3ca_residue','pik3cd_residue','pik3ca_hotspot','pik3cd_hotspot',
               'aligned_hotspot','aligned_qvalue','analysis_type','pik3ca_mutated','pik3cd_mutated','aligned_mutated',
               'total_samples','binomial_parameter')
dd <- dd[order(aligned_pos),]


## annotate which residues are EVER hotspots
pik3ca_hotspots <- unique(hotspots_a$tm[hotspots_a$qvalue < 0.1])
pik3cd_hotspots <- unique(dd$pik3cd_residue[dd$pik3cd_qvalue < 0.1])





pplot(p_pik3cd)


pplot(p)





pplot(p)




d <- data.table(allele=alleles)
d$Reference_Amino_Acid <- substr(d$allele,1,1)
d$Amino_Acid_Position <- as.integer(substr(d$allele,2,(nchar(d$allele)-1)))
d$Variant_Amino_Acid <- substr(d$allele,nchar(d$allele),nchar(d$allele))
d$tm <- paste('PIK3CD',d$allele)
d <- d[order(Amino_Acid_Position),]
d$query_start <- d$Amino_Acid_Position - 15
d$query_end <- d$Amino_Acid_Position + 15








d[is.na(Mutation_Count_pik3cd),Mutation_Count_pik3cd:=0]
d[is.na(Mutation_Count_pik3ca),Mutation_Count_pik3ca:=0]









hotspots <- rbind(hotspots_a,hotspots_d)
hotspots <- hotspots[,c('Hugo_Symbol','qvalue','Mutation_Count','Total_Samples','probability','Analysis_Type','tm','Amino_Acid_Position'),with=F]


(hotspots[,c('tm','Analysis_Type','Mutation_Count'),with=F],idvar='tm',timevar='Analysis_Type',direction='wide')
summarize_position_analysis <- function(hs) {
    

}





