library(rpkg)
library(here)
setwd('~/lab_repos/pi3k')


## load hotspots for pik3ca
pik3ca_hotspots_data <- fread(here('data/output.PIK3CA.txt'))
pik3ca_hotspots <- sortunique(pik3ca_hotspots_data$Amino_Acid_Position[pik3ca_hotspots_data$qvalue < 0.1])
pik3ca_hotspots <- as.integer(pik3ca_hotspots)
pik3ca_hotspots <- pik3ca_hotspots[!is.na(pik3ca_hotspots)]


## load hotspots for pik3cd
pik3cd_hotspots_data <- fread(here('data/output.PIK3CD.txt'))
pik3cd_hotspots <- sortunique(pik3cd_hotspots_data$Amino_Acid_Position[pik3cd_hotspots_data$qvalue < 0.1])
pik3cd_hotspots <- as.integer(pik3cd_hotspots)
pik3cd_hotspots <- pik3cd_hotspots[!is.na(pik3cd_hotspots)]


## load the alignment data from COBALT
al <- read.table('data/cobalt_alignment.fa',sep='\n')
pd <- paste(al$V1[2:20],collapse='')
pa <- paste(al$V1[22:40],collapse='')
al <- data.table(pik3cd=strsplit(pd,'')[[1]], pik3ca=strsplit(pa,'')[[1]])
al$pik3cd_pos <- NA
al$pik3cd_pos[al$pik3cd!='-'] <- 1:sum(al$pik3cd!='-')
al$pik3ca_pos <- NA
al$pik3ca_pos[al$pik3ca!='-'] <- 1:sum(al$pik3ca!='-')
al$pik3cd_tm <- paste0('PIK3CD ',al$pik3cd,al$pik3cd_pos)
al$pik3ca_tm <- paste0('PIK3CA ',al$pik3ca,al$pik3ca_pos)
al[grepl('-NA',pik3cd_tm),pik3cd_tm:='']
al[grepl('-NA',pik3ca_tm),pik3ca_tm:='']
al$i <- 1:nrow(al)

## annotate the alignment data with BLOSUM similarity
b <- fread(here('data/BLOSUM62.txt'))
b <- adt(melt(b))
names(b) <- c('aa1','aa2','score')

al <- merge(al, b, by.x=c('pik3cd','pik3ca'), by.y=c('aa1','aa2'), all.x=T)
setnames(al,'score','AD_score')
al <- merge(al, b[aa1==aa2,], by.x=c('pik3cd'), by.y=c('aa1'), all.x=T)
setnames(al,'score','DD_score')
al[,aa2:=NULL]
al <- al[order(i,decreasing=F),]

alleles <- c('R38C','R88C','G124D','N334K','A414V','C416R','E525K','R636Q','L806M','A835T','R894Q','G971E','E1021K','L1023R','E1025G')
d <- data.table(allele=alleles)
d$Reference_Amino_Acid <- substr(d$allele,1,1)
d$Amino_Acid_Position <- as.integer(substr(d$allele,2,(nchar(d$allele)-1)))
d$Variant_Amino_Acid <- substr(d$allele,nchar(d$allele),nchar(d$allele))
d$tm <- paste('PIK3CD',d$allele)
d <- d[order(Amino_Acid_Position),]
al <- merge(al, d[,c('Amino_Acid_Position','allele'),with=F], by.x='pik3cd_pos', by.y='Amino_Acid_Position', all.x=T)
al <- al[order(i),]
al$AD_score_imputed <- al$AD_score

## impute gaps with the gap opening score = -11
#al[is.na(AD_score_imputed), AD_score_imputed:=0]
k <- 10
i <- 1
n <- nrow(al)
al$sma <- as.numeric(NA)
for(i in 1:n) {
    vals <- (i-k):(i+k)
    vals <- vals[vals > 0 & vals < n]
    al$sma[i] <- mean(al$AD_score_imputed[vals],na.rm=T)
    }

library(ggrepel)
hs <- fread('data/aligned_residue_enrichment.txt')
al$aligned_hotspot <- (
                       al$pik3ca_tm %in% hs$pik3ca_residue[hs$aligned_hotspot==T] | 
                       al$pik3cd_tm %in% hs$pik3cd_residue[hs$aligned_hotspot==T] )

tmp <- al[!is.na(AD_score)] ## exclude gaps

library(ggplot2)
p <- ggplot(tmp, aes(x=aligned_hotspot, y=sma)) +
    geom_point(position=position_jitter(width=0.15,height=0)) +
    #geom_text_repel(aes(label=label)) +
    geom_boxplot(fill=NA)

x <- hs$AD_score[hs$aligned_hotspot==T]
y <- hs$AD_score[hs$aligned_hotspot==F]
wilcox.test(x,y)


al$notNA <- (!is.na(al$AD_score)) + 0
al$group <- cumsum(c(0,diff(al$notNA) != 0)) + 1
al <- al[!is.na(AD_score)]
al$group <- as.integer(factor(al$group, levels=unique(al$group)))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plots of alignment at important residues
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## specify the tested mutations to include in the plot (with a window of +/-15aa around each)
alleles <- c('R38C','R88C','G124D','N334K','A414V','C416R','E525K','R636Q','L806M','A835T','R894Q','G971E','E1021K','L1023R','E1025G')
d <- data.table(allele=alleles)
d$Reference_Amino_Acid <- substr(d$allele,1,1)
d$Amino_Acid_Position <- as.integer(substr(d$allele,2,(nchar(d$allele)-1)))
d$Variant_Amino_Acid <- substr(d$allele,nchar(d$allele),nchar(d$allele))
d$tm <- paste('PIK3CD',d$allele)
d <- d[order(Amino_Acid_Position),]
d$query_start <- d$Amino_Acid_Position - 15
d$query_end <- d$Amino_Acid_Position + 15


summarize_hotspot <- function(i, d, al) {
    allele <- d$allele[i]
    current_range <- seq(d$query_start[i],d$query_end[i])
    pik3ca_range <- al$pik3ca_pos[al$pik3cd_pos %in% current_range]
    tmp <- al[al$pik3cd_pos %in% current_range,]
    tmp <- tmp[order(pik3cd_pos),]
    tmp$pik3ca_tm <- paste0('PIK3CA ',tmp$pik3ca,tmp$pik3ca_pos)
    tmp$pik3cd_hs <- tmp$pik3cd_pos %in% pik3cd_hotspots
    tmp$pik3ca_hs <- tmp$pik3ca_pos %in% pik3ca_hotspots
    tmp$allele <- allele
    tmp$pik3ca_range <- paste0(min(pik3ca_range,na.rm=T),'-',max(pik3ca_range,na.rm=T))
    tmp$pik3cd_range <- paste0(min(current_range,na.rm=T),'-',max(current_range,na.rm=T))
    tmp$x <- 1:nrow(tmp)
    tmp
}

l <- lapply(1:nrow(d), summarize_hotspot, d, al)
ll <- rbindlist(l)
d1 <- melt(ll[,c('x','allele','pik3cd','pik3ca','pik3ca_range','pik3cd_range'),with=F], id.vars=c('x','allele','pik3ca_range','pik3cd_range'))
d2 <- melt(ll[,c('x','allele','DD_score','AD_score'),with=F], id.vars=c('x','allele'))
d3 <- melt(ll[,c('x','allele','pik3cd_hs','pik3ca_hs'),with=F], id.vars=c('x','allele'))

dd <- d1
setnames(dd,'value','AminoAcid')
setnames(dd,'variable','Gene')
dd$score <- d2$value
dd$hotspot <- d3$value
setnames(dd,'x','pos')
dd$Gene <- toupper(dd$Gene)
dd$Gene <- factor(dd$Gene, levels=rev(c('PIK3CD','PIK3CA')))
dd$allele <- paste0('PIK3CD ',dd$allele,' D:',dd$pik3cd_range,' A:',dd$pik3ca_range)
dd$allele <- factor(dd$allele, unique(dd$allele))
dd$`Hotspot residue` <- 'No'
dd[hotspot==T,`Hotspot residue`:='Yes']

cols <- c('#B2182B','white','#79a3cd')
hscols <- c('black','#B2182B'); names(hscols) <- c('No','Yes')

library(ggplot2)
p <- ggplot(dd, aes(x=Gene, y=pos)) +
    geom_tile(aes(fill=score)) +
    scale_fill_gradient2(low=cols[1],mid=cols[2],high=cols[3],midpoint=0,name='BLOSUM62\nsimilarity score') + 
    geom_text(aes(label=AminoAcid,color=`Hotspot residue`)) +
    scale_color_manual(values=hscols) + 
    facet_wrap(facets=~allele, ncol=2, nrow=8) +
    theme_std(base_size=14) +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.line.y=element_blank(), axis.ticks.y=element_blank()) + 
    labs(x=NULL,y=NULL,subtitle='PIK3CD hotspot paralogy') +
    coord_flip() 

ggsave('figures/pik3cd_hotspot_paralogy_figure.pdf',width=10,height=8)






