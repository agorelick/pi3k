## graphically show % sequence identity between various related proteins
## using Needleman-Wunsch global alignment algorithm with the BLOSUM62 subsitution matrix, with gap costs: existence=11, extension=1


setwd('~/lab_repos/pi3k')

library(Biostrings)
library(data.table)
library(ggplot2)
library(cowplot)
library(here)

## load/prep the protein sequences
pep_RAS <- fread('original_data/cobalt/uniprot_RAS_sequences.fa',header=F,sep='\n')
pep_MEK <- fread('original_data/cobalt/uniprot_MEK_sequences.fa',header=F,sep='\n')
pep_EGF <- fread('original_data/cobalt/uniprot_EGF_sequences.fa',header=F,sep='\n')
pep_PI3K <- fread('original_data/cobalt/uniprot_PI3K_sequences.fa',header=F,sep='\n')

pep <- rbind(pep_RAS, pep_MEK, pep_EGF, pep_PI3K)
names(pep) <- 'seq'
pep$header <- grepl('>sp',pep$seq)
pep$protein <- cumsum(pep$header==T)
collapse <- function(pep) {
    header <- pep$seq[pep$header==T]
    txt <- strsplit(header,' ')[[1]]
    symbol <- gsub('GN=','',grep('GN=',txt,value=T))
    aa <- paste(pep$seq[pep$header==F],collapse='')
    list(symbol=symbol,aa=aa)
}
pep <- pep[,collapse(.SD),by=protein]
pep <- pep[!symbol %in% c('HRAS','PIK3CB')]
pep[symbol=='MAP2K1',symbol:='MEK1']
pep[symbol=='MAP2K2',symbol:='MEK2']
pep$symbol <- factor(pep$symbol, levels=unique(pep$symbol))
pep[symbol %in% c('KRAS','NRAS'),family:='RAS']
pep[symbol %in% c('MEK1','MEK2'),family:='MEK']
pep[symbol %in% c('ERBB2','EGFR'),family:='EGF']
pep[symbol %in% c('PIK3CA','PIK3CD'),family:='PI3K']


## for each pair of protein sequences in a family, globally align them
align <- function(pep) {
    proteins <- as.character(pep$symbol)
    aa1 <- pep$aa[pep$symbol==proteins[1]]
    aa2 <- pep$aa[pep$symbol==proteins[2]]

    ## run global peptide alignment using BLOSUM62
    s1 <- AAString(aa1)
    s2 <- AAString(aa2)
    a <- pairwiseAlignment(s1, s2, substitutionMatrix='BLOSUM80', gapOpening=10, gapExtension=1, type='global')
    al1 <- as.character(pattern(a))
    al2 <- as.character(subject(a))
    d <- data.table(al1=strsplit(al1,'')[[1]], al2=strsplit(al2,'')[[1]])
    d <- cbind(i=1:nrow(d),d)
    d$pos1[d$al1!='-'] <- 1:sum(d$al1!='-')
    d$pos2[d$al2!='-'] <- 1:sum(d$al2!='-')
    d$identical <- as.integer(d$al1==d$al2)
    symmetric_moving_average <- function(i,d,n) {
        start <- i-n
        end <- i+n
        region <- seq(start,end)
        region <- region[region >= 1 & region <= nrow(d)]
        dd <- d[region]
        n_identical <- sum(dd$identical)
        region_size <- nrow(dd)
        avg <- n_identical / region_size
        list(i=i,avg=avg,n_identical=n_identical,region_size=region_size)
    }
    l <- lapply(1:nrow(d),symmetric_moving_average,d,10)
    l <- rbindlist(l)
    d <- merge(d, l, by='i', all.x=T)
    d <- d[!is.na(avg),]
    d$percentage <- (1:nrow(d)) / nrow(d)
    d$protein1 <- proteins[1]
    d$protein2 <- proteins[2]
    d
}
d <- pep[,align(.SD),by=family]
d$avg <- d$avg * 100
d$labels <- paste0(d$protein1,' ',d$al1,d$pos1,'/',d$protein2,' ',d$al2,d$pos2)
d[,alignment:=paste0(protein1,'/',protein2)]
d$alignment <- factor(d$alignment, levels=unique(d$alignment))
d$family <- factor(d$family, levels=unique(d$family))
d$pik3cd_domain <- as.character(NA)

x <- d[protein2=='PIK3CD']
x[,pik3cd_domain:=NULL]
x$pik3cd_domain <- 'n/a'
x[is.na(pos2), pik3cd_domain:='Gap']
x[protein2=='PIK3CD' & pos2 %in% 32:107, pik3cd_domain:='p85 BD']
x[protein2=='PIK3CD' & pos2 %in% 175:281, pik3cd_domain:='RAS BD']
x[protein2=='PIK3CD' & pos2 %in% 337:464, pik3cd_domain:='PI3K C2']
x[protein2=='PIK3CD' & pos2 %in% 499:684, pik3cd_domain:='PI3Ka']
x[protein2=='PIK3CD' & pos2 %in% 773:991, pik3cd_domain:='PI3 PI4 kinase']
x$pik3cd_domain_int <- as.integer(factor(x$pik3cd_domain))
x$domain_change <- c(F,diff(x$pik3cd_domain_int) != 0)
x$domain_group <- cumsum(x$domain_change)
x[pik3cd_domain %in% c('n/a','Gap'), pik3cd_domain:=paste(pik3cd_domain, domain_group)]
d[protein2=='PIK3CD', pik3cd_domain:=x$pik3cd_domain]
d$pik3cd_domain <- factor(d$pik3cd_domain, levels=unique(d$pik3cd_domain))
ang::write_tsv(d,here('processed_data/global_pairwise_alignment_nw.txt'))


colors <- c('#2ecf00','#ff5454','#5c5cff','#ebd61c','#ba21e0')
names(colors) <- c('p85 BD','RAS BD','PI3K C2','PI3Ka','PI3 PI4 kinase')
remaining_levels <- levels(d$pik3cd_domain)[!levels(d$pik3cd_domain) %in% names(colors)]
na_toadd <- rep('#babdb5', sum(grepl('n/a',remaining_levels)))
names(na_toadd) <- grep('n/a',remaining_levels,value=T)
gap_toadd <- rep('black', sum(grepl('Gap',remaining_levels)))
names(gap_toadd) <- grep('Gap',remaining_levels,value=T)
colors <- c(colors, na_toadd, gap_toadd)

p <- ggplot(d, aes(x=percentage, y=avg)) +
    scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=c(0,100)) + 
    geom_area(data=d[protein2!='PIK3CD'], aes(x=percentage, y=avg), fill='black') +
    geom_segment(data=d[protein2=='PIK3CD'], aes(x=percentage, xend=percentage, y=0, yend=avg, color=pik3cd_domain), size=0.25) +
    scale_color_manual(values=colors,name='PIK3CD domain') +
    facet_wrap(facets=~alignment, ncol=1) +
    ang::theme_ang(base_size=12) +
    guides(color='none', fill='none') +
    theme(axis.line.x=element_blank(),axis.ticks.x=element_blank(), axis.text.x=element_blank()) + 
    labs(x='Global pairwise alignment',y='% sequence identity\n(Average in 21aa window)')

summarize_identity <- function(d) {
    n_identical <- sum(d$identical)
    n_aligned <- nrow(d)
    pct <- 100 * n_identical / n_aligned
    list(pct=pct)
}
overall <- d[,summarize_identity(.SD),by=alignment]

p_side <- ggplot(overall, aes(x=alignment, y=pct)) +
    scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0,100,by=25)) + 
    geom_bar(stat='identity',fill='black') +
    facet_wrap(facets=~alignment, ncol=1,scale='free_y') +
    ang::theme_ang(base_size=12) +
    theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(), axis.text.y=element_blank()) + 
    coord_flip() +
    labs(y='% sequence identity',x=NULL)

p_combo <- plot_grid(p, p_side, align='h', nrow=1, rel_widths=c(4,1))
ggsave('figures/global_pairwise_alignment.pdf',width=8,height=5)



