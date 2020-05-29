## graphically show % sequence identity between various related proteins
## using BLASTP with default parameters: BLOSUM62, gap costs: existence=11, extension=1



setwd('~/lab_repos/pi3k')

library(Biostrings)
library(data.table)
library(ggplot2)
library(rpkg)
library(cowplot)


## load/prep the protein sequences
pep <- fread('fa/proteins.fa',header=F,sep='\n')
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


## define helper function to align two sequences and get sequence matches
align <- function(pep,protein1,protein2) {
    aa1 <- pep$aa[pep$symbol==protein1]
    aa2 <- pep$aa[pep$symbol==protein2]

    ## peptide global alignment using BLOSUM62
    s1 <- AAString(aa1)
    s2 <- AAString(aa2)
    a <- pairwiseAlignment(s1, s2, substitutionMatrix = 'BLOSUM62', 
                           gapOpening = 11, gapExtension = 1)
    al1 <- as.character(pattern(a))
    al2 <- as.character(subject(a))
    d <- data.table(al1=strsplit(al1,'')[[1]], al2=strsplit(al2,'')[[1]])
    d <- cbind(i=1:nrow(d),d)
    d$pos1[d$al1!='-'] <- 1:sum(d$al1!='-')
    d$pos2[d$al2!='-'] <- 1:sum(d$al2!='-')
    d$identical <- as.integer(d$al1==d$al2)
    overall_sequence_identity <- sum(d$identical) / nrow(d)
    symmetric_moving_average <- function(i,d,n) {
        start <- i-n
        end <- i+n
        if(start < 1 | end > nrow(d)) {
            avg <- as.numeric(NA)
        } else {
            dd <- d[start:end]
            avg <- sum(dd$identical)/nrow(dd)
        }
        list(i=i,avg=avg)
    }
    l <- lapply(1:nrow(d),symmetric_moving_average,d,10)
    l <- rbindlist(l)
    d <- merge(d, l, by='i', all.x=T)
    d <- d[!is.na(avg),]
    d$percentage <- (1:nrow(d)) / nrow(d)
    d$protein1 <- protein1
    d$protein2 <- protein2
    list(data=d,overall_sequence_identity=overall_sequence_identity)
}

l <- align(pep,'PIK3CA','PIK3CD')
pi3k <- l$data
pi3k$alignment <- 'PIK3CA/D'
pi3k_id <- l$overall_sequence_identity

l <- align(pep,'MAP2K1','MAP2K2')
mek <- l$data
mek$alignment <- 'MEK1/2'
mek_id <- l$overall_sequence_identity

l <- align(pep,'KRAS','NRAS')
ras <- l$data
ras$alignment <- 'K/NRAS'
ras_id <- l$overall_sequence_identity

l <- align(pep,'ERBB2','EGFR')
her <- l$data
her$alignment <- 'ERBB2/EGFR'
her_id <- l$overall_sequence_identity




d <- rbind(pi3k, mek, ras, her)
d$avg <- d$avg * 100
d$labels <- paste0(d$protein1,' ',d$al1,d$pos1,'/',d$protein2,' ',d$al2,d$pos2)


## add track for overall sequence identity
dat <- data.table(alignment=c('PIK3CA/D','MEK1/2','K/NRAS','ERBB2/EGFR'),
                  pct=c(pi3k_id,mek_id,ras_id,her_id))
dat <- dat[order(pct,decreasing=F),]
dat$alignment <- factor(dat$alignment,levels=dat$alignment)
dat$pct <- 100*dat$pct
d$alignment <- factor(d$alignment, levels=dat$alignment)

## make plot showing the aligned protein on x-axis and moving average of % sequence identity on Y
plot_alignment <- function(query.alignment, d) { 
    tmp <- d[alignment==query.alignment,]
    mybreaks <- tmp$percentage[c(1,nrow(tmp))]
    mylabels <- tmp$labels[c(1,nrow(tmp))]
    tmp$pos <- factor(1:nrow(tmp))
    minx <- tmp$percentage[1]
    maxx <- tmp$percentage[nrow(tmp)]
    p1 <- ggplot(tmp, aes(x=percentage,y=avg)) + 
        scale_x_continuous(limits=c(minx,maxx),breaks=mybreaks,labels=mylabels) + 
        scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=c(0,100)) + 
        geom_area() + 
        theme_std(base_size=10) +
        theme(axis.line.x=element_blank(),axis.ticks.x=element_blank()) + 
        labs(x=NULL,y=query.alignment)

    tmp2 <- dat[dat$alignment==query.alignment,]
    p2 <- ggplot(tmp2,aes(x=alignment,y=pct)) +
        scale_y_continuous(expand=c(0,0),limits=c(0,100),breaks=seq(0,100,by=25)) + 
        geom_bar(stat='identity') +
        theme_std(base_size=12) +
        theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),
              axis.text.y=element_blank()) + 
        coord_flip() +
        labs(x=NULL,y='% sequence\nidentity')

    p <- plot_grid(p1,p2,align='h',ncol=2)
    p
}


alignments <- rev(levels(d$alignment))
l <- lapply(alignments, plot_alignment, d)
p <- plot_grid(plotlist=l,ncol=1)
p
ggsave('figures/alignment.pdf',width=6,height=8)


