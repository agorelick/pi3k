library(data.table)
library(ggplot2)
library(RColorBrewer)

rm(list=ls())
setwd('~/lab_repos/pi3k')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# annotate COBALT alignment
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

annotate_cobalt_alignment <- function(alignment_fa, aa_table) {
    ## extract the aligned amino acid sequences and format them in a long table
    fa <- fread(alignment_fa, header=F)
    fa$header <- F
    fa[grep('^>',V1),header:=T]
    fa$group <- cumsum(fa$header)
    collapse <- function(fa) {
        label <- fa$V1[fa$header==T]
        sequence <- paste0(fa$V1[fa$header==F],collapse='') 
        sequence <- strsplit(sequence,'')[[1]]
        out <- data.table(label=label,sequence=sequence)
        out$global_position <- 1:nrow(out)
        ind <- which(out$sequence!='-')
        out$local_position[ind] <- 1:length(ind)
        out
    }
    fa <- fa[,collapse(.SD),by=c('group')]
    fa[,group:=NULL]
    fa$label <- gsub('>','',fa$label)
    fa$label <- factor(fa$label, levels=unique(fa$label))
    fa2 <- dcast(global_position ~ label, value.var='sequence', data=fa)
    isoforms <- names(fa2)[2:ncol(fa2)]

    ## load the amino acid table
    aa <- fread(aa_table)[,(2:3),with=F]
    aa <- rbind(aa, data.table(amino_acid='-',side_chain='Gap'))

    ## annotate the local position for each isoform in the global alignment
    for(i in isoforms) {
        ind <- which(fa2[[i]]!='-')
        fa2[[paste0(i,'_pos')]][ind] <- 1:length(ind)
    }

    ## annotate the amino acid side-chain type for each isoform
    for(i in isoforms) {
        fa2 <- merge(fa2, aa, by.x=i, by.y='amino_acid', all.x=T)
        setnames(fa2,'side_chain',paste0(i,'_side_chain'))
    }
    fa2 <- fa2[order(global_position),]

    ## reorder fields in output
    fields <- c('global_position', isoforms, paste0(isoforms,'_pos'), paste0(isoforms,'_side_chain'))
    fa2 <- fa2[,(fields),with=F]
    fa2
}


pairwise_comparison_of_sidechain_types <- function(fa2, isoforms) {

    ## get all unique pairwise combinations of the submitted isoforms
    size <- length(isoforms)
    m <- matrix(nrow=size, ncol=size)
    rownames(m) <- isoforms; colnames(m) <- isoforms
    m[upper.tri(m)] <- 1
    for(i in 1:nrow(m)) m[i,i] <- NA
    m <- as.data.table(reshape2::melt(m))
    m <- m[!is.na(value),]
    f=function(i, m, fa2) {
        from=as.character(m$Var1[i])
        to=as.character(m$Var2[i])
        from_sc <- paste0(from,'_side_chain')
        to_sc <- paste0(to,'_side_chain')
        s1 <- fa2[[from]]
        s2 <- fa2[[to]]
        sc1 <- fa2[[from_sc]]
        sc2 <- fa2[[to_sc]]
        out <- data.table(pos=1:nrow(fa2))
        out[s1=='-' | s2=='-',class:='Gap']
        out[is.na(class) & s1==s2,class:='Identical']
        out[is.na(class) & sc1==sc2,class:='Not-identical, same type']
        out[is.na(class),class:='Not-identical, different type']
        out$isoform1 <- from
        out$isoform2 <- to
        out
    }

    ## run for each unique pairwise combination of isoforms (each row of m)
    l <- lapply(1:nrow(m), f, m, fa2)

    ## return result combined into single data.table
    d <- rbindlist(l)
    d$comparison <- paste0(d$isoform1,':',d$isoform2)
    d
}


plot_pairwise_alignment_sidechain_comparisons <- function(d, fa2, highlighted) {
    cols <- rev(c('black',brewer.pal(9, name='RdBu')[c(3,5,9)]))
    names(cols) <- c('Identical','Not-identical, same type','Not-identical, different type','Gap')
    d$class <- factor(d$class, levels=(names(cols)))
    d$comparison <- factor(d$comparison, levels=rev(unique(d$comparison)))
    levs <- levels(d$comparison)
    d$comparison <- as.integer(d$comparison) - 0.5   
    fa2$label <- ''

    for(h in names(highlighted)) {
        residues <- highlighted[[h]]
        labels <- paste(h,residues)
        pos_field <- paste0(h,'_pos')
        fa2[which(fa2[[pos_field]] %in% residues),label:=paste(label,labels)]
    }

    p <- ggplot(d, aes(x=pos, y=comparison)) +
        scale_y_continuous(expand=c(0,0),limits=c(0,7),breaks=seq(0.5,by=1,length.out=length(levs)),labels=levs) + 
        geom_tile(aes(fill=class)) +
        geom_segment(data=fa2[label!='',], y=3, yend=3.1, aes(x=global_position, xend=global_position)) +
        geom_text(data=fa2[label!='',], y=3.15,aes(x=global_position, label=label),angle=90,hjust=0,vjust=0.5) +
        scale_fill_manual(values=cols,name='Amino acids') +
        ang::theme_ang(base_size=12) +
        theme(legend.position='bottom') +
        labs(x='Multiple sequence alignment position',y=NULL)
    if(nrow(fa2) > 400) {
        p <- p + scale_x_continuous(expand=c(0,0),breaks=seq(0,max(d$pos),by=100))
    } else {
        p <- p + scale_x_continuous(expand=c(0,0),breaks=seq(0,max(d$pos),by=25))   
    }

    p
}


# ~~~~~~~~~~~~~~~~~~~~~~~~
# run for PI3K isoforms
# ~~~~~~~~~~~~~~~~~~~~~~~~

## annotate each residue in the cobalt multiple-alignment, save the result
alignment_fa <- 'original_data/cobalt/cobalt_alignment_PI3K.fa'
aa_table <- 'original_data/ext/amino_acid_types.txt'
fa2 <- annotate_cobalt_alignment(alignment_fa, aa_table)
ang::write_tsv(fa2, 'processed_data/cobalt_alignment_PI3K_annotated.tsv')

## get pairwise comparisons of the aligned side-chain types for each pair of isoforms
d <- pairwise_comparison_of_sidechain_types(fa2, c('PIK3CA','PIK3CB','PIK3CD'))

## highlight these side positions in the global alignment
highlighted <- list(PIK3CD=c(38,124,416,1021))

## generate plot with global alignment indicating side-chain comparisons
p <- plot_pairwise_alignment_sidechain_comparisons(d, fa2, highlighted)
ggsave('figures/multiple_sequence_alignment_comparison_PI3K.pdf',width=8,height=5)


# ~~~~~~~~~~~~~~~~~~~~~~~~
# run for RAS isoforms
# ~~~~~~~~~~~~~~~~~~~~~~~~

## annotate each residue in the cobalt multiple-alignment, save the result
alignment_fa <- 'original_data/cobalt/cobalt_alignment_RAS.fa'
aa_table <- 'original_data/ext/amino_acid_types.txt'
fa2 <- annotate_cobalt_alignment(alignment_fa, aa_table)
ang::write_tsv(fa2, 'processed_data/cobalt_alignment_RAS_annotated.tsv')

## get pairwise comparisons of the aligned side-chain types for each pair of isoforms
d <- pairwise_comparison_of_sidechain_types(fa2, c('HRAS','KRAS','NRAS'))

## highlight these side positions in the global alignment
highlighted <- list(KRAS=c(12,61,146), NRAS=c(12,61,146), HRAS=c(12,61,146))

## generate plot with global alignment indicating side-chain comparisons
p <- plot_pairwise_alignment_sidechain_comparisons(d, fa2, highlighted)
ggsave('figures/multiple_sequence_alignment_comparison_RAS.pdf',width=8,height=5)








## sequence identity between A, B, D
proteins <- c('PIK3CA','PIK3CB','PIK3CD')
m <- matrix(nrow=3, ncol=3)
rownames(m) <- proteins; colnames(m) <- proteins
m[is.na(m)] <- 1
for(i in 1:nrow(m)) m[i,i] <- NA
m <- as.data.table(reshape2::melt(m))
m <- m[!is.na(value),]
f=function(i, m, fa2) {
    from=as.character(m$Var1[i])
    to=as.character(m$Var2[i])
    from_sc <- paste0(from,'_side_chain')
    to_sc <- paste0(to,'_side_chain')
    s1 <- fa2[[from]]
    s2 <- fa2[[to]]
    
    ## length of each protein individually
    n1 <- sum(s1!='-')
    n2 <- sum(s2!='-')

    ## length of alignment (number of residues without gaps in either)
    ind <- which(s1!='-' & s2!='-')
    aligned <- length(ind)

    ## length of alignment with same side chain-type
    sc1 <- fa2[[from_sc]]
    sc2 <- fa2[[to_sc]]
    ind_sametype <- which(s1!='-' & s2!='-' & sc1==sc2)
    sametype <- length(ind_sametype)

    ## length of alignment with same amino acid
    ind_identity <- which(s1!='-' & s2!='-' & s1==s2)
    identity <- length(ind_identity)

    list(v1=from, v2=to, n1=n1, n2=n2, aligned=aligned, identity=identity, sametype=sametype) 
}
l <- lapply(1:nrow(m), f, m, fa2)

d <- rbindlist(l)
d[,diff_sametype:=sametype-identity]
d[,diff_not_sametype:=aligned-sametype]
d[,c('sametype'):=NULL]
d2 <- melt(d, id.vars=c('v1','v2','n1','n2','aligned'))
d2[variable=='identity',variable:='Identical']
d2[variable=='diff_sametype',variable:='Not-identical, same type']
d2[variable=='diff_not_sametype',variable:='Not-identical, different type']
d2$variable <- factor(d2$variable, levels=rev(c('Identical','Not-identical, same type','Not-identical, different type')))
d2[,comparison:=paste0(v1,':',v2)]
d2 <- d2[comparison %in% levs,]
d2$comparison <- factor(d2$comparison, levels=levs[1:3])
p <- ggplot(d2, aes(x=comparison, y=value)) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    geom_bar(stat='identity',aes(fill=variable)) + 
    coord_flip() +
    ang::theme_ang(base_size=12) +
    theme(legend.position='none') +
    labs(y='N residues',x=NULL) +
    scale_fill_manual(values=cols[names(cols) %in% d2$variable], name='Aligned amino acids')
ggsave('figures/multiple_sequence_alignment_comparison_right.pdf',width=4,height=5)

d$prop_identity <- round(100*d$identity / d$aligned, 1)
d$prop_sametype <- round(100*(d$diff_sametype) / d$aligned, 1)
d$prop_difftype <- round(100*(d$diff_not_sametype) / d$aligned, 1)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compare how well-aligned region around important residues is
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

x <- data.table(pos=fa2$global_position, res_a=fa2$PIK3CA, res_d=fa2$PIK3CD, dpos=fa2$PIK3CD_pos,
                same_identity=as.integer(fa2$PIK3CA==fa2$PIK3CD & fa2$PIK3CD!='-' & fa2$PIK3CA!='-'),
                same_type=as.integer(fa2$PIK3CA_side_chain==fa2$PIK3CD_side_chain & fa2$PIK3CD!='-' & fa2$PIK3CA!='-'))
x <- x[res_d!='-']
x$avg <- as.numeric(NA)
n <- nrow(x)

k <- 10
for(i in 1:n) {
    positions <- (i - k) : (i + k)
    positions <- positions[positions %in% x$pos]
    vals <- x$same_type[positions]
    x$avg[i] <- mean(vals, na.rm=T)
}
x <- x[!is.na(avg),]
x$avg <- 100*x$avg

acivating <- c(1023,525,522,88,334,416,38,1025,414,1021,124)
notactivating <- c(806,367,894,971,338,911)
x$Colonies <- 'Unknown'
x[dpos %in% activating, Colonies:='Yes']
x[dpos %in% notactivating, Colonies:='No']

library(dunn.test)
tst <- dunn.test(x=x$avg, g=x$Colonies, kw=T, method='Holm')
pvals <- paste0('P=',prettyNum(tst$P.adjusted,digits=1))
f=function(s) strsplit(s,' - ')[[1]]
comparisons <- lapply(tst$comparisons, f)
x[dpos %in% highlight,label:=paste0(res_d,dpos)]

library(ggrepel)
library(ggsignif)
p <- ggplot(x, aes(Colonies,y=avg)) +
    scale_y_continuous(limits=c(0,120), breaks=seq(0,100,by=25)) + 
    geom_point(position=position_jitter(width=0.15,height=0,seed=42),pch=21,
               fill='#bfbfbf',color='black',stroke=0.25,size=2) +
    geom_text(aes(label=label),position=position_jitter(width=0.15,height=0,seed=42)) +
    geom_boxplot(fill=NA,width=0.8) +
    theme_ang(base_size=12) +
    geom_signif(comparison=comparisons, annotations=pvals ,tip_length=NA, y_position=c(110,105,100)) +
    labs(x='Enriched colonies when mutated',y=paste0('% residues of the\nsame type [pos-',k,', pos+',k,']'))
ggsave(paste0('figures/prop_same_type_boxplot_k',k,'.pdf'),width=4,height=5)




#x <- x[Colonies!='Unknown',]
a <- x$avg[x$Colonies=='Yes']
b <- x$avg[x$Colonies=='No']
plab <- paste0('P=',prettyNum(wilcox.test(a,b,alternative='two.sided')$p.value, digits=2))


tmp <- x[Colonies=='Yes']
tmp <- tmp[order(avg,decreasing=F),]
tmp[,tm:=paste0(res_d,dpos)]
tmp$tm <- factor(tmp$tm, levels=unique(tmp$tm))
p <- ggplot(tmp, aes(x=tm,y=avg)) +
    geom_bar(stat='identity') 






tmp <- melt(fa2[,c('global_position','PIK3CA','PIK3CB','PIK3CD'),with=F],id.vars=c('global_position'))
tmp_sidechain <- melt(fa2[,c('global_position','PIK3CA_side_chain','PIK3CB_side_chain','PIK3CD_side_chain'),with=F],id.vars=c('global_position'))
tmp_sidechain$variable <- gsub('_side_chain','',tmp_sidechain$variable)
tmp_pos <- melt(fa2[,c('global_position','PIK3CA_pos','PIK3CB_pos','PIK3CD_pos'),with=F],id.vars=c('global_position'))
tmp_pos$variable <- gsub('_pos','',tmp_pos$variable)
tmp <- merge(tmp, tmp_sidechain, by=c('global_position','variable'), all=T)
names(tmp) <- c('global_position','isoform','aa','type')
tmp <- merge(tmp, tmp_pos, by.x=c('global_position','isoform'), by.y=c('global_position','variable'), all=T)
setnames(tmp,'value','local_pos')

hs <- fread('data/aligned_residue_enrichment.txt')
hs <- hs[aligned_hotspot==T]
hs$pik3cd_pos <- as.integer(substr(hs$pik3cd_residue,9,nchar(hs$pik3cd_residue)))
hs$pik3ca_pos <- as.integer(substr(hs$pik3ca_residue,9,nchar(hs$pik3ca_residue)))
hotspots_global <- tmp$global_position[tmp$local_pos %in% hs$pik3cd_pos & tmp$isoform=='PIK3CD']
tmp[global_position %in% hotspots_global,focus:=global_position] 

plot_focus <- function(pos, tmp, k) {
    message(pos)
    tmp2 <- tmp[global_position >= pos-k & global_position <= pos+k,]
    get_range <- function(tmp2) {
        centered <- tmp2$local_pos[tmp2$global_position==pos]
        centered_aa <- tmp2$aa[tmp2$global_position==pos]
        list(centered_aa=centered_aa, centered=centered, start=min(tmp2$local_pos,na.rm=T),end=max(tmp2$local_pos,na.rm=T))
    }
    rs <- tmp2[,get_range(.SD),by=isoform]
    rs[,label:=paste0(isoform,' ',centered_aa,centered,' [',start,'-',end,']')]
    tmp2 <- merge(tmp2, rs, by.x='isoform', by.y='isoform', all.x=T)
    tmp2 <- tmp2[order(isoform,global_position),]
    tmp2$label <- factor(tmp2$label, levels=rev(unique(tmp2$label)))

    require(RColorBrewer)
    require(ggplot2)
    cols <- c(brewer.pal(6, 'Accent'),'black')
    names(cols) <- unique(aa$side_chain)
    p <- ggplot(tmp2, aes(x=global_position,y=label)) + 
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        geom_tile(data=tmp2[!global_position %in% pos], aes(fill=type)) +
        geom_tile(data=tmp2[global_position %in% pos], aes(fill=type), color='black', size=1) +
        scale_fill_manual(values=cols,name='Side chain') +
        ang::theme_ang(base_size=12) +
        geom_text(aes(label=aa)) +
        theme(axis.line=element_blank(), axis.ticks=element_blank(), legend.position='bottom') +
        labs(x='Global aligned position', y='Isoform')
    outfile <- paste0('figures/pi3k_alignment_',pos,'.pdf')
    ggsave(outfile,width=8,height=2)
}
l <- lapply(hotspots_global, plot_focus, tmp, k=10)


table_freq <- function (value) {
    if (is.null(value) | length(value) == 0) {
        tbl <- data.table(value = NA, N = NA)
    }
    else {
        tbl <- as.data.table(table(value))
        tbl <- tbl[order(tbl$N, decreasing = T), ]
    }
    tbl
}
tbl <- table_freq(fa$label)
tbl2 <- table_freq(fa$label[fa$sequence!='-'])




