library(rpkg)
#library(ParalogousResidues)

alleles <- c('R88C','G124D','R894Q','C416R','E1021K','A835T','N334K')
d <- data.table(allele=alleles)
d$Reference_Amino_Acid <- substr(d$allele,1,1)
d$Amino_Acid_Position <- as.integer(substr(d$allele,2,(nchar(d$allele)-1)))
d$Variant_Amino_Acid <- substr(d$allele,nchar(d$allele),nchar(d$allele))
d$tm <- paste('PIK3CD',d$allele)
d <- d[order(Amino_Acid_Position),]
d$query_start <- d$Amino_Acid_Position - 15
d$query_end <- d$Amino_Acid_Position + 15



## start here
al <- read.table('G56XA7RA211-alignment.fa',sep='\n')
grep('>',al$V1)
pd <- paste(al$V1[2:20],collapse='')
pb <- paste(al$V1[22:40],collapse='')
pa <- paste(al$V1[42:60],collapse='')

al <- data.table(pik3cd=strsplit(pd,'')[[1]], pik3cb=strsplit(pb,'')[[1]], pik3ca=strsplit(pa,'')[[1]])
al$pik3cd_pos <- NA
al$pik3cd_pos[al$pik3cd!='-'] <- 1:sum(al$pik3cd!='-')
al$pik3cb_pos <- NA
al$pik3cb_pos[al$pik3cb!='-'] <- 1:sum(al$pik3cb!='-')
al$pik3ca_pos <- NA
al$pik3ca_pos[al$pik3ca!='-'] <- 1:sum(al$pik3ca!='-')
al$pik3cd_tm <- paste0('PIK3CD ',al$pik3cd,al$pik3cd_pos)
al[grepl('-NA',pik3cd_tm),pik3cd_tm:='']
al$i <- 1:nrow(al)

hs <- fread('~/luna/ext/resources/hotspots_20170430.tsv')
hs <- hs[Hugo_Symbol %in% c('PIK3CD','PIK3CB','PIK3CA')]
hs <- hs[,c('tm','class'),with=F]
hs <- hs[class=='24k',]

data(BLOSUM62)
b <- adt(melt(BLOSUM62))
names(b) <- c('aa1','aa2','score')

al <- merge(al, b, by.x=c('pik3cd','pik3cb'), by.y=c('aa1','aa2'), all.x=T)
setnames(al,'score','DB_score')
al <- merge(al, b, by.x=c('pik3cd','pik3ca'), by.y=c('aa1','aa2'), all.x=T)
setnames(al,'score','DA_score')
al <- merge(al, b[aa1==aa2,], by.x=c('pik3cd'), by.y=c('aa1'), all.x=T)
setnames(al,'score','DD_score')
al <- al[order(i,decreasing=F),]

summarize_hotspot <- function(i, d, al) {
    allele <- d$allele[i]
    current_range <- seq(d$query_start[i],d$query_end[i])
    tmp <- al[al$pik3cd_pos %in% current_range,]
    tmp <- tmp[order(pik3cd_pos),]
    tmp$pik3cb_tm <- paste0('PIK3CB ',tmp$pik3cb,tmp$pik3cb_pos)
    tmp$pik3ca_tm <- paste0('PIK3CA ',tmp$pik3ca,tmp$pik3ca_pos)
    tmp$pik3cd_hs <- tmp$pik3cd_tm %in% hs$tm    
    tmp$pik3cb_hs <- tmp$pik3cb_tm %in% hs$tm    
    tmp$pik3ca_hs <- tmp$pik3ca_tm %in% hs$tm    
    tmp$allele <- allele
    tmp$x <- 1:nrow(tmp)
    tmp
}

l <- lapply(1:nrow(d), summarize_hotspot, d, al)
ll <- rbindlist(l)
d1 <- melt(ll[,c('x','allele','pik3cd','pik3cb','pik3ca'),with=F], id.vars=c('x','allele'))
d2 <- melt(ll[,c('x','allele','DD_score','DB_score','DA_score'),with=F], id.vars=c('x','allele'))
d3 <- melt(ll[,c('x','allele','pik3cd_hs','pik3cb_hs','pik3ca_hs'),with=F], id.vars=c('x','allele'))

dd <- d1
setnames(dd,'value','AminoAcid')
setnames(dd,'variable','Gene')
dd$score <- d2$value
dd$hotspot <- d3$value
setnames(dd,'x','pos')
dd$Gene <- toupper(dd$Gene)
dd$Gene <- factor(dd$Gene, levels=rev(c('PIK3CD','PIK3CB','PIK3CA')))
dd$allele <- paste('PIK3CD',dd$allele)
dd$allele <- factor(dd$allele, unique(dd$allele))
dd$`Hotspot residue` <- 'No'
dd[hotspot==T,`Hotspot residue`:='Yes']

cols <- c('#B2182B','white','#79a3cd')
hscols <- c('black','#B2182B'); names(hscols) <- c('No','Yes')
p <- ggplot(dd, aes(x=Gene, y=pos)) +
    geom_tile(aes(fill=score)) +
    scale_fill_gradient2(low=cols[1],mid=cols[2],high=cols[3],midpoint=0,name='BLOSUM62\nsimilarity score') + 
    geom_text(aes(label=AminoAcid,color=`Hotspot residue`)) +
    scale_color_manual(values=hscols) + 
    facet_wrap(facets=~allele, ncol=1, nrow=7) +
    theme_std(base_size=14) +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), axis.ticks.x=element_blank(),
        axis.line.y=element_blank(), axis.ticks.y=element_blank()) + 
    labs(x=NULL,y=NULL,subtitle='PIK3CD hotspot paralogy') +
    coord_flip() 
pplot(p)
ggsave('pik3cd_hotspot_paralogy_figure.pdf',width=8,height=8)


geneinfo <- function(ll) {
    min_d <- min(ll$pik3cd_pos,na.rm=T)
    max_d <- max(ll$pik3cd_pos,na.rm=T)
    min_b <- min(ll$pik3cb_pos,na.rm=T)
    max_b <- max(ll$pik3cb_pos,na.rm=T)
    min_a <- min(ll$pik3ca_pos,na.rm=T)
    max_a <- max(ll$pik3ca_pos,na.rm=T)
    list(
         min_d=min_d,max_d=max_d,
         min_b=min_b,max_b=max_b,
         min_a=min_a,max_a=max_a
         )
}

info2 <- ll[,geneinfo(.SD),by=allele]







blastp <- fread('pik3cd_blastp_results_human_proteincoding.txt')
blastp$aa_length <- nchar(blastp$hit_sequence)
blastp[,c('hit_symbols','aa_length'),with=F]

pik3cd <- blastp$hit_sequence[blastp$hit_symbols=='PIK3CD']
d$query_start <- d$Amino_Acid_Position-15
d$query_end <- d$Amino_Acid_Position+15
d$query_seq <- ''
for(i in 1:nrow(d)) d$query_seq[i] <- substr(pik3cd,d$query_start[i],d$query_end[i])
write.tsv(d,'pik3cd_blastp_results_human_proteincoding_annotated.txt')





toadd <- data.table(uniprot_ID='Q9UBF8',Transcript_ID='',Hugo_Symbol='PI4KB',length=816,name='',
           Sequence='MGDTVVEPAPLKPTSEPTSGPPGNNGGSLLSVITEGVGELSVIDPEVAQKACQEVLEKVKLLHGGVAVSSRGTPLELVNGDGVDSEIRCLDDPPAQIREEEDEMGAAVASGTAKGARRRRQNNSAKQSWLLRLFESKLFDISMAISYLYNSKEPGVQAYIGNRLFCFRNEDVDFYLPQLLNMYIHMDEDVGDAIKPYIVHRCRQSINFSLQCALLLGAYSSDMHISTQRHSRGTKLRKLILSDELKPAHRKRELPSLSPAPDTGLSPSKRTHQRSKSDATASISLSSNLKRTASNPKVENEDEELSSSTESIDNSFSSPVRLAPEREFIKSLMAIGKRLATLPTKEQKTQRLISELSLLNHKLPARVWLPTAGFDHHVVRVPHTQAVVLNSKDKAPYLIYVEVLECENFDTTSVPARIPENRIRSTRSVENLPECGITHEQRAGSFSTVPNYDNDDEAWSVDDIGELQVELPEVHTNSCDNISQFSVDSITSQESKEPVFIAAGDIRRRLSEQLAHTPTAFKRDPEDPSAVALKEPWQEKVRRIREGSPYGHLPNWRLLSVIVKCGDDLRQELLAFQVLKQLQSIWEQERVPLWIKPYKILVISADSGMIEPVVNAVSIHQVKKQSQLSLLDYFLQEHGSYTTEAFLSAQRNFVQSCAGYCLVCYLLQVKDRHNGNILLDAEGHIIHIDFGFILSSSPRNLGFETSAFKLTTEFVDVMGGLDGDMFNYYKMLMLQGLIAARKHMDKVVQIVEIMQQGSQLPCFHGSSTIRNLKERFHMSMTEEQLQLLVEQMVDGSMRSITTKLYDGFQYLTNGIM',
            id='PI4KB'
           )
peptide_sequences <- rbind(peptide_sequences,toadd)

f <- function(symbol) {
    s <- strsplit(symbol,';')[[1]]
    out <- data.table(symbol=symbol,alias=s)
    out$have_tx <- out$alias %in% peptide_sequences$Hugo_Symbol
    out
}
l <- lapply(blastp$hit_symbol,f)
l <- rbindlist(l)
l <- l[have_tx==T,]

## get PRs for hotspots in each of these genes (published 24K hotspots only for now)
run_for_gene <- function(hit_gene) { 
    run_for_residue <- function(position, pattern_gene, hit_gene) { 
        message(pattern_gene,' residue=',position)
        result <- getPR(pattern_gene=pattern_gene,pattern_pos=position,subject_gene=hit_gene,cores=1,padding=15)
        result <- result[order(evalue,decreasing=F),]
        result 
    }
    l <- lapply(d$Amino_Acid_Position, run_for_residue, 'PIK3CD', hit_gene)
    results <- rbindlist(l)
    results
}

l2 <- mclapply(l$alias, run_for_gene, mc.cores=12)
l2 <- rbindlist(l2)




results <- results[evalue < 1e-5,]

bm <- adt(melt(BLOSUM62,id.var='-'))
names(bm) <- c('aa1','aa2','value')

results <- merge(results, bm, by.x=c('pattern_aa','subject_aa'),by.y=c('aa1','aa2'),all.x=T)
results$subject_tm <- paste0(results$subject_gene,' ',results$subject_aa,results$subject_pos)
results <- merge(results, hotspots[,c('tm','class'),with=F], by.x='subject_tm', by.y='tm', all.x=T)
write.tsv(results,'~/luna/projects/collaborations/arijh/PIK3CD/hotspot_sequence_paralogy.txt')








