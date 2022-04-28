
library(data.table)
a <- fread('data/aaw2872-data-s2.txt')
setnames(a,'d','relrate_exp')

d <- fread('data/output.PIK3CD.txt')
d <- d[Analysis_Type=='pancan',]
check_ref_aa <- function(s) strsplit(s,'[:]')[[1]][1]
d$Ref_AA_cleaned <- sapply(d$Ref_AA, check_ref_aa)
d <- d[nchar(Ref_AA_cleaned)==1]

parse_position <- function(i, d) {
    pos <- d$Genomic_Position[i]
    tm <- d$tm[i]
    s <- strsplit(pos,'[|]')[[1]]
    parse2 <- function(s) {
        str <- strsplit(s,'[:]')[[1]]
        chr <- str[1]
        str <- strsplit(str[2],'_')[[1]]
        pos <- str[1]
        n <- str[2]
        list(chr=chr, pos=pos, n=n)
    }
    l <- lapply(s, parse2)
    out <- rbindlist(l)
    out$Genomic_Position <- pos
    out$tm <- tm
    out
}
l <- lapply(1:nrow(d), parse_position, d)
dd <- rbindlist(l)
dd$chr <- as.integer(dd$chr)
dd$pos <- as.integer(dd$pos)
dd$n <- as.integer(dd$n)

dd <- merge(dd, a, by=c('chr','pos'), all.x=T)
out <- dd[!is.na(relrate_exp)]
out <- out[,c('chr','pos','tm','relrate_exp')]




