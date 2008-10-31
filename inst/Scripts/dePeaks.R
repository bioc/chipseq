
library("lattice")
library("chipseq")
library("geneplotter")
   library("BSgenome.Mmusculus.UCSC.mm9")

if( FALSE) {
if (file.exists("alignInfo.rda")) load("alignInfo.rda") else
{

 
    ##lanes 1, 3, 6 are Myoblasts
    ##lanes 2, 4, 7 are Myotubes
    ##lane 8 is a reference lane

    lanes <- c(1, 2, 3, 4, 6, 7, 8)

    reads = vector("list", length = length(lanes))
    names(reads) = as.character(lanes)

    for (i in seq_along(lanes) ) {
        lane <- lanes[i]
        message("Starting Lane ", lane)
        pat <- paste("s_", lane, ".map", sep="")
        ## we drop the sex chromosomes and mitochondria.
        reads[[i]] <- readAndClean("/home/jdavison/ycao/26-06-2008/binary",
                                   pattern = pat, exclude = "[MXY]|rand")
    }

    lreads <- lapply(reads, as.list) # same info (no quality etc) as nested list

    ## mouse has 19 chromosomes

    chrom.list <- paste("chr", c(1:19), sep = "")

    ## nchrom <- length(chrom.list)
    ## chromLens = rep(NA, nchrom)
    ## names(chromLens) = chrom.list
    ## for( i in 1:nchrom) 
    ##     chromLens[i] = nchar(unmasked(Mmusculus[[chrom.list[i]]])) 

    chromLens <-
        sapply(chrom.list,
               function(chr) {
                   nchar(unmasked(Mmusculus[[chr]]))
                   ## same as 'length(Mmusculus[[chr]])' ?
               },
               simplify = TRUE)


    save(lreads, chromLens, file = "alignInfo.rda")

    ## system.time(seqRanges.old <- lapply(reads, growSeqs), gcFirst=TRUE)
}
}

load("myodMyo.rda")

chromLens <- seqlengths(Mmusculus)

seqRanges <- lapply(myodMyo, growSeqs)

cblasts = combineLanes(seqRanges[c(1,3,5)])
ctubes = combineLanes(seqRanges[c(2,4,6)])

peakSummary.blasts.wrt.tubes <-
    diffPeakSummary(obs.ranges = cblasts, ref.ranges = ctubes,
                    chrom.lens = chromLens, lower = 10)

peakSummary.tubes.wrt.blasts <-
    diffPeakSummary(obs.ranges = ctubes, ref.ranges = cblasts,
                    chrom.lens = chromLens, lower = 10)



## Robust regression (note that lmrob() in robustbase doesn't seem to
## handle weights).  The log2 responses seem reasoanbly homoskedastic.
## Our model of y ~ c.x translates to log(y) ~ log(c) + log(x), so all
## we need to do is estimate log(c). A robust estimate is
## median(log(y)-log(x)).  So residuals are
##  [ (log(y)-log(x)) - median(log(y)-log(x)) ]


peakSummary.blasts.wrt.tubes <- 
    within(peakSummary.blasts.wrt.tubes,
       {
           diffs <- log2(obs.sums)-log2(ref.sums)
           resids <-  # prefer per-chromosome?
               (diffs - median(diffs)) / mad(diffs)
           resids.mean <- diffs - mean(diffs[is.finite(diffs)])
       })

peakSummary.tubes.wrt.blasts <- 
    within(peakSummary.tubes.wrt.blasts,
       {
           diffs <- log2(obs.sums)-log2(ref.sums)
           resids <-  # prefer per-chromosome?
               (diffs - median(diffs)) / mad(diffs)
           resids.mean <- diffs - mean(diffs[is.finite(diffs)])
       })

xyplot(log2(obs.sums) ~ log2(ref.sums) | chromosome,
       data = peakSummary.tubes.wrt.blasts, auto.key = TRUE,
       subset = (chromosome %in% c("chr1", "chr2") &
                 obs.sums > 0),
       panel = function(x, y, ...) {
           panel.smoothScatter(x, y, ...)
           ## panel.xyplot(x, y, ...)
           panel.abline(median(y - x), 1)
       },
       par.settings = simpleTheme(pch = ".", cex = 3),
       type = c("p", "g", "r"), col.line = "black", aspect = "iso")

bwplot(chromosome ~ resids, data = peakSummary.blasts.wrt.tubes)
bwplot(chromosome ~ resids, data = peakSummary.tubes.wrt.blasts)

toppeaks <- subset(peakSummary.tubes.wrt.blasts, abs(resids) > 4)
rownames(toppeaks) <- NULL
toppeaks[rev(order(abs(toppeaks$resids))), ]

save(peakSummary.blasts.wrt.tubes, peakSummary.tubes.wrt.blasts,
     file = "peakSummary.rda")

if (FALSE)
{

    load("peakSummary.rda")

    peakSummary.blasts.wrt.tubes <- 
        peakSummary.blasts.wrt.tubes[c("chromosome", "start", "end",
                                       "ref.sums", "obs.sums",
                                       "ref.maxs", "obs.maxs")]
    rownames(peakSummary.blasts.wrt.tubes) <- NULL
    peakSummary.tubes.wrt.blasts <- 
        peakSummary.tubes.wrt.blasts[c("chromosome", "start", "end",
                                       "ref.sums", "obs.sums",
                                       "ref.maxs", "obs.maxs")]
    rownames(peakSummary.tubes.wrt.blasts) <- NULL

    peakSummary.blasts.wrt.tubes <- 
        within(peakSummary.blasts.wrt.tubes,
           {
               diffs <- log2(obs.sums)-log2(ref.sums)
               resids <-  (diffs - median(diffs)) / mad(diffs)
               rm(diffs)
           })

    peakSummary.tubes.wrt.blasts <- 
        within(peakSummary.tubes.wrt.blasts,
           {
               diffs <- log2(obs.sums)-log2(ref.sums)
               resids <-  (diffs - median(diffs)) / mad(diffs)
               rm(diffs)
           })

    write.csv(peakSummary.blasts.wrt.tubes, file="peakSummary_blasts_wrt_tubes.csv", quote=FALSE, row.names = FALSE)
    write.csv(peakSummary.tubes.wrt.blasts, file="peakSummary_tubes_wrt_blasts.csv", quote=FALSE, row.names = FALSE)
    
}

