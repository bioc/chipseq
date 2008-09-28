
library("lattice")
library("chipseq")
library("geneplotter")

if (file.exists("alignInfo.rda")) load("alignInfo.rda") else
{

    library("BSgenome.Mmusculus.UCSC.mm9")

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


## basically same, but retains order of chromosomes
system.time(seqRanges <- lapply(lreads, growSeqs), gcFirst=TRUE)

cblasts = combineLanes(seqRanges[c(1,3,5)])
ctubes = combineLanes(seqRanges[c(2,4,6)])

peakSummary.blasts.wrt.tubes <-
    diffPeakSummary(obs.ranges = cblasts, ref.ranges = ctubes,
                    chrom.lens = chromLens, lower = 10)



## Robust regression (note that lmrob() in robustbase doesn't seem to
## handle weights).  The log2 responses seem reasoanbly homoskedastic.
## Our model of y ~ c.x translates to log(y) ~ log(c) + log(x), so all
## we need to do is estimate log(c). A robust estimate is
## median(log(y)-log(x)).  So residuals are
##  [ (log(y)-log(x)) - median(log(y)-log(x)) ]


peakSummary.rob <- 
    within(peakSummary.blasts.wrt.tubes,
       {
           diffs <- log2(obs.sums)-log2(ref.sums)
           resids <-  # prefer per-chromosome?
               (diffs - median(diffs)) / mad(diffs)
           resids.mean <- diffs - mean(diffs[is.finite(diffs)])
       })

xyplot(log2(obs.sums) ~ log2(ref.sums) | chromosome,
       data = peakSummary.blasts.wrt.tubes, auto.key = TRUE,
       subset = (chromosome %in% c("chr1", "chr2", "chr3", "chr4") &
                 obs.sums > 0),
       panel = function(...) {
           panel.smoothScatter(...)
           ## panel.xyplot(...)
           panel.abline(median(diffs), 1)
       },
       par.settings = simpleTheme(pch = ".", cex = 3),
       type = c("p", "g", "r"), col.line = "black", aspect = "iso")

bwplot(chromosome ~ resids, data = peakSummary.rob)

toppeaks <- subset(peakSummary.rob, abs(resids) > 4)
rownames(toppeaks) <- NULL
toppeaks[rev(order(abs(toppeaks$resids))), 1:6]



