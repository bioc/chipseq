
## following up on an idea in the Kharchenko et al paper, try to
## estimate average fragment length by strand shifting


library(chipseq)

library(BSgenome.Mmusculus.UCSC.mm9)

load("myodMyo.rda")
load("myodFibro.rda")


gsample <- function(x, ...)
{
    if (length(x) == 0) x
    else sample(x, ...)
}

combineLaneReads <- function(laneList, chromList = names(laneList[[1]])) {
    names(chromList) = chromList ##to get the return value named
    lapply(chromList,
           function(chr) { ## sample() to make order random
               list("+" = gsample(unlist(lapply(laneList, function(x) x[[chr]][["+"]]), use.names = FALSE)),
                    "-" = gsample(unlist(lapply(laneList, function(x) x[[chr]][["-"]]), use.names = FALSE)))
           })
}

combinedLanes <- 
    list(cblasts = combineLaneReads(myodMyo[c("1","3","6")]),
         ctubes = combineLaneReads(myodMyo[c("2","4","7")]),
         cfibro = combineLaneReads(myodFibro[c("1","3","6")]),
         cfibromyod = combineLaneReads(myodFibro[c("2","4","7")]))



## some measure of similarity.  Works on data from a single
## chromosome, in the form of a list("+" = ..., "-" = ...)


## latest attempt: measure goodness by number of bases covered after
## shifting.  Best match would have this be low.


computeCovered <-
    function(x, chr, shift = 0, 
             chrlens = seqlengths(Mmusculus))
{
    xchr <- x[[chr]]
    xchr[["-"]] <- xchr[["-"]] - shift
    n <- chrlens[chr]
    cov <- coverage(extendReads(xchr, readLen = 35, seqLen = 35),
                    1, n)
    s <- slice(cov, lower = 1)
    sum(as.numeric(width(s)))
}





sum.xrle <- function(x)
{
    sum(as.numeric(as.integer(x@values)) * as.numeric(as.integer(x@lengths)))
}

setMethod("sum", "XRleInteger",
          function(x, ..., na.rm = FALSE) sum.xrle(x))



## cross-correlation

similarity.cc <- function(pos, neg)
{
    ## pos, neg are external-vector integers
    sum(pos * neg)
}


## a version of uncentered correlation.  Other thoughts: contingency
## table indices

similarity.corr <- function(pos, neg)
{
    ## pos, neg are external-vector integers
    sum(pos * neg) / sqrt(sum(pos * pos) * sum(neg * neg))
}



computeSimilarity <-
    function(x, chr, shift = 0, chrlens = seqlengths(Mmusculus),
             similarity)
{
    xchr <- x[[chr]]
    xchr[["-"]] <- xchr[["-"]] - shift
    n <- chrlens[chr]
    covpos <-
        coverage(extendReads(xchr, readLen = 35, seqLen = 35, strand = "+"),
                 1, n)
    covneg <-
        coverage(extendReads(xchr, readLen = 35, seqLen = 35, strand = "-"),
                 1, n)
##     ans <- numeric(100)
##     for (i in seq_along(ans))
##         ans[i] <-
##             similarity(subseq(covpos, start = 1, end = n - i + 1),
##                        subseq(covneg, start = i, end = n))
##     ans
    similarity(covpos, covneg)
}

simdf <-
    expand.grid(shift = seq(from = 0, to = 150, by = 5),
                chr = factor(sprintf("chr%s", 1:19), levels = sprintf("chr%s", 1:19)),
                sample = factor(c("cfibromyod", "cblasts", "ctubes"), levels = c("cfibromyod", "cblasts", "ctubes")),
                sim = NA_real_,
                KEEP.OUT.ATTRS = FALSE)

## for (i in seq_len(nrow(simdf)))
## {
##     with(simdf, message(chr[i], "\t", shift[i], "\t", sample[i]))
##     simdf$sim[i] <-
##         computeSimilarity(combinedLanes[[as.character(simdf$sample[i])]],
##                           chr = as.character(simdf$chr[i]),
##                           shift = simdf$shift[i],
##                           similarity = similarity.corr)
## }

for (i in seq_len(nrow(simdf)))
{
    with(simdf, message(chr[i], "\t", shift[i], "\t", sample[i]))
    simdf$sim[i] <-
        computeCovered(combinedLanes[[as.character(simdf$sample[i])]],
                       chr = as.character(simdf$chr[i]),
                       shift = simdf$shift[i])
}

library(latticeExtra)

xyp <-
    xyplot(sim ~ shift | chr + sample, simdf, type = "o",
           subset = (sample == "cfibromyod"),
           panel = function(...) {
               panel.abline(v = 60, col = "grey", lwd = 3)
               panel.xyplot(...)
           })

xyplot(sim ~ shift | chr, simdf, type = "o",
       subset = (sample == "cfibromyod"),
       scales = list(relation = "free", draw = FALSE),
       panel = function(...) {
           panel.abline(v = 105, col = "grey", lwd = 3)
           panel.xyplot(...)
       })


pdf("strand-shift.pdf", width = 11, height = 8)
useOuterStrips(xyp[1:5, ])
useOuterStrips(xyp[6:10, ])
useOuterStrips(xyp[11:15, ])
useOuterStrips(xyp[16:19, ])
dev.off()



