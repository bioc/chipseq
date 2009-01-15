
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
## shifting.  For the best mean shift, this would be low.


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


simdf <-
    expand.grid(shift = seq(from = 0, to = 150, by = 5),
                chr = factor(sprintf("chr%s", 1:19), levels = sprintf("chr%s", 1:19)),
                sample = factor(c("cfibromyod", "cblasts", "ctubes"), levels = c("cfibromyod", "cblasts", "ctubes")),
                sim = NA_real_,
                KEEP.OUT.ATTRS = FALSE)

for (i in seq_len(nrow(simdf))[1:10])
{
    with(simdf, message(chr[i], "\t", shift[i], "\t", sample[i]))
    simdf$sim[i] <-
        computeCovered(combinedLanes[[as.character(simdf$sample[i])]],
                       chr = as.character(simdf$chr[i]),
                       shift = simdf$shift[i])
}



## A slightly different and more useful question is: how much should
## be extend.  This needs some way to evaluate an extension.  Here's
## one idea: for a given extension, take covpos + covneg, and
## pmax(covpos, covneg).  The difference in their sum() (area under
## coverage) is a measure of overlap.  



computeOverlap <-
    function(x, chr, seqLen = 100,
             chrlens = seqlengths(Mmusculus))
{
    xchr <- x[[chr]]
    n <- chrlens[chr]
    covpos <- coverage(extendReads(xchr, readLen = 35, seqLen = seqLen, strand = "+"), 1L, n)
    covneg <- coverage(extendReads(xchr, readLen = 35, seqLen = seqLen, strand = "-"), 1L, n)
    cov <- covpos + covneg
    c(total = sum(cov), diff = sum(cov - pmax(covpos, covneg)))
}


simdf$total <- NA_real_
simdf$diff <- NA_real_

for (i in seq_len(nrow(simdf))[1:10])
{
    with(simdf, message(chr[i], "\t", shift[i], "\t", sample[i]))
    ansi <-
        computeOverlap(combinedLanes[[as.character(simdf$sample[i])]],
                       chr = as.character(simdf$chr[i]),
                       seqLen = 35L + simdf$shift[i])
    simdf$total[i] <- ansi["total"]
    simdf$diff[i] <- ansi["diff"]
}

frag.size <- simdf

save(frag.size, file = "frag.size.rda")





pdf("strand-shift.pdf", width = 11, height = 8)

xyplot(sim ~ (shift+35) | chr, simdf, type = "o",
       subset = (sample == "cfibromyod"),
       scales = list(y = list(relation = "free", draw = FALSE)),
       panel = function(...) {
           panel.abline(v = 140, col = "grey", lwd = 3)
           panel.xyplot(...)
       })

xyplot(sim ~ (shift+35) | chr, simdf, type = "o",
       subset = (sample == "cblasts"),
       scales = list(y = list(relation = "free", draw = FALSE)),
       panel = function(...) {
           panel.abline(v = 90, col = "grey", lwd = 3)
           panel.xyplot(...)
       })

xyplot((diff/total) ~ (shift + 35) | chr, simdf, 
       subset = (sample == "cfibromyod"),
       ## scales = list(y = list(relation = "free", draw = FALSE)),
       type = c("l", "g"))


xyplot((diff) ~ (shift + 35) | chr, simdf, 
       subset = (sample == "cfibromyod"),
       scales = list(y = list(relation = "free", draw = FALSE)),
       type = c("l", "g"),
       prepanel = function(x, y, ...) {
           prepanel.default.xyplot(x[-1], diff(y), ...)
       },
       panel = function(x, y, ...) {
           panel.xyplot(x[-1], diff(y), ...)
       })




dev.off()




stop("End of file")

## older experiments


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

## for (i in seq_len(nrow(simdf)))
## {
##     with(simdf, message(chr[i], "\t", shift[i], "\t", sample[i]))
##     simdf$sim[i] <-
##         computeSimilarity(combinedLanes[[as.character(simdf$sample[i])]],
##                           chr = as.character(simdf$chr[i]),
##                           shift = simdf$shift[i],
##                           similarity = similarity.corr)
## }

xyp <-
    xyplot(sim ~ shift | chr + sample, simdf, type = "o",
           subset = (sample == "cfibromyod"),
           panel = function(...) {
               panel.abline(v = 60, col = "grey", lwd = 3)
               panel.xyplot(...)
           })




