
## following up on an idea in the Kharchenko et al paper, try to
## estimate average fragment length by strand shifting


library(chipseq)

library(BSgenome.Mmusculus.UCSC.mm9)

load("myodMyo.rda")
load("myodFibro.rda")
load("solexa54.rda")

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
         cfibromyod = combineLaneReads(myodFibro[c("2","4","7")]),
         methyl0 = combineLaneReads(solexa54[c("1","2")]),
         methyl96 = combineLaneReads(solexa54[c("3","4")]),
         realtube = combineLaneReads(solexa54[c("7","8")]))





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
    expand.grid(shift = seq(from = 0, to = 200, by = 5),
                chr = factor(sprintf("chr%s", 1:19), levels = sprintf("chr%s", 1:19)),
                sample = factor(rev(names(combinedLanes)), levels = rev(names(combinedLanes))),
                sim = NA_real_,
                KEEP.OUT.ATTRS = FALSE)

for (i in seq_len(nrow(simdf)))
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

for (i in seq_len(nrow(simdf)))
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



diff0Plot <- function(which = "cfibromyod")
{
    xyplot(diff/total ~ (shift + 35) | chr, frag.size, 
           subset = (sample == which), main = which, 
           scales = list(y = list(relation = "free", draw = FALSE)),
           type = c("l", "g"),
           xlab = "Amount of extension (bases)",
           ylab = "Difference between total and pairwise max coverage / Total")
}

diff1Plot <- function(which = "cfibromyod")
{
    xyplot(diff ~ (shift + 35) | chr, frag.size, 
           subset = (sample == which), main = which, 
           ## scales = list(y = list(relation = "free", draw = FALSE)),
           type = c("l", "g"),
           xlab = "Amount of extension (bases)",
           ylab = "First derivative of Difference",
           prepanel = function(x, y, ...) {
               prepanel.default.xyplot(head(x, -1), diff(y), ...)
           },
           panel = function(x, y, ...) {
               panel.xyplot(head(x, -1), diff(y), ...)
               panel.abline(h = 0)
           })
}

diff2Plot <- function(which = "cfibromyod")
{
    xyplot(diff ~ (shift + 35) | chr, frag.size, 
           subset = (sample == which), main = which, 
           scales = list(y = list(relation = "free", draw = FALSE)),
           type = c("l", "g"),
           xlab = "Amount of extension (bases)",
           ylab = "Second derivative of Difference",
           prepanel = function(x, y, ...) {
               prepanel.default.xyplot(head(x, -2), diff(diff(y)), ...)
           },
           panel = function(x, y, ...) {
               panel.xyplot(head(x, -2), diff(diff(y)), ...)
               panel.abline(h = 0)
           })
}


coveredPlot <- function(which = "cfibromyod")
{
    xyplot(sim ~ (shift+35) | chr, frag.size, type = "o",
           subset = (sample == which), main = which, 
           scales = list(y = list(relation = "free", draw = FALSE)),
           xlab = "Estimated mean fragment size",
           ylab = "Second derivative of Difference",
           panel = function(...) {
               panel.abline(v = 140, col = "grey", lwd = 3)
               panel.xyplot(...)
           })
}



pdf("strand-shift.pdf", width = 11, height = 8)


coveredPlot("cfibromyod")
coveredPlot("cblasts")
coveredPlot("ctubes")
coveredPlot("cfibro")
coveredPlot("methyl0")
coveredPlot("methyl96")
coveredPlot("realtube")

diff0Plot("cfibromyod")
diff0Plot("cblasts")
diff0Plot("ctubes")
diff0Plot("cfibro")
diff0Plot("methyl0")
diff0Plot("methyl96")
diff0Plot("realtube")

diff1Plot("cfibromyod")
diff1Plot("cblasts")
diff1Plot("ctubes")
diff1Plot("cfibro")
diff1Plot("methyl0")
diff1Plot("methyl96")
diff1Plot("realtube")

diff2Plot("cfibromyod")
diff2Plot("cblasts")
diff2Plot("ctubes")
diff2Plot("cfibro")
diff2Plot("methyl0")
diff2Plot("methyl96")
diff2Plot("realtube")


dev.off()




stop("End of file")


## A demo of this

frag.mu <- 150
frag.sd <- 5

peak.locs <- c(0, 50, 350)
nreads <- 150
strand <- sample(c(-1, 1), nreads, replace = TRUE)
peak <- sample(peak.locs, nreads, replace = TRUE)

x <- round(peak - strand * runif(nreads, min = 0, pmax(0, rnorm(nreads, mean = frag.mu, sd = frag.sd)))) - ifelse(strand > 0, 0, 35)

stripplot(factor(strand) ~ x, jitter = TRUE,
          panel = function(...) {
              panel.abline(v = peak.locs)
              panel.stripplot(...)
          })

reads <- split(x, strand)
names(reads) <- c("-", "+")



plotOverlap <- function(seqLen = 100, main = as.character(seqLen), ...)
{
    rng <- as.integer(range(reads)) + c(-300L, 300L)
    cov.pos <- coverage(extendReads(reads, readLen = 35, seqLen = seqLen, strand = "+"), rng[1], rng[2])
    cov.neg <- coverage(extendReads(reads, readLen = 35, seqLen = seqLen, strand = "-"), rng[1], rng[2])
    cov.total <- cov.pos + cov.neg
    cov.max <- pmax(cov.pos, cov.neg)
    df <- data.frame(coverage = c(as.numeric(cov.total), as.numeric(cov.max)),
                     location = seq(from = rng[1], to = rng[2]),
                     which = gl(2, length(cov.max), labels = c("total", "pmax")))
    if (require(mosaiq))
    {
        print(mosaiq.xyplot(coverage ~ location, data = df, groups = which,
                            type = "l", main = main, ...))
    }
    c(total = sum(cov.total), diff = sum(cov.total) - sum(cov.max), depth = max(cov.total))
}

foo <- data.frame(seqlen = seq(10, 400, by = 5))
foo <- cbind(foo, t(sapply(foo$seqlen, plotOverlap)))

xyplot(diff + total + (diff/total) + depth ~ seqlen, foo,
       outer = TRUE, scales = list(y = list(relation = "free", rot = 0)),
       panel = function(...) {
           panel.xyplot(...)
           panel.abline(v = frag.mu + c(-2, 0, 2) * frag.sd)
       },
       type = c("l", "g"), as.table = TRUE, layout = c(1, 4))

xyplot(diff(diff) ~ seqlen[-1], foo,
       panel = function(...) {
           panel.xyplot(...)
           panel.abline(v = frag.mu + c(-2, 0, 2) * frag.sd)
       },
       type = c("l", "g"))


xyplot(diff(diff(diff)) ~ seqlen[-(1:2)], foo,
       panel = function(...) {
           panel.xyplot(...)
       },
       type = c("l", "g"))



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




