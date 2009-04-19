


library(lattice)
library(chipseq)
library(chipseqData)
library(BSgenome.Mmusculus.UCSC.mm9)

data(myodMyo)
data(myodFibro)
data(pairedReads)

data(CpG.mm9)
mouse.chromlens <- seqlengths(Mmusculus)

data(geneMouse)
gregions <- genomic_regions(genes = geneMouse, proximal = 1000)
gregions <- subset(gregions, chrom %in% paste("chr", 1:19, sep = ""))
gregions$chrom <- gregions$chrom[drop = TRUE]

gpromoters <- gregions[c("chrom", "promoter.start", "promoter.end")]
names(gpromoters) <- c("chr", "start", "end")
gpromoters.split <- split(gpromoters, gpromoters$chr)


computeOverlap <- function(chr, start, end, tchr, tstart, tend)
{
    tsplit <- split(IRanges(tstart, tend), tchr)
    ans <- logical(length(chr))
    for (chrom in names(tsplit))
    {
        id <- which(chr == chrom)
        ans[id] <-
            !is.na(overlap(tsplit[[chrom]], 
                           IRanges(start[id], end[id]),
                           multiple = FALSE))
    }
    ans
}    

sqrt.scale.comps <- function (axis = c("x", "y"))
{
    axis <- match.arg(axis)
    switch(axis, x = function(...) {
        ans <- xscale.components.default(...)
        ans$bottom$labels$labels <- (ans$bottom$labels$at)^2
        ans
    }, y = function(...) {
        ans <- yscale.components.default(...)
        ans$left$labels$labels <- (ans$left$labels$at)^2
        ans
    })
}



ctubes <- combineLaneReads(myodMyo[c("2","4","7")])
cblast <- combineLaneReads(myodMyo[c("1","3","6")])
ptubes <- pairedReads[["2"]]

ctubes.ext <- gdapply(ctubes, extendReads, seqLen = 200)
cblast.ext <- gdapply(cblast, extendReads, seqLen = 200)
ptubes.ext <- gdapply(ptubes, extendReads, seqLen = 200)



peakSummaryTubeBlast <-
    diffPeakSummary(ctubes.ext, cblast.ext,
                    chrom.lens = mouse.chromlens,
                    lower = 15, islands = FALSE, merge = 20L)

## with(peakSummaryTubeBlast,
##      splom(data.frame(asinh(overlap1), asinh(maxs1), asinh(sums1),
##                       asinh(overlap2), asinh(maxs2), asinh(sums2)),
##            panel = panel.smoothScatter))

## with(peakSummaryTubeBlast,
##  {
##      xyplot(z ~ log2(comb.overlap) | CpG + promoter,
##             panel = function(...) {
##                 panel.smoothScatter(...)
##                 panel.grid(h = -1, v = -1)
##                 panel.loess(..., col = "yellow")
##             })
##  })
     


peakSummaryTubeBlast <-
    within(peakSummaryTubeBlast,
       {
           promoter <-
               ifelse(computeOverlap(chromosome, start, end,
                                     gpromoters$chr, 
                                     gpromoters$start, 
                                     gpromoters$end),
                      "In promoter", "Not in promoter")
           CpG <-
               ifelse(computeOverlap(chromosome, start, end,
                                     CpG.mm9$chr, 
                                     CpG.mm9$start, 
                                     CpG.mm9$end),
                      "In CpG island", "Not in CpG island")
           comb.overlap <- overlap1 + overlap2
           p <- sum(overlap2) / sum(comb.overlap)
           ## p <- median(overlap2 / comb.overlap)
           z <- ((overlap2 / comb.overlap) - p) / sqrt(p * (1-p) / comb.overlap)
           reg.z <- 
               local({
                   d <- asinh(sums2) - asinh(sums1)
                   ## dkeep <- d[CpG == "In CpG island"]
                   dkeep <- d
                   mu <- print(median(dkeep))
                   sigma <- print(mad(dkeep))
                   (d-mu)/sigma
               })
       })

## concordant?
with(peakSummaryTubeBlast, xtabs(~ (abs(z) > 2) + (abs(reg.z) > 2) + CpG))

## differences are related to coverage
subset(peakSummaryTubeBlast,
       (z > 5) & (reg.z < 2) & (CpG == "In CpG island"))[4:10]

pdf("reg-binom-comparison.pdf", width = 11, height = 8)

with(peakSummaryTubeBlast, 
 {
     comb.depth <- cut(comb.max, breaks = quantile(comb.max, c(0, 0.25, 0.5, 0.75, 1)))
     xyplot(z ~ reg.z | comb.depth + CpG,
            panel = function(...) {
                panel.smoothScatter(...)
                panel.abline(0, 1)
            },
            strip = strip.custom(strip.names = TRUE),
            xlab = "regression z-score",
            ylab = "binomial z-score",
            aspect = "iso")
 })

dev.off()


peakSummaryPTubeBlast <-
    diffPeakSummary(ptubes.ext, cblast.ext,
                    chrom.lens = mouse.chromlens,
                    lower = 15, islands = FALSE, merge = 20L)

peakSummaryPTubeBlast <-
    within(peakSummaryPTubeBlast,
       {
           promoter <-
               ifelse(computeOverlap(chromosome, start, end,
                                     gpromoters$chr, 
                                     gpromoters$start, 
                                     gpromoters$end),
                      "In promoter", "Not in promoter")
           CpG <-
               ifelse(computeOverlap(chromosome, start, end,
                                     CpG.mm9$chr, 
                                     CpG.mm9$start, 
                                     CpG.mm9$end),
                      "In CpG island", "Not in CpG island")
           comb.overlap <- overlap1 + overlap2
           ## p <- sum(overlap2) / sum(comb.overlap)
           p <- median(overlap2 / comb.overlap)
           z <- ((overlap2 / comb.overlap) - p) / sqrt(p * (1-p) / comb.overlap)
       })


peakSummaryTubeTube <-
    diffPeakSummary(ctubes.ext, ptubes.ext,
                    chrom.lens = mouse.chromlens,
                    lower = 15, islands = FALSE, merge = 20L)


peakSummaryTubeTube <-
    within(peakSummaryTubeTube,
       {
           promoter <-
               ifelse(computeOverlap(chromosome, start, end,
                                     gpromoters$chr, 
                                     gpromoters$start, 
                                     gpromoters$end),
                      "In promoter", "Not in promoter")
           CpG <-
               ifelse(computeOverlap(chromosome, start, end,
                                     CpG.mm9$chr, 
                                     CpG.mm9$start, 
                                     CpG.mm9$end),
                      "In CpG island", "Not in CpG island")
           comb.overlap <- overlap1 + overlap2
           p <- sum(overlap2) / sum(comb.overlap)
           ## p <- mean(overlap2 / comb.overlap)
           z <- ((overlap2 / comb.overlap) - p) / sqrt(p * (1-p) / comb.overlap)
       })


pdf("cpg-effect.pdf", width = 8, height = 11)

xyplot(sqrt(maxs2) ~ sqrt(maxs1) | CpG + promoter,
       data = peakSummaryTubeBlast,
       xlab = "Max depth in myotubes",
       ylab = "Max depth in myoblasts",
       panel = panel.smoothScatter,
       xscale.components = sqrt.scale.comps("x"),
       yscale.components = sqrt.scale.comps("y"),
       main = "Combined peaks (depth >= 15, merged if gap <= 20)",
       aspect = "iso")

## xyplot(asinh(overlap2) ~ asinh(overlap1) | CpG + promoter,
##        data = peakSummaryTubeBlast,
##        xlab = "asinh(number of overlapping reads) in Myotube",
##        ylab = "asinh(number of overlapping reads) in Fibroblast+MyoD",
##        panel = panel.smoothScatter,
##        main = "Fibroblast+MyoD and Myotube combined peaks\n(depth >= 15, merged if gap <= 20)",
##        aspect = "iso")

xyplot(z ~ log2(comb.overlap) | CpG + promoter, peakSummaryTubeBlast,
       main = "Myotube and myoblast combined peaks\n(depth >= 15, merged if gap <= 20)",
       panel = function(...) {
           panel.smoothScatter(...)
           panel.grid(h = -1, v = -1)
       })



xyplot(sqrt(maxs2) ~ sqrt(maxs1) | CpG + promoter,
       data = peakSummaryTubeTube,
       xlab = "Max depth in combined myotubes",
       ylab = "Max depth in paired-end myotubes",
       panel = panel.smoothScatter,
       xscale.components = sqrt.scale.comps("x"),
       yscale.components = sqrt.scale.comps("y"),
       main = "Combined peaks (depth >= 15, merged if gap <= 20)",
       aspect = "iso")

## xyplot(asinh(sums2) ~ asinh(sums1) | CpG + promoter,
##        data = peakSummaryTubeTube,
##        xlab = "asinh(peak coverage) in combined myotubes",
##        ylab = "asinh(peak coverage) in paired-end myotubes",
##        panel = panel.smoothScatter,
##        main = "Combined and paired-end myotube combined peaks\n(depth >= 15, merged if gap <= 20)",
##        aspect = "iso")

xyplot(z ~ log2(comb.overlap) | CpG + promoter, peakSummaryTubeTube,
       main = "Combined and paired-end myotube combined peaks\n(depth >= 15, merged if gap <= 20)",
       panel = function(...) {
           panel.smoothScatter(...)
           panel.grid(h = -1, v = -1)
       })




xyplot(sqrt(maxs2) ~ sqrt(maxs1) | CpG + promoter,
       data = peakSummaryPTubeBlast,
       xlab = "Max depth in paired-end Myotube",
       ylab = "Max depth in myoblasts",
       panel = panel.smoothScatter,
       xscale.components = sqrt.scale.comps("x"),
       yscale.components = sqrt.scale.comps("y"),
       main = "Combined peaks(depth >= 15, merged if gap <= 20)",
       aspect = "iso")

## xyplot(asinh(sums2) ~ asinh(sums1) | CpG + promoter,
##        data = peakSummaryPTubeBlast,
##        xlab = "asinh(peak coverage) in paired-end Myotube",
##        ylab = "asinh(peak coverage) in Fibroblast+MyoD",
##        panel = panel.smoothScatter,
##        main = "Fibroblast+MyoD and paired-end Myotube combined peaks\n(depth >= 15, merged if gap <= 20)",
##        aspect = "iso")

xyplot(z ~ log2(comb.overlap) | CpG + promoter, peakSummaryPTubeBlast,
       main = "Myoblast and paired-end myotube combined peaks\n(depth >= 15, merged if gap <= 20)",
       panel = function(...) {
           panel.smoothScatter(...)
           panel.grid(h = -1, v = -1)
       })


dev.off()



with(peakSummaryTubeBlast,
     xyplot(z ~ start | factor(chromosome, levels = unique(chromosome)),
            subset = (CpG == "Not in CpG island"),
            layout = c(1, 19), strip = FALSE, strip.left = TRUE,
            main = "Myoblast and myotube combined peaks\n(depth >= 15, merged if gap <= 20)",
            panel = function(x, y, ...) {
                panel.grid(h = -1, v = -1)
                panel.lines(x, y, ...)
                print(acf(y, plot = FALSE)$acf[2])
            }))




with(peakSummaryTubeBlast,
     xyplot(z ~ start | factor(chromosome, levels = unique(chromosome)),
            subset = (CpG == "Not in CpG island"),
            subscripts = TRUE,
            main = "Myoblast and myotube combined peaks\n(depth >= 15, merged if gap <= 20)",
            prepanel = function(x, y) {
                rng <- range(y, finite = TRUE)
                list(xlim = rng, ylim = rng)
            },            
            panel = function(x, y, ...) {
                ord <- order(x)
                y <- y[ord]
                n <- length(y)
                panel.smoothScatter(y[-n], y[-1], ...)
                ## panel.xyplot(y[-n], y[-1], pch = ".", cex = 2)
                ## panel.grid(h = -1, v = -1)
                panel.lmline(y[-n], y[-1], col = "black")
            }))


pdf("chromatin-effect.pdf", width = 8, height = 11)

with(peakSummaryTubeBlast,
     xyplot(z ~ start | factor(chromosome, levels = unique(chromosome)),
            subset = (CpG != "Not in CpG island"),
            subscripts = TRUE,
            main = "Myoblast and myotube combined peaks\n(depth >= 15, merged if gap <= 20)",
            xlab = "Distance from start of last peak",
            ylab = expression(Z[i] %*% Z[i-1]),
            prepanel = function(x, y) {
                ord <- order(x)
                x <- x[ord]
                y <- y[ord]
                n <- length(y)
                xx <- diff(x)
                yy <- y[-n] * y[-1]
                keep <- xx > 200 & xx < 3000
                list(xlim = range(xx[keep]), ylim = range(yy[keep]))
            },            
            panel = function(x, y, ...) {
                ord <- order(x)
                x <- x[ord]
                y <- y[ord]
                n <- length(y)
                xx <- diff(x)
                yy <- y[-n] * y[-1]
                keep <- xx < 3000
                ## panel.smoothScatter(xx[keep], yy[keep], ...)
                panel.grid(h = -1, v = -1)
                panel.xyplot(xx[keep], yy[keep], pch = ".", cex = 2, ...)
            },
            ylim = c(-10, 10)))


with(peakSummaryTubeBlast,
     xyplot(z ~ start | factor(chromosome, levels = unique(chromosome)),
            subset = (CpG == "Not in CpG island"),
            subscripts = TRUE,
            main = "Myoblast and myotube combined peaks\n(depth >= 15, merged if gap <= 20)",
            sub = "Non-CpG Peaks with last peak start site within 3000 bases", 
            xlab = "Z-score for last peak",
            ylab = "Z-score for current peak",
            prepanel = function(x, y) {
                ord <- order(x)
                x <- x[ord]
                y <- y[ord]
                n <- length(y)
                keep <- diff(x) > 200 & diff(x) < 3000
                xx <- y[-n]
                yy <- y[-1]
                list(xlim = range(xx[keep]), ylim = range(yy[keep]))
            },            
            panel = function(x, y, ...) {
                ord <- order(x)
                x <- x[ord]
                y <- y[ord]
                n <- length(y)
                keep <- diff(x) > 200 & diff(x) < 3000
                xx <- y[-n]
                yy <- y[-1]
                ## panel.smoothScatter(xx[keep], yy[keep], ...)
                panel.grid(h = -1, v = -1)
                panel.xyplot(xx[keep], yy[keep], pch = ".", cex = 2, ...)
            }))


dev.off()







with(peakSummaryTubeBlast,
     tapply(asinh(sums1)-asinh(sums2), CpG, length))
with(peakSummaryTubeBlast,
     tapply(asinh(sums1)-asinh(sums2), promoter, length))


with(peakSummaryTubeBlast,
     tapply(asinh(sums1)-asinh(sums2), CpG, mad))


with(peakSummaryTubeBlast,
     tapply(asinh(sums1)-asinh(sums2), list(chromosome, promoter), length))
foo <-
    with(peakSummaryTubeBlast,
         tapply(asinh(sums1)-asinh(sums2), list(chromosome, promoter), mad))


with(peakSummaryTubeBlast, tapply(end-start, list(promoter), length))
