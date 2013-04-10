
## regression: ``DE peaks''

library(chipseq)
library(lattice)
library(latticeExtra)
library(geneplotter)

library(chipseqData)


data(myodMyo)
data(myodFibro)

library(BSgenome.Mmusculus.UCSC.mm9)
mouse.chromlens <- seqlengths(Mmusculus)


computeOverlap <- function(chr, start, end, tchr, tstart, tend)
{
    tsplit <- split(IRanges(tstart, tend), tchr)
    ans <- logical(length(chr))
    for (chrom in names(tsplit))
    {
        id <- which(chr == chrom)
        ans[id] <- IRanges(start[id], end[id]) %over% tsplit[[chrom]]
    }
    ans
}    




data(geneMouse)
gregions <- genomic_regions(genes = geneMouse, proximal = 1000)
gregions <- subset(gregions, chrom %in% paste("chr", 1:19, sep = ""))
gregions$chrom <- gregions$chrom[drop = TRUE]
gpromoters <- gregions[c("chrom", "promoter.start", "promoter.end")]
names(gpromoters) <- c("chr", "start", "end")

data(CpG.mm9)






combinedMyo <- 
    GenomeDataList(list(cblasts = combineLaneReads(myodMyo[c("1","3","6")]),
                        ctubes = combineLaneReads(myodMyo[c("2","4","7")])))

combinedFibro <- 
    GenomeDataList(list(cfibro = combineLaneReads(myodFibro[c("1","3","6")]),
                        cfibromyod = combineLaneReads(myodFibro[c("2","4","7")])))

### DE islands/peaks

extRangesMyo <- gdApply(combinedMyo, extendReads, seqLen = 200)
extRangesFibro <- gdApply(combinedFibro, extendReads, seqLen = 200)

peakSummaryTubeBlast <-
    diffPeakSummary(extRangesMyo$ctubes, extRangesMyo$cblasts,
                    chrom.lens = mouse.chromlens,
                    lower = 20, islands = FALSE)

peakSummaryTubeFibro <-
    diffPeakSummary(extRangesMyo$ctubes, extRangesFibro$cfibromyod,
                    chrom.lens = mouse.chromlens,
                    lower = 20, islands = FALSE)

peakSummaryBlastFibro <-
    diffPeakSummary(extRangesMyo$cblasts, extRangesFibro$cfibromyod,
                    chrom.lens = mouse.chromlens, lower = 8)


peakSummaryTubeFibro <-
    within(peakSummaryTubeFibro,
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
       })

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



xyplot(sqrt(maxs2) ~ sqrt(maxs1) | CpG + promoter,
       data = peakSummaryTubeFibro,
       xlab = "Max depth in Myotube",
       ylab = "Max depth in Fibroblast+MyoD",
       panel = panel.smoothScatter,
       xscale.components = sqrt.scale.comps("x"),
       yscale.components = sqrt.scale.comps("y"),
       main = "Fibroblast+MyoD and Myotube combined peaks, depth >= 20",
       aspect = "iso")




pdf("diffPeaks.pdf")

xyplot(asinh(sums2) ~ asinh(sums1) | chromosome,
       data = peakSummaryTubeBlast, 
       subset = (chromosome %in% c("chr1", "chr2", "chr3", "chr4")),
       panel = function(x, y, ...) {
           ## panel.xyplot(x, y, ...)
           panel.smoothScatter(x, y, ...)
           panel.abline(median(y - x), 1)
       },
       main = "Myoblast vs Myotube",
       type = c("p", "g"), alpha = 0.5, aspect = "iso")

xyplot(asinh(sums2) ~ asinh(sums1) | chromosome,
       data = peakSummaryTubeFibro, 
       subset = (chromosome %in% c("chr1", "chr2", "chr3", "chr4")),
       panel = function(x, y, ...) {
           ## panel.xyplot(x, y, ...)
           panel.smoothScatter(x, y, ...)
           panel.abline(median(y - x), 1)
       },
       main = "Fibroblast+MyoD vs Myotube",
       type = c("p", "g"), alpha = 0.5, aspect = "iso")

xyplot(asinh(sums2) ~ asinh(sums1) | chromosome,
       data = peakSummaryBlastFibro, 
       subset = (chromosome %in% c("chr1", "chr2", "chr3", "chr4")),
       panel = function(x, y, ...) {
           ## panel.xyplot(x, y, ...)
           panel.smoothScatter(x, y, ...)
           panel.abline(median(y - x), 1)
       },
       main = "Fibroblast+MyoD vs Myoblast",
       type = c("p", "g"), alpha = 0.5, aspect = "iso")

dev.off()



stop()



peakSummaryTubeFibro <- 
    within(peakSummaryTubeFibro,
       {
           diffs <- asinh(sums2) - asinh(sums1)
           resids <- (diffs - median(diffs)) / mad(diffs)
           up <- resids > 2
           down <- resids < -2
       })

xyplot(asinh(sums2) ~ asinh(sums1), # | (up + 2 * down),
       data = peakSummaryTubeFibro, 
       panel = panel.smoothScatter,
       ## pch = ".", #alpha = 0.5,
       main = "Fibroblast+MyoD vs Myotube",
       aspect = "iso")



peakSummary <- 
    within(peakSummary,
       {
           diffs <- asinh(sums2) - asinh(sums1)
           resids <- (diffs - median(diffs)) / mad(diffs)
           up <- resids > 2
           down <- resids < -2
       })



head(peakSummary)

data(geneMouse)
gregions <- genomic_regions(genes = geneMouse, proximal = 500)
gregions$gene <- as.character(gregions$gene)
str(gregions)

ans <- contextDistribution(peakSummary, gregions)
head(ans)

sumtab <- 
    as.data.frame(rbind(total = xtabs(total ~ type, ans),
                        promoter = xtabs(promoter ~ type, ans),
                        threeprime = xtabs(threeprime ~ type, ans),
                        upstream = xtabs(upstream ~ type, ans),
                        downstream = xtabs(downstream ~ type, ans),
                        gene = xtabs(gene ~ type, ans)))

cbind(sumtab, ratio = round(sumtab[, "down"] /  sumtab[, "up"], 3))

