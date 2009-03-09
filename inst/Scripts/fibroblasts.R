
pdf("fibroblast_comparison.pdf", width = 12, height = 10)

## are fibroblasts a good model system (for what)?

## real C2C12 mouse muscles is the gold standard
##  - C2C12 cell line myotube is a cell-line equivalent [?]
##  - fibroblasts (MEF) cell-line is a preferable model system.
##  - is it sufficiently good?

## Strategy:

##  - find real C2C12 peaks
##  - look at coverage under myotube (cell line) and fibroblasts (cell line)

library(lattice)
library(chipseq)
library(chipseqData)

data(solexa54)
data(myodMyo)
data(myodFibro)

## promoter information

data(geneMouse)
gregions <- genomic_regions(genes = geneMouse, proximal = 1000)
gregions <- subset(gregions, chrom %in% paste("chr", 1:19, sep = ""))
gregions$chrom <- gregions$chrom[drop = TRUE]

gpromoters <- gregions[c("chrom", "promoter.start", "promoter.end")]
names(gpromoters) <- c("chr", "start", "end")
gpromoters.split <- split(gpromoters, gpromoters$chr)

## samples being compared

combined <- combineLaneReads(c(solexa54[c("7", "8")],
                               myodMyo[c("2", "4", "7")],
                               myodFibro[c("2", "4", "7")]))

all.reads <- GenomeDataList(list(realMouse.6975 = solexa54[["7"]],
                                 realMouse.6196 = solexa54[["8"]],
                                 ctubes = combineLaneReads(myodMyo[c("2", "4", "7")]),
                                 cfibromyod = combineLaneReads(myodFibro[c("2", "4", "7")]),
                                 combined = combined))

ereads <- gdApply(all.reads,
                  function(x, seqLen = 200) {
                      sort(extendReads(x, seqLen = seqLen))
                  })

## lengths

do.call(cbind, lapply(ereads[1:2], function(x) unlist(lapply(x, length)) / 1e3))

countPeaks <- function(x, lower = c(10))
{
    cov <- coverage(x, 1, max(end(x)) + 500) 
    sapply(lower, function(i) length(slice(cov, lower = i)))
}

do.call(cbind, lapply(ereads[1:2], function(x) unlist(lapply(x, countPeaks, lower = 10))))




## 


summarizeData <- function(edata, peak.ref, peak.cutoff = 6, ref = peak.ref, ref.cutoff = peak.cutoff,
                          include = names(edata))
{
    peaks <-
        gdApply(edata[[peak.ref]],
                function(g, cutoff = peak.cutoff) {
                    print(length(g))
                    IntervalTree(slice(coverage(g, 1L, max(end(g)) + 100L), lower = cutoff))
                })
    ## accumulate per-peak information
    peakSummary <-
        sapply(names(peaks),
               function(chr) {
                   print(chr)
                   chrpeaks <- peaks[[chr]]
                   in.promoter <- !is.na(overlap(with(gpromoters.split[[chr]], IRanges(start, end)),
                                                 chrpeaks, multiple = FALSE))
                   countOverlapping <- function(x)
                   {
                       as.numeric(as.table(t(overlap(chrpeaks, edata[[x]][[chr]], multiple = TRUE))))
                   }
                   ans <- data.frame(start = start(chrpeaks),
                                     end = end(chrpeaks),
                                     promoter = in.promoter)
                   for (nm in names(edata))
                       ans[[nm]] <- countOverlapping(nm)
                   ans
               }, simplify = FALSE)
    peakSummary.df <- do.call(make.groups, peakSummary)
    rownames(peakSummary.df) <- NULL
    computeRates <- function(cutoff = 5)
    {
        mytab <- function(x) table(factor(x, levels = c(FALSE, TRUE)))
        dsub <- peakSummary.df ## [peakSummary.df[[ref]] >= ref.cutoff, ]
        dsub.promoter <- subset(dsub, promoter)
        props <- c(sapply(include, function(x) prop.table(mytab(dsub[[x]] >= cutoff))[2]),
                   sapply(include, function(x) prop.table(mytab(dsub.promoter[[x]] >= cutoff))[2]))
        counts <- c(sapply(include, function(x) mytab(dsub[[x]] >= cutoff)[2]),
                    sapply(include, function(x) mytab(dsub.promoter[[x]] >= cutoff)[2]))
        data.frame(cutoff = cutoff, ref = ref, ref.cutoff = ref.cutoff,
                   promoter = rep(c("All", "Promoter"), each = length(include)),
                   sample = rep(include, 2),
                   proportion = props, counts = counts,
                   stringsAsFactors = FALSE)
    }
    props <- do.call(rbind, lapply(1:15, computeRates))
    list(peakSummary = peakSummary.df, props = props)
}

fibroPeakSummary.6 <- summarizeData(ereads, peak.ref = "cfibromyod", peak.cutoff = 6,
                                    include = setdiff(names(ereads), c("combined")))

xyplot(proportion ~ cutoff | promoter, fibroPeakSummary.6$props, type = c("g", "o"), 
       groups = sample, auto.key = list(lines = TRUE, points = FALSE, columns = 2),
       ylab = "Proportion of peaks with >= cutoff reads",
       main = "fibroblasts+MyoD peaks, depth >= 6")


fibroPeakSummary.10 <- summarizeData(ereads, peak.ref = "cfibromyod", peak.cutoff = 10,
                                     include = setdiff(names(ereads), c("combined")))

ctubePeakSummary.10 <- summarizeData(ereads, peak.ref = "ctubes", peak.cutoff = 10,
                                     include = setdiff(names(ereads), c("combined")))



png("sample_comparison_%03d.png", width = 600, height = 400)

xyplot(proportion ~ cutoff | promoter, fibroPeakSummary.10$props, type = c("g", "o"), 
       groups = sample, auto.key = list(lines = TRUE, points = FALSE, columns = 2),
       ylab = "Proportion of peaks with number of overlapping reads >= cutoff ",
       main = "fibroblasts+MyoD peaks, depth >= 10")

## xyplot(counts ~ cutoff | promoter, fibroPeakSummary.10$props, type = c("g", "o"), 
##        groups = sample, auto.key = list(lines = TRUE, points = FALSE, columns = 2),
##        ylab = "Proportion of peaks with number of overlapping reads >= cutoff ",
##        main = "fibroblasts+MyoD peaks, depth >= 10")


xyplot(proportion ~ cutoff | promoter, ctubePeakSummary.10$props, type = c("g", "o"), 
       groups = sample, auto.key = list(lines = TRUE, points = FALSE, columns = 2),
       ylab = "Proportion of peaks with number of overlapping reads >= cutoff ",
       main = "Myotube peaks, depth >= 10")

dev.off()






library(geneplotter)

splom(~data.frame(asinh(realMouse), asinh(ctubes), asinh(cfibromyod)),
      data = realMousePeakSummary.df,
      varnames = c("realMouse", "myotubes", "fibroblasts"),
      main = "asinh(number of reads overlapping real_mouse peaks with depth >= 6)",
      panel = panel.smoothScatter)

splom(~data.frame(asinh(realMouse), asinh(ctubes), asinh(cfibromyod)),
      data = realMousePeakSummary.df,
      subset = promoter,
      varnames = c("realMouse", "myotubes", "fibroblasts"),
      main = "asinh(number of reads overlapping real_mouse peaks with depth >= 6, promoters only)",
      panel = panel.smoothScatter)



if (FALSE)
{
    require(geneplotter)
    FUN <- sqrt
    FUN <- log1p
    FUN <- asinh
    xyplot(FUN(ctubes) + FUN(cfibromyod) ~ FUN(realMouse), data = realMousePeakSummary.df,
           outer = TRUE,
           panel = panel.smoothScatter)
    xyplot(FUN(ctubes) + FUN(cfibromyod) ~ FUN(realMouse), data = realMousePeakSummary.df,
           outer = TRUE, pch = ".", cex = 2)
}




## #####


## computeRates <-
##     function(data,
##              ref = c("realMouse", "ctubes", "cfibromyod"),
##              cutoff = 5, ref.cutoff = 6)
## {
##     ref <- match.arg(ref)
##     rest <- setdiff(c("realMouse", "ctubes", "cfibromyod"),
##                     ref)
##     dsub <- data[data[[ref]] >= ref.cutoff, ]
##     dsub.promoter <- subset(dsub, promoter)
##     props <- c(sapply(rest, function(x) prop.table(table(dsub[[x]] >= cutoff))[2]),
##                sapply(rest, function(x) prop.table(table(dsub.promoter[[x]] >= cutoff))[2]))
##     counts <- c(sapply(rest, function(x) table(dsub[[x]] >= cutoff)[2]),
##                 sapply(rest, function(x) table(dsub.promoter[[x]] >= cutoff)[2]))
##     data.frame(cutoff = cutoff, ref = ref, ref.cutoff = ref.cutoff,
##                promoter = rep(c("All", "Promoter"), each = length(rest)),
##                sample = rep(rest, 2),
##                proportion = props, counts = counts,
##                stringsAsFactors = FALSE)
## }




