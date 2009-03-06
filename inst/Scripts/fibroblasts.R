
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

all.reads <- GenomeDataList(list(realMouse = combineLaneReads(solexa54[c("7", "8")]),
                                 ctubes = combineLaneReads(myodMyo[c("2", "4", "7")]),
                                 cfibromyod = combineLaneReads(myodFibro[c("2", "4", "7")])))

ereads <- gdApply(all.reads,
                  function(x, seqLen = 200) {
                      sort(extendReads(x, seqLen = seqLen))
                  })

realMousePeaks <-
    gdApply(ereads$realMouse,
            function(g, cutoff = 6) {
                print(length(g))
                IntervalTree(slice(coverage(g, 1L, max(end(g)) + 100L), lower = cutoff))
            })


## promoter information

data(geneMouse)
gregions <- genomic_regions(genes = geneMouse, proximal = 1000)
gregions <- subset(gregions, chrom %in% paste("chr", 1:19, sep = ""))
gregions$chrom <- gregions$chrom[drop = TRUE]

gpromoters <- gregions[c("chrom", "promoter.start", "promoter.end")]
names(gpromoters) <- c("chr", "start", "end")
gpromoters.split <- split(gpromoters, gpromoters$chr)

## accumulate per-peak information

realMousePeakSummary <-
    sapply(names(realMousePeaks),
           function(chr) {
                print(chr)
                peaks <- realMousePeaks[[chr]]
                in.promoter <- !is.na(overlap(with(gpromoters.split[[chr]], IRanges(start, end)),
                                              peaks, multiple = FALSE))
                data.frame(start = start(peaks),
                           end = end(peaks),
                           promoter = in.promoter, 
                           realMouse = as.numeric(as.table(t(overlap(peaks, ereads$realMouse[[chr]], multiple = TRUE)))),
                           ctubes = as.numeric(as.table(t(overlap(peaks, ereads$ctubes[[chr]], multiple = TRUE)))),
                           cfibromyod = as.numeric(as.table(t(overlap(peaks, ereads$cfibromyod[[chr]], multiple = TRUE)))))
            }, simplify = FALSE)


realMousePeakSummary.df <- do.call(make.groups, realMousePeakSummary)
rownames(realMousePeakSummary.df) <- NULL



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




## ## proportion of Real Mouse peaks with >= 5 C2C12 Myotube extended reads, in 

## rbind("Real Mouse >= 6"                 = prop.table(xtabs(~(ctubes >= 5), realMousePeakSummary.df, subset = realMouse >= 6))[2],
##       "Real Mouse >= 10"                = prop.table(xtabs(~(ctubes >= 5), realMousePeakSummary.df, subset = realMouse >= 10))[2],
##       "Real Mouse >= 6, promoter only"  = prop.table(xtabs(~(ctubes >= 5), realMousePeakSummary.df, subset = promoter & realMouse >= 6))[2],
##       "Real Mouse >= 10, promoter only" = prop.table(xtabs(~(ctubes >= 5), realMousePeakSummary.df, subset = promoter & realMouse >= 10))[2])


## ## proportion of Real Mouse peaks with >= 5 Fibroblast+MyoD extended reads, in

## rbind("Real Mouse >= 6"                 = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = realMouse >= 6))[2],
##       "Real Mouse >= 10"                = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = realMouse >= 10))[2],
##       "Real Mouse >= 6, promoter only"  = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = promoter & realMouse >= 6))[2],
##       "Real Mouse >= 10, promoter only" = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = promoter & realMouse >= 10))[2])

## ## proportion of Real Mouse peaks with >= 5 Fibroblast+MyoD extended reads, in

## rbind("Myotubes >= 6"                 = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = ctubes >= 6))[2],
##       "Myotubes >= 10"                = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = ctubes >= 10))[2],
##       "Myotubes >= 6, promoter only"  = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = promoter & ctubes >= 6))[2],
##       "Myotubes >= 10, promoter only" = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = promoter & ctubes >= 10))[2])



## ## proportion of Real Mouse peaks with >= 5 C2C12 Myotube extended reads, in

## Real Mouse >= 6                 0.4416646
## Real Mouse >= 10                0.6982119
## Real Mouse >= 6, promoter only  0.7014439
## Real Mouse >= 10, promoter only 0.8786127


## ## proportion of Real Mouse peaks with >= 5 Fibroblast+MyoD extended reads, in

## Real Mouse >= 6                 0.1504344
## Real Mouse >= 10                0.2735935
## Real Mouse >= 6, promoter only  0.4099446
## Real Mouse >= 10, promoter only 0.5662331


## ## proportion of Real Mouse peaks with >= 5 Fibroblast+MyoD extended reads, in

## Myotubes >= 6                 0.3537711
## Myotubes >= 10                0.4740258
## Myotubes >= 6, promoter only  0.5945051
## Myotubes >= 10, promoter only 0.7028571




## ## number of Real Mouse peaks with >= 5 C2C12 Myotube extended reads, in 

## rbind("Real Mouse >= 6"                 = (xtabs(~(ctubes >= 5), realMousePeakSummary.df, subset = realMouse >= 6))[],
##       "Real Mouse >= 10"                = (xtabs(~(ctubes >= 5), realMousePeakSummary.df, subset = realMouse >= 10))[],
##       "Real Mouse >= 6, promoter only"  = (xtabs(~(ctubes >= 5), realMousePeakSummary.df,
##                                                  subset = promoter & realMouse >= 6))[],
##       "Real Mouse >= 10, promoter only" = (xtabs(~(ctubes >= 5), realMousePeakSummary.df,
##                                                  subset = promoter & realMouse >= 10))[])

                    
## ## number of Real Mouse peaks with >= 5 Fibroblast+MyoD extended reads, in

## rbind("Real Mouse >= 6"                 = (xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = realMouse >= 6))[],
##       "Real Mouse >= 10"                = (xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = realMouse >= 10))[],
##       "Real Mouse >= 6, promoter only"  = (xtabs(~(cfibromyod >= 5), realMousePeakSummary.df,
##                                                  subset = promoter & realMouse >= 6))[],
##       "Real Mouse >= 10, promoter only" = (xtabs(~(cfibromyod >= 5), realMousePeakSummary.df,
##                                                  subset = promoter & realMouse >= 10))[])


#####

## proportion of Real Mouse peaks with >= 5 C2C12 Myotube extended reads, in 

computeRates <- function(data, cutoff = 5, real.cutoff = 6)
{
    props <- c(prop.table(xtabs(~(ctubes >= cutoff), data, subset = realMouse >= real.cutoff))[2],
               prop.table(xtabs(~(ctubes >= cutoff), data, subset = promoter & realMouse >= real.cutoff))[2],
               prop.table(xtabs(~(cfibromyod >= cutoff), data, subset = realMouse >= real.cutoff))[2],
               prop.table(xtabs(~(cfibromyod >= cutoff), data, subset = promoter & realMouse >= real.cutoff))[2])
    counts <- c(xtabs(~(ctubes >= cutoff), data, subset = realMouse >= real.cutoff)[2],
                xtabs(~(ctubes >= cutoff), data, subset = promoter & realMouse >= real.cutoff)[2],
                xtabs(~(cfibromyod >= cutoff), data, subset = realMouse >= real.cutoff)[2],
                xtabs(~(cfibromyod >= cutoff), data, subset = promoter & realMouse >= real.cutoff)[2])
    data.frame(cutoff = cutoff, real.cutoff = real.cutoff,
               promoter = factor(c("All", "Promoter", "All", "Promoter")),
               sample = factor(c("Myotubes", "Myotubes", "Fibroblast+MyoD", "Fibroblast+MyoD")),
               proportion = props, counts = counts)
}

####

props.6 <-
    do.call(rbind,
            lapply(1:10,
                   function(i) computeRates(realMousePeakSummary.df,
                                            cutoff = i,
                                            real.cutoff = 6)))
props.10 <-
    do.call(rbind,
            lapply(1:10,
                   function(i) computeRates(realMousePeakSummary.df,
                                            cutoff = i,
                                            real.cutoff = 10)))


xyplot(proportion ~ cutoff | promoter, props.6, type = "o", ylim = c(0, 1),
       groups = sample, auto.key = list(lines = TRUE, points = FALSE),
       ylab = "Proportion of peaks with >= cutoff reads",
       main = "Real mouse peaks, depth >= 6")

xyplot(proportion ~ cutoff | promoter, props.10, type = "o", ylim = c(0, 1),
       groups = sample, auto.key = list(lines = TRUE, points = FALSE),
       ylab = "Proportion of peaks with depth >= cutoff reads",
       main = "Real mouse peaks, >= 10")

xyplot(counts ~ cutoff | promoter, props.6, type = "o",
       scales = list(relation = "free", rot = 0),
       groups = sample, auto.key = list(lines = TRUE, points = FALSE),
       ylab = "Number of peaks with depth >= cutoff reads",
       main = "Real mouse peaks, >= 6")

xyplot(counts ~ cutoff | promoter, props.10, type = "o", 
       scales = list(relation = "free", rot = 0),
       groups = sample, auto.key = list(lines = TRUE, points = FALSE),
       ylab = "Number of peaks with >= cutoff reads",
       main = "Real mouse peaks, depth >= 10")



##### similar calculations with fibroblast+MyoD peaks


fibromyodPeaks <-
    gdApply(ereads$cfibromyod,
            function(g, cutoff = 6) {
                print(length(g))
                IntervalTree(slice(coverage(g, 1L, max(end(g)) + 100L), lower = cutoff))
            })


## accumulate per-peak information

fibromyodPeakSummary <-
    sapply(names(fibromyodPeaks),
           function(chr) {
                print(chr)
                peaks <- fibromyodPeaks[[chr]]
                in.promoter <- !is.na(overlap(with(gpromoters.split[[chr]], IRanges(start, end)),
                                              peaks, multiple = FALSE))
                data.frame(start = start(peaks),
                           end = end(peaks),
                           promoter = in.promoter, 
                           realMouse = as.numeric(as.table(t(overlap(peaks, ereads$realMouse[[chr]], multiple = TRUE)))),
                           ctubes = as.numeric(as.table(t(overlap(peaks, ereads$ctubes[[chr]], multiple = TRUE)))),
                           cfibromyod = as.numeric(as.table(t(overlap(peaks, ereads$cfibromyod[[chr]], multiple = TRUE)))))
            }, simplify = FALSE)

fibromyodPeakSummary.df <- do.call(make.groups, fibromyodPeakSummary)
rownames(fibromyodPeakSummary.df) <- NULL


splom(~data.frame(asinh(realMouse), asinh(ctubes), asinh(cfibromyod)),
      data = fibromyodPeakSummary.df,
      varnames = c("realMouse", "myotubes", "fibroblast+MyoD"),
      main = "asinh(number of reads overlapping fibroblast+MyoD peaks with depth >= 6)",
      panel = panel.smoothScatter)

splom(~data.frame(asinh(realMouse), asinh(ctubes), asinh(cfibromyod)),
      data = fibromyodPeakSummary.df,
      subset = promoter,
      varnames = c("realMouse", "myotubes", "fibroblast+MyoD"),
      main = "asinh(number of reads overlapping fibroblast+MyoD peaks with depth >= 6, promoters only)",
      panel = panel.smoothScatter)



computeRates <- function(data, cutoff = 5, fibro.cutoff = 6)
{
    props <- c(prop.table(xtabs(~(ctubes >= cutoff), data, subset = cfibromyod >= fibro.cutoff))[2],
               prop.table(xtabs(~(ctubes >= cutoff), data, subset = promoter & cfibromyod >= fibro.cutoff))[2],
               prop.table(xtabs(~(realMouse >= cutoff), data, subset = cfibromyod >= fibro.cutoff))[2],
               prop.table(xtabs(~(realMouse >= cutoff), data, subset = promoter & cfibromyod >= fibro.cutoff))[2])
    counts <- c(xtabs(~(ctubes >= cutoff), data, subset = cfibromyod >= fibro.cutoff)[2],
                xtabs(~(ctubes >= cutoff), data, subset = promoter & cfibromyod >= fibro.cutoff)[2],
                xtabs(~(realMouse >= cutoff), data, subset = cfibromyod >= fibro.cutoff)[2],
                xtabs(~(realMouse >= cutoff), data, subset = promoter & cfibromyod >= fibro.cutoff)[2])
    data.frame(cutoff = cutoff, fibro.cutoff = fibro.cutoff,
               promoter = factor(c("All", "Promoter", "All", "Promoter")),
               sample = factor(c("Myotubes", "Myotubes", "Real Mouse", "Real Mouse")),
               proportion = props, counts = counts)
}

####

props.6 <-
    do.call(rbind,
            lapply(1:10,
                   function(i) computeRates(fibromyodPeakSummary.df,
                                            cutoff = i,
                                            fibro.cutoff = 6)))
props.10 <-
    do.call(rbind,
            lapply(1:10,
                   function(i) computeRates(fibromyodPeakSummary.df,
                                            cutoff = i,
                                            fibro.cutoff = 10)))


xyplot(proportion ~ cutoff | promoter, props.6, type = "o", ylim = c(0, 1),
       groups = sample, auto.key = list(lines = TRUE, points = FALSE),
       ylab = "Proportion of peaks with >= cutoff reads",
       main = "Fibroblast+MyoD peaks, depth >= 6")

xyplot(proportion ~ cutoff | promoter, props.10, type = "o", ylim = c(0, 1),
       groups = sample, auto.key = list(lines = TRUE, points = FALSE),
       ylab = "Proportion of peaks with depth >= cutoff reads",
       main = "Fibroblast+MyoD peaks, >= 10")

xyplot(counts ~ cutoff | promoter, props.6, type = "o",
       scales = list(relation = "free", rot = 0),
       groups = sample, auto.key = list(lines = TRUE, points = FALSE),
       ylab = "Number of peaks with depth >= cutoff reads",
       main = "Fibroblast+MyoD peaks, >= 6")

xyplot(counts ~ cutoff | promoter, props.10, type = "o", 
       scales = list(relation = "free", rot = 0),
       groups = sample, auto.key = list(lines = TRUE, points = FALSE),
       ylab = "Number of peaks with >= cutoff reads",
       main = "Fibroblast+MyoD peaks, depth >= 10")

dev.off()



