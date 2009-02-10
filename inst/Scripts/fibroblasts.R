

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






require(geneplotter)


FUN <- sqrt
FUN <- log1p
FUN <- asinh

xyplot(FUN(ctubes) + FUN(cfibromyod) ~ FUN(realMouse), data = realMousePeakSummary.df,
       outer = TRUE,
       panel = panel.smoothScatter)

xyplot(FUN(ctubes) + FUN(cfibromyod) ~ FUN(realMouse), data = realMousePeakSummary.df,
       outer = TRUE, pch = ".", cex = 2)




## proportion of Real Mouse peaks with >= 5 C2C12 Myotube extended reads, in 

rbind("Real Mouse >= 6"                 = prop.table(xtabs(~(ctubes >= 5), realMousePeakSummary.df, subset = realMouse >= 6))[2],
      "Real Mouse >= 10"                = prop.table(xtabs(~(ctubes >= 5), realMousePeakSummary.df, subset = realMouse >= 10))[2],
      "Real Mouse >= 6, promoter only"  = prop.table(xtabs(~(ctubes >= 5), realMousePeakSummary.df, subset = promoter & realMouse >= 6))[2],
      "Real Mouse >= 10, promoter only" = prop.table(xtabs(~(ctubes >= 5), realMousePeakSummary.df, subset = promoter & realMouse >= 10))[2])


## proportion of Real Mouse peaks with >= 5 Fibroblast+MyoD extended reads, in

rbind("Real Mouse >= 6"                 = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = realMouse >= 6))[2],
      "Real Mouse >= 10"                = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = realMouse >= 10))[2],
      "Real Mouse >= 6, promoter only"  = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = promoter & realMouse >= 6))[2],
      "Real Mouse >= 10, promoter only" = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = promoter & realMouse >= 10))[2])

## proportion of Real Mouse peaks with >= 5 Fibroblast+MyoD extended reads, in

rbind("Myotubes >= 6"                 = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = ctubes >= 6))[2],
      "Myotubes >= 10"                = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = ctubes >= 10))[2],
      "Myotubes >= 6, promoter only"  = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = promoter & ctubes >= 6))[2],
      "Myotubes >= 10, promoter only" = prop.table(xtabs(~(cfibromyod >= 5), realMousePeakSummary.df, subset = promoter & ctubes >= 10))[2])



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




                    
