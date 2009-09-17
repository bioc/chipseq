

pdf("antibody_comparison.pdf", width = 12, height = 10)


## are the three antibodies the same?  

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

all.reads <- GenomeDataList(list(realMouse.6975 = solexa54[["7"]],
                                 realMouse.6196 = solexa54[["8"]],
                                 myotubes.7311 = myodMyo[["2"]],
                                 myotubes.6975 = myodMyo[["4"]],
                                 myotubes.6196 = myodMyo[["7"]],
                                 combined = combineLaneReads( c(solexa54[c("7", "8")], myodMyo[c("2", "4", "7")]) )))


ereads <- gdApply(all.reads,
                  function(x, seqLen = 200) {
                      sort(extendReads(x, seqLen = seqLen))
                  })

## > cbind(unlist(lapply(gdApply(ereads, length), function(x) sum(unlist(x)))))
##                    [,1]
## realMouse.6975  4344784
## realMouse.6196  4183563
## myotubes.7311   2364119
## myotubes.6975   1534874
## myotubes.6196   1901992
## combined       14329332



combinedPeaks <-
    gdApply(ereads$combined,
            function(g, cutoff = 6) {
                print(length(g))
                IntervalTree(slice(coverage(g, width = max(end(g)) + 100L), lower = cutoff))
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



PeakSummary <-
    sapply(names(combinedPeaks),
           function(chr) {
                print(chr)
                peaks <- combinedPeaks[[chr]]
                in.promoter <- !is.na(findOverlaps(peaks,
                    with(gpromoters.split[[chr]], IRanges(start, end)),
                    multiple = FALSE))
                countOverlapping <- function(x)
                {
                    as.numeric(as.table(t(findOverlaps(ereads[[x]][[chr]], peaks, multiple = TRUE))))
                }
                data.frame(start = start(peaks),
                           end = end(peaks),
                           promoter = in.promoter, 
                           realMouse.6975 = countOverlapping("realMouse.6975"),
                           realMouse.6196 = countOverlapping("realMouse.6196"),
                           myotubes.7311 = countOverlapping("myotubes.7311"),
                           myotubes.6975 = countOverlapping("myotubes.6975"),
                           myotubes.6196 = countOverlapping("myotubes.6196"))
            }, simplify = FALSE)

PeakSummary.df <- do.call(make.groups, PeakSummary)
rownames(PeakSummary.df) <- NULL

####

splom(~data.frame(asinh(realMouse.6975), asinh(realMouse.6196),
                  asinh(myotubes.7311), asinh(myotubes.6975), asinh(myotubes.6196)),
      data = PeakSummary.df, ##pch = ".", cex = 2, alpha = 0.3)
      varnames = c("realMouse.6975", "realMouse.6196",
                  "myotubes.7311", "myotubes.6975", "myotubes.6196"),
      main = "asinh(number of reads overlapping combined peaks with depth >= 6)",
      panel = panel.smoothScatter)

## splom(~data.frame(asinh(realMouse.6975), asinh(realMouse.6196),
##                   asinh(myotubes.7311), asinh(myotubes.6975), asinh(myotubes.6196)),
##       data = PeakSummary.df, ##pch = ".", cex = 2, alpha = 0.3)
##       subset = promoter,
##       varnames = c("realMouse.6975", "realMouse.6196",
##                   "myotubes.7311", "myotubes.6975", "myotubes.6196"),
##       main = "asinh(number of reads overlapping combined peaks with depth >= 6)",
##       panel = panel.smoothScatter)



#### Rates of overlap


computeRates <-
    function(data,
             ref = c("realMouse.6975", "realMouse.6196",
                     "myotubes.7311", "myotubes.6975",
                     "myotubes.6196"),
             cutoff = 5, ref.cutoff = 6)
{
    ref <- match.arg(ref)
    rest <- setdiff(c("realMouse.6975", "realMouse.6196", "myotubes.7311", "myotubes.6975", "myotubes.6196"),
                    ref)
    dsub <- data[data[[ref]] >= ref.cutoff, ]
    dsub.promoter <- subset(dsub, promoter)
    props <- c(sapply(rest, function(x) prop.table(table(dsub[[x]] >= cutoff))[2]),
               sapply(rest, function(x) prop.table(table(dsub.promoter[[x]] >= cutoff))[2]))
    counts <- c(sapply(rest, function(x) table(dsub[[x]] >= cutoff)[2]),
                sapply(rest, function(x) table(dsub.promoter[[x]] >= cutoff)[2]))
    data.frame(cutoff = cutoff, ref = ref, ref.cutoff = ref.cutoff,
               promoter = rep(c("All", "Promoter"), each = length(rest)),
               sample = rep(rest, 2),
               proportion = props, counts = counts,
               stringsAsFactors = FALSE)
}


####

props.6 <-
    do.call(rbind,
            lapply(c("realMouse.6975", "realMouse.6196",
                     "myotubes.7311", "myotubes.6975",
                     "myotubes.6196"),
                   function(ref) 
               {
                   do.call(rbind,
                           lapply(1:10,
                                  function(i) computeRates(PeakSummary.df,
                                                           cutoff = i,
                                                           ref = ref,
                                                           ref.cutoff = 6)))
               }))

props.10 <-
    do.call(rbind,
            lapply(c("realMouse.6975", "realMouse.6196",
                     "myotubes.7311", "myotubes.6975",
                     "myotubes.6196"),
                   function(ref) 
               {
                   do.call(rbind,
                           lapply(1:10,
                                  function(i) computeRates(PeakSummary.df,
                                                           cutoff = i,
                                                           ref = ref,
                                                           ref.cutoff = 10)))
               }))

####



xyplot(proportion ~ cutoff | ref + promoter, props.6, type = c("l", "g"), 
       groups = sample, auto.key = list(lines = TRUE, points = FALSE, columns = 2),
       ylab = "Proportion of reference peaks with >= cutoff overlapping reads",
       main = "Combined peaks, depth >= 6.  Reference (orange strip) has >= 6 overlapping reads")

xyplot(proportion ~ cutoff | ref + promoter, props.10, type = c("l", "g"), 
       groups = sample, auto.key = list(lines = TRUE, points = FALSE, columns = 2),
       ylab = "Proportion of reference peaks with >= cutoff overlapping reads",
       main = "Combined peaks, depth >= 6.  Reference (orange strip) has >= 10 overlapping reads")



xyplot(counts ~ cutoff | ref + promoter, props.6, type = c("l", "g"),
       scales = list(relation = "free", rot = 0),
       groups = sample, auto.key = list(lines = TRUE, points = FALSE, columns = 2),
       ylab = "Number of reference peaks with >= cutoff overlapping reads",
       main = "Combined peaks, depth >= 6.  Reference (orange strip) has >= 6 overlapping reads")

xyplot(counts ~ cutoff | ref + promoter, props.10, type = c("l", "g"),
       scales = list(relation = "free", rot = 0),
       groups = sample, auto.key = list(lines = TRUE, points = FALSE, columns = 2),
       ylab = "Number of reference peaks with >= cutoff overlapping reads",
       main = "Combined peaks, depth >= 6.  Reference (orange strip) has >= 10 overlapping reads")


dev.off()

