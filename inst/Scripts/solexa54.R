


## code related to the paper

## regression: ``DE peaks''

## Note: Zizhen suggested a permutation test, which we should try out
## (haven't yet)


library(chipseq)
library(lattice)
library(latticeExtra)
library(geneplotter)


load("solexa54.rda")

library(BSgenome.Mmusculus.UCSC.mm9)
mouse.chromlens <- seqlengths(Mmusculus)

## function to compute GC content

gcContent <- function(isochore = isochore.df,
                      chrom = "chr1")
{
    v <- with(isochore[[chrom]],
              Views(Mmusculus[[chrom]], start = start, end = end))
    rowSums(alphabetFrequency(v, baseOnly = TRUE, freq = TRUE)[, c("G", "C")])
}

combined54 <- 
    list(h0 = combineLaneReads(solexa54[c("1","2")]),
         h96 = combineLaneReads(solexa54[c("3","4")]))


## Step 3: Get ``peaks'' in combined data.

### DE islands/peaks

extRanges54 <- lapply(combined54, extendReads)


if (FALSE)  ## dump coverage
{ 

cov54 <- lapply(extRanges54, function(x) { 
    sapply(names(x), 
           function(chr) { 
               message(chr)
               coverage(x[[chr]], 1, mouse.chromlens[chr])
           }, simplify = FALSE)
    })

cov2df <- function(x, chr = "") {
    vals <- as.integer(x@values)
    lens <- as.integer(x@lengths)
    ends <- cumsum(lens)
    starts <- c(1L, head(ends, -1) + 1L)
    data.frame(chr = chr, start = starts, end = ends, coverage = vals, 
               stringsAsFactors = FALSE)
}

h0cov <- 
    do.call(rbind, 
            sapply(names(cov54$h0), 
                   function(chr) {
                       message(chr)
                       cov2df(cov54$h0[[chr]], chr = chr)
                   }, simplify = FALSE))

h96cov <- 
    do.call(rbind, 
            sapply(names(cov54$h96), 
                   function(chr) {
                       message(chr)
                       cov2df(cov54$h96[[chr]], chr = chr)
                   }, simplify = FALSE))

write.csv(h0cov, row.names = FALSE, quote = FALSE, file = "h0cov.csv")
write.csv(h96cov, row.names = FALSE, quote = FALSE, file = "h96cov.csv")

}



peakSummaryMethyl <-
    diffPeakSummary(extRanges54$h0, extRanges54$h96,
                    chrom.lens = mouse.chromlens, lower = 8, islands = FALSE)


with(peakSummaryMethyl, 
     xyplot(asinh(sums2) ~ asinh(sums1) | chromosome,
            subset = (chromosome %in%
                      c("chr1", "chr2", "chr3", "chr4",
                        "chr5", "chr6", "chr7", "chr8", "chr9")),
            panel = function(x, y, ...) {
                ## panel.xyplot(x, y, ...)
                panel.smoothScatter(x, y, ...)
                panel.abline(median(y - x), 1)
            },
            main = "Fibroblast+MyoD vs Myotube",
            type = c("p", "g"), alpha = 0.5, aspect = "iso"))


## Add columns to peak-summary data frame giving GC content and
## promoter information.  

## get contxt information so that we can identify promoter peaks

data(geneMouse)
gregions <- genomic_regions(genes = geneMouse, proximal = 500)
gregions <- subset(gregions, chrom %in% paste("chr", 1:19, sep = ""))
gregions$chrom <- gregions$chrom[drop = TRUE]

gpromoters <- gregions[c("chrom", "promoter.start", "promoter.end")]
names(gpromoters) <- c("chr", "start", "end")



### FIXME: isochore data had problems (no chr11 and above), lift-over
### needs updating

overlap.index <-
    function(object, query)
{
    overlap(with(object, IRanges(start, end)),
            with(query, IRanges(start, end)),
            multiple = FALSE)
}

updatePeakSummary <- function(x, promoters = gpromoters)
{
    x.split <- split(x, x$chr)
    promoter.split <- split(promoters, promoters$chr)
    for (chr in names(x.split))
    {
        message("Processing ", chr)
        ## peak GC content
        x.split[[chr]]$peakGC <- gcContent(x.split, chrom = chr)
        ## In promoter?
        oo <- overlap.index(promoter.split[[chr]], x.split[[chr]])
        x.split[[chr]]$promoter <- !is.na(oo)
    }
    ans <- do.call(rbind, x.split)
    rownames(ans) <- NULL
    ans$promoter <- factor(ans$promoter, levels = c(FALSE, TRUE),
                           labels = c("Outside promoter", "In promoter"))
    ans <-
        within(ans,
           {
               diffs <- asinh(sums2) - asinh(sums1)
           })
    ans
}


peakSummaryMethyl <- updatePeakSummary(peakSummaryMethyl)
peakSummaryMethyl <- na.omit(peakSummaryMethyl)

fmx.control <- lm(diffs ~ 1 + promoter + peakGC, peakSummaryMethyl)
anova(fmx.control)
## summary(fmx.control)

peakSummaryMethyl$fitted <- fitted(fmx.control)
peakSummaryMethyl$rstandard <- rstandard(fmx.control)
peakSummaryMethyl$rstudent <- rstudent(fmx.control)

qqmath(~rstandard | chromosome, peakSummaryMethyl, aspect = "iso",
       ## f.value = ppoints(1000),
       type = c("p", "g"), pch = ".")

qqmath(~rstandard, peakSummaryMethyl, aspect = "iso",
       f.value = ppoints(1000),
       type = c("p", "g"))

foo <- subset(peakSummaryMethyl, abs(rstandard) > 2.5) 

splom(foo[c("fitted", "rstandard", "rstudent")])


write.csv(foo[, 1:3], row.names = FALSE, file = "methyl_DE.csv")



