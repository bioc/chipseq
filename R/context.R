### Function for obtaining regions for genomic context analysis

### These functions do not seem that generally useful. Incorporating
### overlap information (e.g. as a logical vector) into a RangedData
### is not that difficult. From there, an analysis might take many
### paths.

setGeneric("genomicContext", function(x, ...) standardGeneric("genomicContext"))

setGeneric("genomicContext<-",
           function(x, ...) standardGeneric("genomicContext<-"))


contextDistributionPeakSet <- function(peaks, gregions)
{
    query <- with(peaks, IRanges(start, end))
    subject <- values(gregions)[[1]]
    subject$gene <- ranges(gregions)[[1]]
    c(total = length(query),
      sapply(subject, function(x) sum(query %in% x)))
}


## FIXME: need more flexible interface to define subgroups such as
## 'up' and 'down'

contextDistributionByChr <- function(chr, peaks, gregions)
{
    gregions.sub <- gregions[chr]
    peaks.sub <- subset(peaks, chromosome == chr)
    all <- contextDistributionPeakSet(peaks.sub, gregions.sub)
    up <- 
        if ("up" %in% names(peaks))
            contextDistributionPeakSet(subset(peaks.sub, up), gregions.sub)
        else NULL
    down <- 
        if ("down" %in% names(peaks))
            contextDistributionPeakSet(subset(peaks.sub, down), gregions.sub)
        else NULL
    ans <- as.data.frame(rbind(all = all, up = up, down = down))
    ans <- cbind(type = factor(rownames(ans), levels = unique(rownames(ans))), 
                 ans)
    ans
}

contextDistribution <-
    function(peaks, gregions,
             chroms = unique(as.character(peaks[["chromosome"]])),
             ...)
{
    .Deprecated("the GenomicFeatures package directly")
    if (!is.data.frame(peaks))
      stop("'peaks' must be a data.frame")
    if (!all(c("chromosome", "start", "end") %in% names(peaks)))
      stop("'peaks' must have column names 'chromosome', 'start', and 'end'")
    if (!is(gregions, "RangedData"))
      stop("'gregions' must be a RangedData as returned by ",
           "'transcripts_deprecated' in the GenomicFeatures package")
    if (!is.character(chroms) || (length(chroms) == 0))
      stop("'chroms' must be a character vector")
    if (!all(chroms %in% peaks[["chromosome"]]))
      stop("'chroms' must be a character vector containing values",
           " from 'peaks[[\"chromosome\"]]'")
    if (!all(chroms %in% names(gregions)))
      stop("'chroms' must be a character vector containing values",
           " from 'names(gregions)'")
    ans <-
        do.call(lattice::make.groups,
                sapply(chroms, function(chr) {
                    contextDistributionByChr(chr, peaks, gregions)
                }, simplify = FALSE))
    names(ans)[names(ans) == "which"] <- "chromosome"
    rownames(ans) <- NULL
    ans
}
