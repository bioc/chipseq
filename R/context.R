### Function for obtaining regions for genomic context analysis

### These functions do not seem that generally useful. Incorporating
### overlap information (e.g. as a logical vector) into a RangedData
### is not that difficult. From there, an analysis might take many
### paths.

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

### FIXME: could push up to GenomicFeatures

.nearestTss <- function(x, txdb) {
  if (!is(txdb, "TranscriptDb"))
    stop("'txdb' must be a 'TranscriptDb' object")
  tss <- resize(transcripts(txdb), 1L)
  seqlevels(tss) <- seqlevels(x)
  n <- nearest(x, tss, ignore.strand = TRUE)
  DataFrame(tss.id = values(tss)$tx_id[n], tss.distance = ranges(tss)[n])
}

setGeneric("addNearestTss", function(x, ...) standardGeneric("addNearestTss"))

setMethod("addNearestTss", "GenomicRanges", function(x, txdb) {
  values(x) <- cbind(values(x), .nearestTss(x, txdb))
  x
})

.genomicContext <- function(x, txdb, proximal = 2000L, distal = 10000L) {
  tx <- transcripts(txdb)
  tx_ext <- tx + proximal
  stream_width <- distal - proximal
  tx_no_strand <- tx
  strand(tx_no_strand) <- "*"
  x_no_strand <- x
  strand(x_no_strand) <- "*"
  DataFrame(in.cds = x %in% cds_tx,
            in.3utr = x %in% threeUTRsByTranscript(txdb),
            in.5utr = x %in% fiveUTRsByTranscript(txdb),
            in.intron = x %in% intronsByTranscript(txdb),
            in.promoter = x %in% flank(tx, start=TRUE, width = proximal),
            in.3prime = x %in% flank(tx, start=FALSE, width = proximal),
            in.upstream = x %in% flank(tx_ext, start=TRUE,
              width=stream_width),
            in.downstream = x %in% flank(tx_ext, start=FALSE,
              width=stream_width),
            in.intergenic = x_no_strand %in% gaps(tx_no_strand + distal))
}

setGeneric("addGenomicContext",
           function(x, ...) standardGeneric("addGenomicContext"))

setMethod("addGenomicContext", "GenomicRanges",
          function(x, txdb, proximal = 2000L, distal = 10000L) {
            context <- .genomicContext(x, txdb, proximal = proximal,
                                       distal = distal)
            tss <- .nearestTss(x, txdb)
            values(x) <- cbind(values(x), context, tss)
            x
          })
