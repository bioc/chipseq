### Function for obtaining regions for genomic context analysis

### FIXME: could push up nearestTss to GenomicFeatures
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
  DataFrame(in.cds = x %in% cdsBy(txdb),
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
