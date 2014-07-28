### Function for obtaining regions for genomic context analysis

### FIXME: could push up nearestTss to GenomicFeatures
.nearestTss <- function(x, txdb) {
  if (!is(txdb, "TxDb"))
    stop("'txdb' must be a 'TxDb' object")
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
  DataFrame(in.cds = x %over% cdsBy(txdb),
            in.3utr = x %over% threeUTRsByTranscript(txdb),
            in.5utr = x %over% fiveUTRsByTranscript(txdb),
            in.intron = x %over% intronsByTranscript(txdb),
            in.promoter = x %over% flank(tx, start=TRUE, width = proximal),
            in.3prime = x %over% flank(tx, start=FALSE, width = proximal),
            in.upstream = x %over% flank(tx_ext, start=TRUE,
              width=stream_width),
            in.downstream = x %over% flank(tx_ext, start=FALSE,
              width=stream_width),
            in.intergenic = x_no_strand %over% gaps(tx_no_strand + distal))
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
