### Function for obtaining regions for genomic context analysis

## 'genes' should be data.frame from dump of UCSC's "knownGenes" table
## Eventually, that should be available from an annotation package
## Mark added the gene ends for this release, but exons apparently won't make it
## 'proximal' is the number of bases on either side of TSS and 3'-end for
##   the promoter and end region, respectively.
## 'distal' is the number of bases on either side for upstream/downstream,
##   i.e. enhancer/silencer regions.
genomic_regions <- function(genes, proximal = 500L, distal = 10000L) {

  .Deprecated("transcripts", "chipseq",
              "Deprecated, use 'transcripts' from the GenomicFeatures package")

  genomic_regions_chrom <- function(sub_genes) {
    gc()

    ## direction of transcription depends on strand (but UCSC is only by +)
    starts <- ifelse(sub_genes$strand == "+", sub_genes$txStart,
                     sub_genes$txEnd)
    ends <- ifelse(sub_genes$strand == "+", sub_genes$txEnd, sub_genes$txStart)
    proximal <- ifelse(sub_genes$strand == "+", proximal, -proximal)
    distal <- ifelse(sub_genes$strand == "+", distal, -distal)
    offset <- ifelse(sub_genes$strand == "+", 1L, -1L)

    rangebind <- function(x, y) cbind(start = pmin.int(x,y), end = pmax.int(x,y))

    ## [start-proximal,start+proximal]    
    promoter <- rangebind(as.integer(starts - proximal),
                          as.integer(starts + proximal))
    ## [end-proximal,end+proximal]
    threeprime <- rangebind(as.integer(ends - proximal),
                            as.integer(ends + proximal))
    ## [start-distal,start-1]
    upstream <- rangebind(as.integer(starts - distal),
                          as.integer(starts - proximal - offset))
    ## [end+1, end+distal]
    downstream <- rangebind(as.integer(ends + proximal + offset),
                            as.integer(ends + distal))
    ## [start, end]
    instream <- cbind(start = as.integer(sub_genes$txStart),
                      end = as.integer(sub_genes$txEnd))

    data.frame(chrom = sub_genes$chrom, gene = sub_genes$name,
               promoter = promoter, threeprime = threeprime,
               upstream = upstream, downstream = downstream,
               gene = instream)
  }

  do.call(rbind, as.list(by(genes, genes$chrom, genomic_regions_chrom)))
}

genomic_exons <- function(genes) {

  .Deprecated("transcripts", "chipseq",
              "Deprecated, use 'exons' from the GenomicFeatures package")

  splitPos <- function(pos) {
    as.integer(unlist(strsplit(as.character(pos), ",")))
  }
  genomic_exons_chrom <- function(sub_genes) {
    ## [exon:start, exon:end]
    data.frame(chrom = sub_genes$chrom[1],
               gene = rep(sub_genes$name, sub_genes$exonCount),
               start = splitPos(sub_genes$exonStarts),
               end = splitPos(sub_genes$exonEnds))
  }
  do.call(rbind, as.list(by(genes, genes$chrom, genomic_exons_chrom)))
}


## RG thinks that there is always one more exon than intron
## and that the first intron comes after the first exon
##
genomic_introns <- function(genes) {

  .Deprecated("transcripts", "chipseq",
              "Deprecated, use 'introns' from the GenomicFeatures package")

  ##a couple of helper functions - defined here so we
  ##don't instantiate every time throught he loop
  splitEnd <- function(pos) {
     sps = strsplit(as.character(pos), ",")
     lp = lapply(sps, function(x) x[-length(x)])
     as.integer(unlist(lp))
  }
  splitStart <- function(pos) {
      sps = strsplit(as.character(pos), ",")
      lp = lapply(sps, function(x) x[-1])
      as.integer(unlist(lp))
  }
  genomic_introns_chrom <- function(sub_genes) {
    ## introns are defined by [exon:end, exon:start]
    sub_genes = sub_genes[sub_genes$exonCount > 1,]
    starts = splitEnd(sub_genes$exonEnds)
    ends = splitStart(sub_genes$exonStarts)
    chrom <- 
        if (nrow(sub_genes) > 0) 
            sub_genes$chrom[1]
        else
            character()
    data.frame(chrom = chrom,
               gene = rep(sub_genes$name, sub_genes$exonCount-1),
               start = starts,
               end = ends)
  }
  do.call(rbind, as.list(by(genes, genes$chrom, genomic_introns_chrom)))
}



contextDistributionPeakSet <- function(peaks, gregions)
{
    query <- with(peaks, IRanges(start, end))
    subject <- values(gregions)[[1]]
    subject$gene <- ranges(gregions)[[1]]
    c(total = length(query),
      sapply(subject, function(x) countOverlaps(query, x)))
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
    if (!is.data.frame(peaks))
      stop("'peaks' must be a data.frame")
    if (!all(c("chromosome", "start", "end") %in% names(peaks)))
      stop("'peaks' must have column names 'chromosome', 'start', and 'end'")
    if (!is(gregions, "RangedData"))
      stop("'gregions' must be a RangedData as returned by 'transcripts'",
           " in the GenomicFeatures package")
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
