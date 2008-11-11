### Function for obtaining regions for genomic context analysis

## 'genes' should be data.frame from dump of UCSC's "knownGenes" table
## Eventually, that should be available from an annotation package
## Mark added the gene ends for this release, but exons apparently won't make it
## 'proximal' is the number of bases on either side of TSS and 3'-end for
##   the promoter and end region, respectively.
## 'distal' is the number of bases on either side for upstream/downstream,
##   i.e. enhancer/silencer regions.
genomic_regions <- function(genes, proximal = 500, distal = 10000) {

  genomic_regions_chrom <- function(sub_genes) {
    gc()

    ## direction of transcription depends on strand (but UCSC is only by +)
    starts <- ifelse(sub_genes$strand == "+", sub_genes$txStart,
                     sub_genes$txEnd)
    ends <- ifelse(sub_genes$strand == "+", sub_genes$txEnd, sub_genes$txStart)
    proximal <- ifelse(sub_genes$strand == "+", proximal, -proximal)
    distal <- ifelse(sub_genes$strand == "+", distal, -distal)
    offset <- ifelse(sub_genes$strand == "+", 1, -1)

    rangebind <- function(x, y) cbind(start = pmin(x,y), end = pmax(x,y))
    
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
               upstream = upstream, downstream = downstream, gene = instream)
  }
  
  do.call("rbind", as.list(by(genes, genes$chrom, genomic_regions_chrom)))
}

genomic_exons <- function(genes) {
  genomic_exons_chrom <- function(sub_genes) {
    ## [exon:start, exon:end]
    splitPos <- function(pos) {
      as.integer(unlist(strsplit(as.character(pos), ",")))
    }
    data.frame(chrom = sub_genes$chrom[1],
               gene = rep(sub_genes$name, sub_genes$exonCount),
               start = splitPos(sub_genes$exonStarts),
               end = splitPos(sub_genes$exonEnds))
  }
  do.call("rbind", as.list(by(genes, genes$chrom, genomic_exons_chrom)))
}



countHits <- function(subject, query)
{
    sum(!is.na(overlap(subject, query, multiple = FALSE)))
}

contextDistributionPeakSet <- function(peaks, gregions)
{
    query <- with(peaks, IRanges(start, end))
    irangeByType <-
        function(type = c("promoter", "threeprime",
                          "upstream", "downstream", "gene"))
        {
            type <- match.arg(type)
            istarts <- sprintf("%s.start", type)
            iends <- sprintf("%s.end", type)
            keep <- !duplicated(gregions[[istarts]]) ## what's the right thing to do here???
            IRanges(start = gregions[[istarts]][keep],
                    end = gregions[[iends]][keep])
        }
    subject <-
        list(promoter = irangeByType("promoter"),
             threeprime = irangeByType("threeprime"),
             upstream = irangeByType("upstream"),
             downstream = irangeByType("downstream"),
             gene = irangeByType("gene"))
    c(total = length(query),
      sapply(subject, countHits, query = query))
}


## FIXME: need more flexible interface to define subgroups such as
## 'up' and 'down'

contextDistributionByChr <- function(chr, peaks, gregions)
{
    gregions.sub <- subset(gregions, chrom == chr)
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
             chroms = levels(peaks$chromosome),
             ...)
{
    ans <-
        do.call(lattice::make.groups,
                sapply(chroms, function(chr) {
                    contextDistributionByChr(chr, peaks, gregions)
                }, simplify = FALSE))
    names(ans)[names(ans) == "which"] <- "chromosome"
    rownames(ans) <- NULL
    ans
}




