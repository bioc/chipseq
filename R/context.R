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
