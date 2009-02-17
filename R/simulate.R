
## quality scores needs to be a 3way array,
## you can create one using something like this
## 
# sp <- SolexaPath("/mnt/fred/solexa/ycao/080623_HWI-EAS88_0001")
# rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")
# abc <- alphabetByCycle(quality(rfq))
## this gives a three way array, rows areq


seqQScores <- function(inScores) {
    ##set up the quality scores
    ##we need to translate to integer from ascii,
    ##then average the integer codes and translate back
    q = rawToChar(as.raw(32:(32+93))) # these are valid quality encodings
    intQ = as(SFastqQuality(q), "numeric") + 33L # actually, 'integer' :(

    ##we start with a 3-way array (letter by quality by position)
    ##and drop to a two-way (letter by position) containing the
    ##average quality for that letter and position
    x1 = apply(abc, c(1,3), function(x)
               if(sum(x) > 0) round(sum(x*intQ)/sum(x)) else 0)

    ##drop all the letters with zero scores
    x1 = x1[rowSums(x1) != 0, ]
    ##fix up the N's in the first 12, as they will have scores of
    ##zero, since Solexa requires the first 12 to be unambiguous
    ## calls
    x1[5, 1:12] = 36

    ##and now, translate back from integer codes to ascii
    ss = strsplit(q, "")[[1]]
    x2 = ss[x1]
    dim(x2) = dim(x1)
    dimnames(x2) = dimnames(x1)
    x2
}

## assume size is a named vector with names correspond to chromosomes
## whose values are the number of reads to simulate from that chromosome,
## genome is the name of a genome, eg Mmusculus, assume BSgenome
## is loaded, 
## readLength is the length of the reads
## qualityScores is a matrix of 5 (ACGTN) by number of cycles
## replace is a logical for sampling with replacement
simulateReads <- function(size, genome, readLength, qualityScores,
                          replace=FALSE, verbose=FALSE) {
    if(ncol(qualityScores) != readLength)
        stop("ncol(qualityScores) != readLength")

    if(all(size < 1))
        size = size * seqlengths(genome)[names(size)]
    simValues =
      lapply(seq_len(length(size)), function(i) {
        if(verbose) cat("Iteration: ", i, " ", sep="")
        chr = names(size)[i]
        seq = unmasked(genome[[chr]])
        if(verbose) cat(".")
        chrlen = seqlengths(genome)[[chr]] - readLength
        starts = sample(chrlen, size[i], replace=replace)
        simReads =
          as.character(Views(seq, IRanges(start=starts, width=readLength)))
        if(verbose) cat(".")
        simChars = do.call(rbind, strsplit(simReads, ""))
        simQ = do.call(paste,
                       c(lapply(seq_len(readLength), function(i)
                                qualityScores[simChars[,i], i]),
                         sep=""))
        if(verbose) cat(".done\n")
        return(list(sread = simReads, quality = simQ,
                    id = paste(chr, starts, sep = ",")))
    })
    return(ShortReadQ(sread =
                      DNAStringSet(unlist(lapply(simValues, "[[", "sread"))),
                      quality =
                      SFastqQuality(unlist(lapply(simValues, "[[", "quality"))),
                      id = BStringSet(unlist(lapply(simValues, "[[", "id")))))
}
