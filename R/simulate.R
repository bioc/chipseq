
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

##assume num is a named vector, names correspond to chromosomes
## values are the number of reads to simulate from that chromosome
## genome is the name of a genome, eg Mmusculus, assume BSgenome
## is loaded, 
## readLen is the length of the reads
## qualityScores is a matrix of 5 (ACGTN) by number of cycles
simulateReads <- function(num, genome, readLen, qualityScores, verbose=FALSE) {
    ##now we sample 
    if(ncol(qualityScores) != readLen )
        stop("need qualityScores for each position")

    ans = vector("list", length=length(num))
    prop = all(num < 1)
    names(ans) = names(num)
    for(i in 1:length(num)){
        chr = names(num)[i]
        seq = unmasked(genome[[chr]])
        nc = nchar(seq) - readLen ##don't run over the end
        ##set either proportions or numbers
        if(prop) myreads = sample(1:nc, num[i]*nc) else
            myreads = sample(1:nc, num[i])	
        myV = Views(seq, start=myreads, end = myreads+(readLen - 1))
        SimChars = as.character(myV)
        Char2Row = sapply(SimChars,
                          function(x) chartr("ACGTN", "12345", x))
        SimQ = sapply(Char2Row,
                      function(x) {
		                  r = as.integer(unlist(strsplit(x, "")))
                          return(qualityScores[matrix(c(r, 1:35), nc=2)])
                      })
        ans[[i]] = list(chromosome = chr,
                        start = myreads,
                        seqs = SimChars,
                        quality = apply(SimQ, 2, paste, 
                                        sep="", collapse=""))
        if(verbose) print(paste("iteration", i))
    }
    return(ans)
}

printSim <- function(x, filename="SimOut.txt") {
    #we need to turn the data into a fastQ format

    head = "@RGsim"
    head2 = "+RGsim"

    for(seqs in x) {
        numChars = length(seqs$seqs)
        ans = character(length=4*numChars)
        ans[seq(1, by=4, length.out=numChars)] =
          paste(head, ":", seqs$chr, ",", seqs$start, sep="")
        ans[seq(2, by=4, length.out=numChars)] = seqs$seqs
        ans[seq(3, by=4, length.out=numChars)] =
          paste(head2, ":", seqs$chr, ",", seqs$start, sep="")
        ans[seq(4, by=4, length.out=numChars)] = seqs$quality

        cat(ans, file=filename, sep="\n", append=TRUE)
    }
    return(filename)
}
