
 splitbyChr <- function(x) {
    split(x@position, paste(x@chromosome, x@strand, sep=""))
}

##take input MAQ mapped reads and grow them appropriately
 growSeqs <- function(reads, readLen=35, seqLen=200 )
{
    s1 <- splitbyChr(reads)
    chrs = unique(gsub("\\+|-", "", names(s1)))
    byCIR <- vector("list", length=length(chrs))
    names(byCIR) <- chrs
    for(i in chrs) {
        ipos = grep(paste("^", i, "\\+$", sep=""), names(s1))
        ineg = grep(paste("^", i, "-$", sep=""), names(s1))
        byCIR[[i]] = IRanges(c(s1[[ipos]], s1[[ineg]]-seqLen-readLen),
              c(s1[[ipos]]+seqLen-1, s1[[ineg]]+readLen-1))
        byCIR
    }
    byCIR
}

combineRanges <- function(rlist)
{
    IRanges(start = unlist(lapply(rlist, start), use.names = FALSE),
            end = unlist(lapply(rlist, end), , use.names = FALSE))
}


##FIXME: a lot of this functionality will become available in 
##ReadAligned from the ShortRead package
readAndClean <- function(maqDir, pattern = pat, exclude="random", 
   dropDups = TRUE, minScore = 15 ) 
{
    s1 <- readAligned(maqDir, pattern = pat, type = "MAQMap")
    exChr = grep(exclude, s1@chromosome)
    s1 = s1[-exChr]
    keep = (s1@alignQuality@quality >= minScore)
    s2 = s1[keep]
    if( dropDups )
        return(s2[!srduplicated(sread(s2))])
    s2
}

