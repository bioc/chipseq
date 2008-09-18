
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

##take two or more lanes and combine the data
combineLanes <- function(laneList, chromList = names(laneList[[1]])) {
    names(chromList) = chromList ##to get the return value named
    lapply(chromList, function(x) combineRanges(laneList, x))
}

##one might also want a one argument version, but for now, this seems
##to be somewhat effective
combineRanges <- function(rlist, byChr)
{
    IRanges(start = unlist(lapply(rlist, function(x)
            start(x[[byChr]])), use.names = FALSE),
            end = unlist(lapply(rlist, function(x)  end(x[[byChr]])),
            use.names = FALSE))
}

##subsampling two "lanes", so that on a per chromosome basis they have
##the same number of reads; let's not do this if we are reasonably
##close - when fudge
laneSubsample <- function(lane1, lane2, fudge = 0.05)
{
    chromList = names(lane1) ##lane2 should have the same names
    l1Len = sapply(lane1, length)
    l2Len = sapply(lane2, length)
    for(i in 1:length(l1Len) ) {
        if(abs(l1Len[i]-l2Len[i])/l1Len[i] < fudge) next
        if(l1Len[i] < l2Len[i])
            lane2[[i]] = sample(lane2[[i]], l1Len[i])
        if(l1Len[i] > l2Len[i])
            lane1[[i]] = sample(lane2[[i]], l2Len[i])
    }
    return(list(lane1=lane1, lane2=lane2))
}

##Fixme: a lot of this functionality will become available in 
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

