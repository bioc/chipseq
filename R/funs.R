##FIXME: exclude should be used, and the chromosomeFilter

readReads <-
    function(srcdir, lane, ...,
             include = "chr[0-9]+$", type = "MAQMapShort",
             simplify = TRUE, minScore=15)
{
    filt <-
        compose(strandFilter(strandLevels=c("-", "+")),
                chromosomeFilter(regex = include),
                uniqueFilter(withSread = FALSE),
                alignQualityFilter(minScore),
                ...)
    message(sprintf("reading data from lane %s [%s], using filter %s", lane, 
         srcdir, name(filt)))
    ans <- readAligned(srcdir, lane, type = type, filter = filt)
    if (simplify) as(ans, "GenomeData")
    else ans
}


## summarize an "AlignedRead" as a nested list; may be useful as an
## intermediate form for storage etc.  Return value is of the form:
## 
## List of 20
##  $ chr1        :List of 2
##   ..$ -: int [1:88596] 3011060 3023510 3024278 ...
##   ..$ +: int [1:88333] 3025629 3034277 3058436 ...
##  $ chr2        :List of 2
##   ..$ -: int [1:94326] 3006779 3010130 3015418 ...
##   ..$ +: int [1:95154] 3009910 3010789 3017043 ...
##  $ ...

setMethod("as.list", "AlignedRead",
          function(x, ...) {
              readStart <- ifelse(strand(x) == "-",
                                  position(x) + width(x) - 1L,
                                  position(x))
              alignLocs <-
                  split(data.frame(position = readStart, strand = strand(x)),
                        chromosome(x)[drop=TRUE])
              lapply(alignLocs,
                     function(df) with(df, split(position, strand))[c("-", "+")])
          })

setAs("AlignedRead", "GenomeData",
      function(from) {
          GenomeData(as.list(from))
      })




## take aligned reads and grow them appropriately.  FIXME: how should
## we support list version? methods?

splitbyChr <- function(x) {
    split(position(x), paste(chromosome(x), strand(x), sep=""))
}

extendReads <- function(reads, seqLen=200, strand = c("+", "-"))
{
    if (is(reads, "GenomeDataList") || is(reads, "GenomeData"))
        gdapply(reads, extendReads, seqLen = seqLen, strand = strand)
    else if (is(reads, "AlignedRead")) 
    {
        rng <- IRanges(ifelse(strand(reads) == "+",
                              position(reads),
                              position(reads) + width(reads) - seqLen),
                       ifelse(strand(reads) == "+",
                              position(reads) + seqLen - 1L,
                              position(reads) + width(reads) - 1L))
        strandIdx <- strand(reads) %in% strand
        split(rng[strandIdx], chromosome(reads)[strandIdx])
    }
    else if (is.list(reads)) 
    {
        seqLen <- as.integer(seqLen)
        if (all(c("+", "-") %in% names(reads))) 
        {
            reads <- reads[strand]
            starts <- as.integer(c(reads[["+"]], reads[["-"]] - seqLen + 1L))
            new("IRanges", start = starts, width = rep(seqLen, length.out = length(starts)), NAMES = NULL)
            ## IRanges(start = c(reads[["+"]], reads[["-"]] - seqLen + 1L),
            ##        end = c(reads[["+"]] + seqLen - 1L, reads[["-"]]))
            ## width = seqLen) # error when start=numeric(0)
        }
    }
    else stop("Invalid value for 'reads'")
}



## take some sort of a view and copy it to a new vector.  FIXME: can
## we improve the interface?


copyIRanges <- function(IR1, newX) Views(newX, IR1)

copyIRangesbyChr <- function(IR1, newX) {
    nms = names(IR1)
    ans = vector("list", length=length(nms))
    names(ans) = nms
    for(i in nms)
        ans[[i]] = Views(newX[[i]], IR1[[i]])
    ans
}





## Subsampling two "lanes", so that on a per chromosome basis they
## have the same number of reads; let's not do this if we are
## reasonably close - when fudge

laneSubsample <- function(lane1, lane2, fudge = 0.05)
{
    chromList = names(lane1) ##lane2 should have the same names
    l1Len = unlist(lapply(lane1, length)) # sapply doesn't work for "GenomeData"
    l2Len = unlist(sapply(lane2, length)) 
    for(i in seq_len(length(l1Len)))
    {
        if(abs(l1Len[i]-l2Len[i])/l1Len[i] < fudge) next
        if(l1Len[i] < l2Len[i])
            lane2[[i]] <- sample(lane2[[i]], l1Len[i])
        if(l1Len[i] > l2Len[i])
            lane1[[i]] <- sample(lane1[[i]], l2Len[i])
    }
    GenomeDataList(list(lane1=lane1, lane2=lane2))
}




## Take an IRanges object and merge all Ranges that are less than
## maxgap apart.  FIXME: make this a merge() method?

merge <- function(IR, maxgap)
{
    if (!(isNormal(IR))) stop("IR must be a normal IRanges object")
    width(IR) <- width(IR) + maxgap
    reduce(IR) 
}

## Take two or more lanes and combine the read positions.  FIXME: This
## should be unified with combineLanes below (which works on extended
## reads)

combineLaneReads <- function(laneList, chromList = names(laneList[[1]]), keep.unique = FALSE)
{
    my.unlist <-
        if(keep.unique) function(x, ...) { unique(unlist(x, ...)) }
        else unlist
    combineReads <- function(chr)
    {
        list("+" = my.unlist(lapply(laneList, function(x) x[[chr]][["+"]]), use.names = FALSE),
             "-" = my.unlist(lapply(laneList, function(x) x[[chr]][["-"]]), use.names = FALSE))
    }
    names(chromList) = chromList ##to get the return value named
    GenomeData(lapply(chromList, combineReads))
}

## Take two or more lanes worth of extended reads and combine.

combineLaneRanges <- function(laneList, chromList = names(laneList[[1]]), keep.unique = FALSE)
{
    combineRanges <- function(chr)
    {
        ans <-
            IRanges(start = unlist(lapply(laneList, function(x) start(x[[chr]])), use.names = FALSE),
                    end = unlist(lapply(laneList, function(x) end(x[[chr]])), use.names = FALSE))
        if (keep.unique) ans[!duplicated(ans)]
        else ans
    }
    names(chromList) = chromList ##to get the return value named
    GenomeData(lapply(chromList, combineRanges))
}


## exported interface
combineLanes <- function(x, chromList = names(x[[1]]), keep.unique = FALSE)
{
    isRange <- is(x[[1]][[chromList[1]]], "IRanges")
    if (isRange)
        combineLaneRanges(laneList = x, chromList = chromList, keep.unique = keep.unique)
    else
        combineLaneReads(laneList = x, chromList = chromList, keep.unique = keep.unique)
}



