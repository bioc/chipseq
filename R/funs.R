
##Fixme: a lot of this functionality will become available in 
##ReadAligned from the ShortRead package
readAndClean <-
    function(maqDir, pattern, exclude="random", 
             dropDups = TRUE, minScore = 15,
             type = "MAQMap") 
{
    s1 <- readAligned(maqDir, pattern = pattern, type = type)
    exChr = grep(exclude, s1@chromosome)
    s1 = s1[-exChr]
    keep = (s1@alignQuality@quality >= minScore)
    s2 = s1[keep]
    if( dropDups )
        return(s2[!srduplicated(sread(s2))])
    s2
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
              alignLocs <-
                  split(data.frame(position = x@position, strand = x@strand),
                        x@chromosome[drop=TRUE])
              lapply(alignLocs,
                     function(df) with(df, split(position, strand))[c("-", "+")])
          })




splitbyChr <- function(x) {
    split(x@position, paste(x@chromosome, x@strand, sep=""))
}


## take input MAQ mapped reads and grow them appropriately.
## FIXME: how should we support list version? methods?
extendReads <- function(reads, readLen=35, seqLen=200,
                        strand = c("+", "-"))
{
    if (is(reads, "AlignedRead")) 
    {
        s1 <- splitbyChr(reads)
        chrs <- unique(gsub("\\+|-", "", names(s1)))
        byCIR <- vector("list", length=length(chrs))
        names(byCIR) <- chrs
        for(i in chrs) {
            ipos <-
                if ("+" %in% strand) grep(paste("^", i, "\\+$", sep=""), names(s1))
                else character(0)
            ineg <-
                if ("-" %in% strand) grep(paste("^", i, "-$", sep=""), names(s1))
                else character(0)
            byCIR[[i]] = IRanges(c(s1[[ipos]], s1[[ineg]]-seqLen+readLen),
                                 c(s1[[ipos]]+seqLen-1, s1[[ineg]]+readLen-1))
            byCIR
        }
        byCIR
        ## or, extendReads(as.list(reads), readLen=readLen, seqLen=seqLen)
    }
    else if (is.list(reads)) 
    {
        if (all(c("+", "-") %in% names(reads))) 
        {
            reads <- reads[strand]
            IRanges(start = c(reads[["+"]], reads[["-"]] - seqLen + readLen),
                    end = c(reads[["+"]] + seqLen - 1, reads[["-"]] + readLen - 1))
            ## width = seqLen) # error when start=numeric(0)
        }
        else lapply(reads, extendReads)
    }
    else stop("Invalid value for 'reads'")
}

## alias, remove later
growSeqs <- function(...) extendReads(...)

## IRanges(start = c(x[["+"]], x[["-"]] - seqLen + readLen),
##         end = c(x[["+"]] + seqLen - 1, x[["-"]] + readLen - 1),
##         width = seqLen)








## FIXME: make this a merge() method?
#take an IRanges object and merge all Ranges that
#are less than maxgap apart.
merge <- function(IR, maxgap)
{
  if(!(isNormal(IR))) stop("IR must be a normal IRanges object")
  width(IR) = width(IR)+maxgap
  reduce(IR) 
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
    IRanges(start = unlist(lapply(rlist, function(x) start(x[[byChr]])),
                           use.names = FALSE),
            end = unlist(lapply(rlist, function(x) end(x[[byChr]])),
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
            lane1[[i]] = sample(lane1[[i]], l2Len[i])
    }
    return(list(lane1=lane1, lane2=lane2))
}


laneCoverage <- function(lane, chromLens) {
    sapply(names(lane),
           function(chr) {
               coverage(lane[[chr]], 1, chromLens[chr])
           }, 
           simplify = FALSE)
}

islands <- function(x) lapply(x, slice, lower = 1)

readsPerIsland <- function(lane, ntperread = 200L) 
        lapply(lane,
               function(x) viewSums(x) / ntperread)

maxHeightPerIsland <- function(lane) lapply(lane, viewMaxs)

islandSummary <- function(x, ntperread = 200L)
{
    a1 = sapply(names(x), function(chr) {
        ans <- data.frame(start = start(x[[chr]]),
                          end = end(x[[chr]]),
                          width = width(x[[chr]]),
                          clones = viewSums(x[[chr]])/ntperread,
                          maxdepth = viewMaxs(x[[chr]]))
        if (is.unsorted(ans$start))
            stop("Starts not sorted! Check code.")
        n <- nrow(ans)
        ans$inter <- with(ans, c(NA, start[-1] - end[-n]))
        ans },
         simplify = FALSE)
    names(a1) = names(x)
    a1
}

## take some sort of a view and copy it to a new vector
copyIRanges <- function(IR1, newX)
    Views(newX, start=start(IR1), end=end(IR1))

copyIRangesbyChr <- function(IR1, newX) {
    nms = names(IR1)
    ans = vector("list", length=length(nms))
    names(ans) = nms
    for(i in nms)
        ans[[i]] = copyIRanges(IR1[[i]], newX[[i]])
    ans
}
