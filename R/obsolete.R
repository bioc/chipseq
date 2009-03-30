

if (FALSE)
{
    
    coverageperChromosome <- function (x, end, strand = c("+", "-")) 
    {
        g <- extendReads(x, strand = strand)
        if( missing(end) || is.na(end) ) max = max(end(g)) + 400L
        else max = end
        coverage(g, 1, max)
    }

    ## specialized coverage methods for our containers

    setMethod("coverage", "GenomeList",
              function(x, start, end, strand) {
                  ans = lapply(x, coverageperChromosome, end=end,
                               strand = strand)
                  new("GenomeList", ans, genome=x@genome)
              })

    setMethod("coverage", "AlignedList",
              function(x, start, end, strand) {
                  ans = lapply(x, coverage, start, end, strand)
                  new("AlignedList", ans)
              })


    ## common infrastructure to summarize reads in the form of nested lists: 
    ##  reads.list=list("1" = list("chr1" = list("+"=..., "-"=...), 
    ##                             "chr2"=...), 
    ##                  "2" = list(...))

    ## NOTE: summary.fun must return a data frame in the current
    ## implementation (and shouldn't be very big at that)

    summarizeLane <- function(clist, summary.fun, ...)
    {
        ## clist is a list at the lane level, with one list("+"=, "-"=) for each chromsome
        ans <- do.call(lattice::make.groups, lapply(clist, summary.fun, ...))
        names(ans)[names(ans) == "which"] <- "chromosome"
        ## cbind(chr = factor(colnames(ans), levels = colnames(ans)), as.data.frame(t(ans)))
        ans
    }

    ## Note: uses of summarizeLane should be replaced by
    ## as(gdApply(clist, FUN = summary.fun, ...), "data.frame")

    summarizeLane2 <- function(clist, summary.fun, ..., seqlen)
    {
        ## clist is a list at the lane level, with one list("+"=, "-"=) for each chromsome
        stopifnot(all(names(clist) %in% names(seqlen)))
        seqlen <- seqlen[names(clist)]
        mapply(summary.fun, clist, seqlen, ..., SIMPLIFY = FALSE)
    }

    summarizeReads <- 
        function(reads.list, lanes = names(reads.list), ..., verbose = FALSE)
        {
            if (verbose) cat(paste("Processing lanes", paste(lanes, collapse = ",")), fill = TRUE)
            ans <- do.call(lattice::make.groups, lapply(reads.list[lanes], summarizeLane, ...))
            names(ans)[names(ans) == "which"] <- "lane"
            ans
        }

    sliceSummary <- 
        function(x, lower = 1,
                 viewSummary = list(sums = viewSums, maxs = viewMaxs))
            ## x is a list at the lane->chromosome level, with components "+" and "-"
        {
            g <- growSeqs(x)
            cov <- coverage(g, 1, max(end(g) + 400L))
            s <- slice(cov, lower = lower)
            ans <- data.frame(start = start(s), end = end(s))
            if (is.list(viewSummary)) 
            {
                for (nm in names(viewSummary))
                    ans[[nm]] <- viewSummary[[nm]](s)
            }
            else ans[["summary"]] <- viewSummary(s)
            ans
        }

    ## ???
    count.singletons <- function(x) ## x <- summarizeReads(.)
    {
        ans <-
            do.call(make.groups,
                    lapply(x,
                           function(x) {
                               as.data.frame.table(sapply(x, length))
                           }))
        names(ans)[names(ans) == "Freq"] <- "count"
        names(ans)[names(ans) == "Var1"] <- "chromosome"
        names(ans)[names(ans) == "which"] <- "lane"
        ans
    }



    readAndClean <-
        function(maqDir, pattern, exclude="random", 
                 dropDups = TRUE, minScore = 15,
                 type = "MAQMap") 
        {
            .Deprecated("readReads")
            s1 <- readAligned(maqDir, pattern = pattern, type = type)
            exChr = grep(exclude, s1@chromosome)
            s1 = s1[-exChr]
            keep = (s1@alignQuality@quality >= minScore)
            s2 = s1[keep]
            if( dropDups )
                return(s2[!srduplicated(sread(s2))])
            s2
        }

    ## used in MyoD.R
    ##takes a GenomeList as input, and returns one

    islands <- function(x) {
        ans = lapply(x, slice, lower = 1)
        if( is(x, "GenomeList") )
            new("GenomeList", ans, genome = x@genome)
        else ans
    }

    
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
    
    
}
