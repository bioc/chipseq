

## functions to identify "differential peaks". See ../inst/Scripts/poisson.R

## symmetric version
## RG thinks it would be good to also have a diffSummary list that one
## could specify which summaries of the differences are wanted

## ML: it would be more flexible to have a function that takes a set
## of "interesting" intervals and summarizes them over the coverage of
## any number of samples. Actually, 'coverage' could be generalized to
## any quantitative track. Really, that would amount to summarizing a
## list of Views(List) objects and putting it all in one data frame.
## Ideally, we could have a MultiViews object, i.e., a single set of
## intervals on multiple subjects. Calling viewSums() would return a
## DataFrame, instead of a vector. Anyway, let us wait for a use case,
## for now, the chipseq package will do this:

setGeneric("diffPeakSummary",
           function(ranges1, ranges2,
                    viewSummary = list(sums = viewSums, maxs = viewMaxs))
           standardGeneric("diffPeakSummary"))

setMethod("diffPeakSummary", c("RleViewsList", "RleViewsList"), 
          function(ranges1, ranges2, 
                   viewSummary = list(sums = viewSums, maxs = viewMaxs))
{
    cov1 <- subject(ranges1)
    cov2 <- subject(ranges2)
    all.peaks <- union(ranges(ranges1), ranges(ranges2))
    peaks1 <- Views(cov1, all.peaks)
    peaks2 <- Views(cov2, all.peaks)
    peaks.comb <- Views(cov1 + cov2, all.peaks)

    ans <- GRanges(all.peaks)
    
    ans$comb.max <- unlist(viewMaxs(peaks.comb))
    
    if (is.list(viewSummary)) {
      for (nm in names(viewSummary)) {
          mcols(ans)[[paste0(nm, "1")]] <- unlist(viewSummary[[nm]](peaks1))
          mcols(ans)[[paste0(nm, "2")]] <- unlist(viewSummary[[nm]](peaks2))
      }
    }
    else {
      mcols(ans)[["summary1"]] <- unlist(viewSummary(peaks1))
      mcols(ans)[["summary2"]] <- unlist(viewSummary(peaks2))
    }

    ans
})


## version with respect to a reference

## diffPeakSummaryRef <-
##     function(obs.ranges, ref.ranges, chrom.lens,
##              lower = 10, extend = 0,
##              peak.fun = function(x) slice(x, lower = lower),
##              viewSummary = list(sums = viewSums, maxs = viewMaxs))

##     ## 'extend' is unused.  The intent is to extend the peaks by this
##     ## amount before summarizing

## {
##     if (!is(obs.ranges, "list")) obs.ranges <- as(obs.ranges, "list")
##     if (!is(ref.ranges, "list")) ref.ranges <- as(ref.ranges, "list")

##     obs.cov <- laneCoverage(obs.ranges, chrom.lens)
##     ref.cov <- laneCoverage(ref.ranges, chrom.lens)
##     ref.peaks <- lapply(ref.cov, peak.fun)
##     ref.peaks.in.obs <- copyIRangesbyChr(ref.peaks, obs.cov)

##     peakSummary <-
##         do.call(rbind,
##                 sapply(names(ref.peaks),
##                        function(chr) {
##                            ans <-
##                                data.frame(start = start(ref.peaks[[chr]]),
##                                           end = end(ref.peaks[[chr]]),
##                                           stringsAsFactors = FALSE)
##                            if (is.list(viewSummary))
##                            {
##                                for (nm in names(viewSummary))
##                                {
##                                    ans[[paste("ref", nm, sep = ".")]] <- viewSummary[[nm]](ref.peaks[[chr]])
##                                    ans[[paste("obs", nm, sep = ".")]] <- viewSummary[[nm]](ref.peaks.in.obs[[chr]])
##                                }
##                            }
##                            else 
##                            {
##                                ans[["ref.summary"]] <- viewSummary(ref.peaks[[chr]])
##                                ans[["obs.summary"]] <- viewSummary(ref.peaks.in.obs[[chr]])
##                            }
##                            ans
##                        },
##                        simplify = FALSE))
##     rownames(peakSummary) <- NULL
##     peakSummary
## }

## Subsampling two "lanes", so that on a per chromosome basis they
## have the same number of reads; let's not do this if we are
## reasonably close - when fudge

laneSubsample <- function(lane1, lane2, fudge = 0.05)
{
    idx1 <- split(seq_along(lane1), as.character(seqnames(lane1)))
    idx2 <- split(seq_along(lane2), as.character(seqnames(lane2)))
    keep <- intersect(names(idx1), names(idx2))
    idx1 <- idx1[keep]
    idx2 <- idx2[keep]

    len1 <- sapply(idx1, length)
    len2 <- sapply(idx2, length)

    samp <- names(idx1)[abs(len1 - len2) >= fudge * len1]
    for (i in samp) {
        if (len1[i] < len2[i])
            idx2[[i]] <- sample(idx2[[i]], len1[i])
        else
            idx1[[i]] <- sample(idx1[[i]], len2[i])
    }

    GRangesList(lane1=lane1[sort(unlist(idx1, use.names=FALSE))],
                lane2=lane2[sort(unlist(idx2, use.names=FALSE))])
}
