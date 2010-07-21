

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
          ## 
          ## 'extend' is unused.  The intent is to extend the peaks by this
          ## amount before summarizing
          ##
{
    cov1 <- subject(ranges1)
    cov2 <- subject(ranges2)
    all.peaks <- union(ranges1, ranges2)
    peaks1 <- Views(cov1, all.peaks)
    peaks2 <- Views(cov2, all.peaks)
    peaks.comb <- Views(cov1 + cov2, all.peaks)

    ans <- RangedData(all.peaks)
    
    ans$comb.max <- unlist(viewMaxs(peaks.comb))
    
    if (is.list(viewSummary)) {
      for (nm in names(viewSummary)) {
        ans[[paste(nm, "1", sep = "")]] <- unlist(viewSummary[[nm]](peaks1))
        ans[[paste(nm, "2", sep = "")]] <- unlist(viewSummary[[nm]](peaks2))
      }
    }
    else {
      ans[["summary1"]] <- unlist(viewSummary(peaks1))
      ans[["summary2"]] <- unlist(viewSummary(peaks2))
    }

    ans
})


## version with respect to a reference

diffPeakSummaryRef <-
    function(obs.ranges, ref.ranges, chrom.lens,
             lower = 10, extend = 0,
             peak.fun = function(x) slice(x, lower = lower),
             viewSummary = list(sums = viewSums, maxs = viewMaxs))

    ## 'extend' is unused.  The intent is to extend the peaks by this
    ## amount before summarizing

{
    if (!is(obs.ranges, "list")) obs.ranges <- as(obs.ranges, "list")
    if (!is(ref.ranges, "list")) ref.ranges <- as(ref.ranges, "list")

    obs.cov <- laneCoverage(obs.ranges, chrom.lens)
    ref.cov <- laneCoverage(ref.ranges, chrom.lens)
    ref.peaks <- lapply(ref.cov, peak.fun)
    ref.peaks.in.obs <- copyIRangesbyChr(ref.peaks, obs.cov)

    peakSummary <-
        do.call(rbind,
                sapply(names(ref.peaks),
                       function(chr) {
                           ans <-
                               data.frame(start = start(ref.peaks[[chr]]),
                                          end = end(ref.peaks[[chr]]),
                                          stringsAsFactors = FALSE)
                           if (is.list(viewSummary))
                           {
                               for (nm in names(viewSummary))
                               {
                                   ans[[paste("ref", nm, sep = ".")]] <- viewSummary[[nm]](ref.peaks[[chr]])
                                   ans[[paste("obs", nm, sep = ".")]] <- viewSummary[[nm]](ref.peaks.in.obs[[chr]])
                               }
                           }
                           else 
                           {
                               ans[["ref.summary"]] <- viewSummary(ref.peaks[[chr]])
                               ans[["obs.summary"]] <- viewSummary(ref.peaks.in.obs[[chr]])
                           }
                           ans
                       },
                       simplify = FALSE))
    rownames(peakSummary) <- NULL
    peakSummary
}

## Subsampling two "lanes", so that on a per chromosome basis they
## have the same number of reads; let's not do this if we are
## reasonably close - when fudge

laneSubsample <- function(lane1, lane2, fudge = 0.05)
{
  .quickAndDirtyCoercionFromGRangesToGenomeData <- function(x)
  {
    if ("*" %in% unique(strand(x)))
        stop("cannot coerce 'x' if 'strand(x)' contains \"*\"")
    y <- split(x, seqnames(x))
    listData <- lapply(seq_len(length(y)),
                    function(i) {
                        split(start(y[[i]]),
                              as.character(strand(y[[i]])))
                    })
    names(listData) <- names(y)
    GenomeData(listData)
  }
  if (is(lane1, "GRanges"))
    lane1 <- .quickAndDirtyCoercionFromGRangesToGenomeData(lane1)
  if (is(lane2, "GRanges"))
    lane2 <- .quickAndDirtyCoercionFromGRangesToGenomeData(lane2)
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
                    
