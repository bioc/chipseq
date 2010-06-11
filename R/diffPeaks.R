

## functions to identify "differential peaks". See ../inst/Scripts/poisson.R

## symmetric version
## RG thinks it would be good to also have a diffSummary list that one
## could specify which summaries of the differences are wanted


laneCoverage <- function(lane, chromLens) {
    sapply(names(lane),
           function(chr) {
               coverage(lane[[chr]], width = chromLens[chr])
           }, 
           simplify = FALSE)
}


countOverlappingReads <- function(peaks, reads)
{
    as.numeric(as.table(t(findOverlaps(sort(reads), peaks))))
}



diffPeakSummary <-
    function(ranges1, ranges2, chrom.lens,
             lower = 10, extend = 0, 
             peak.fun = NULL, merge = 0L, islands = FALSE,
             viewSummary = list(sums = viewSums, maxs = viewMaxs))
    ## 
    ## 'extend' is unused.  The intent is to extend the peaks by this
    ## amount before summarizing
    ##
{
    if (!is(ranges1, "list")) ranges1 <- as(ranges1, "list")
    if (!is(ranges2, "list")) ranges2 <- as(ranges2, "list")
    if (is.null(peak.fun)) 
        peak.fun <- function(x)
        {
            peaks <-
                if (islands)
                {
                    s <- slice(x, lower = 1)
                    s[viewMaxs(s) >= lower]
                }
                else 
                    slice(x, lower = lower)
            if (merge > 0)
            {
                end(peaks) <- end(peaks) + merge
                peaks <- reduce(peaks)
                end(peaks) <- end(peaks) - merge
            }
            peaks
        }
    combined <- combineLanes(list(ranges1, ranges2))
    comb.cov <- laneCoverage(combined, chrom.lens)
    comb.peaks <- lapply(comb.cov, peak.fun)
    cov1 <- laneCoverage(ranges1, chrom.lens)
    cov2 <- laneCoverage(ranges2, chrom.lens)
    peaks1 <- copyIRangesbyChr(comb.peaks, cov1)
    peaks2 <- copyIRangesbyChr(comb.peaks, cov2)
    peakSummary <-
        do.call(rbind,
                lapply(names(comb.peaks),
                       function(chr) {
                           if (length(peaks1[[chr]]) == 0) return(NULL)
                           ans <-
                               data.frame(chromosome = chr,
                                          start = start(comb.peaks[[chr]]),
                                          end = end(comb.peaks[[chr]]),
                                          comb.max = viewMaxs(comb.peaks[[chr]]),
                                          overlap1 = countOverlappingReads(comb.peaks[[chr]], ranges1[[chr]]),
                                          overlap2 = countOverlappingReads(comb.peaks[[chr]], ranges2[[chr]]),
                                          stringsAsFactors = FALSE)
                           if (is.list(viewSummary))
                           {
                               for (nm in names(viewSummary))
                               {
                                   ans[[paste(nm, "1", sep = "")]] <- viewSummary[[nm]](peaks1[[chr]])
                                   ans[[paste(nm, "2", sep = "")]] <- viewSummary[[nm]](peaks2[[chr]])
                               }
                           }
                           else 
                           {
                               ans[["summary1"]] <- viewSummary(peaks1[[chr]])
                               ans[["summary2"]] <- viewSummary(peaks2[[chr]])
                           }
                           ans
                       }))
    ## names(peakSummary)[names(peakSummary) == "which"] <- "chromosome"
    rownames(peakSummary) <- NULL
    peakSummary
}


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
                    
