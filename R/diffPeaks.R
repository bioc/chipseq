

## functions to identify "differential peaks". See ../inst/Scripts/poisson.R

## symmetric version
## RG thinks it would be good to also have a diffSummary list that one
## could specify which summaries of the differences are wanted
diffPeakSummary <-
    function(ranges1, ranges2, chrom.lens,
             lower = 10, extend = 0, islands = TRUE,
             viewSummary = list(sums = viewSums, maxs = viewMaxs))

    ## 'extend' is unused.  The intent is to extend the peaks by this
    ## amount before summarizing

{
    doSlice <- function(x, lower)
    {
        if (islands)
        {
            s <- slice(x, lower = 1)
            s[viewMaxs(s) >= lower]
        }
        else 
            slice(x, lower = lower)
    }
    
    combined <- combineLanes(list(ranges1, ranges2))
    comb.cov <- laneCoverage(combined, chrom.lens)
    comb.peaks <- lapply(comb.cov, doSlice, lower = lower)

    cov1 <- laneCoverage(ranges1, chrom.lens)
    cov2 <- laneCoverage(ranges2, chrom.lens)
    peaks1 <- copyIRangesbyChr(comb.peaks, cov1)
    peaks2 <- copyIRangesbyChr(comb.peaks, cov2)

    peakSummary <-
        do.call(rbind,
                lapply(names(peaks1),
                       function(chr) {
                           ans <-
                               data.frame(chromosome = chr,
                                          start = start(peaks1[[chr]]),
                                          end = end(peaks1[[chr]]),
                                          comb.max = viewMaxs(comb.peaks[[chr]]),
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
             viewSummary = list(sums = viewSums, maxs = viewMaxs))

    ## 'extend' is unused.  The intent is to extend the peaks by this
    ## amount before summarizing

{
    obs.cov <- laneCoverage(obs.ranges, chrom.lens)
    ref.cov <- laneCoverage(ref.ranges, chrom.lens)
    ref.peaks <- lapply(ref.cov, slice, lower = lower) # extend=extend (FIXME: unused)
    ref.peaks.in.obs <- copyIRangesbyChr(ref.peaks, obs.cov)

    peakSummary <-
        do.call(make.groups,
                sapply(names(ref.peaks),
                       function(chr) {
                           ans <-
                               data.frame(start = start(ref.peaks[[chr]]),
                                          end = end(ref.peaks[[chr]]))
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
    names(peakSummary)[names(peakSummary) == "which"] <- "chromosome"
    peakSummary
}


                    
