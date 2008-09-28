

## functions to identify "differential peaks". See ../inst/Scripts/poisson.R

diffPeakSummary <-
    function(obs.ranges, ref.ranges, chrom.lens,
             lower = 10, extend = 0,
             viewSummary = list(sums = viewSums, maxs = viewMaxs))

    ## 'extend' is unused.  The intent is to extend the peaks by this
    ## amount before summarizing

{
    obs.cov <- laneCoverage(obs.ranges, chromLens)
    ref.cov <- laneCoverage(ref.ranges, chromLens)
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


                    
