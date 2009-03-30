

## functions to identify "differential peaks". See ../inst/Scripts/poisson.R

## symmetric version
## RG thinks it would be good to also have a diffSummary list that one
## could specify which summaries of the differences are wanted


laneCoverage <- function(lane, chromLens) {
    sapply(names(lane),
           function(chr) {
               coverage(lane[[chr]], 1, chromLens[chr])
           }, 
           simplify = FALSE)
}


diffPeakSummary <-
    function(ranges1, ranges2, chrom.lens,
             lower = 10, extend = 0, islands = FALSE,
             peak.fun = NULL,
             viewSummary = list(sums = viewSums, maxs = viewMaxs))

    ## 'extend' is unused.  The intent is to extend the peaks by this
    ## amount before summarizing

{
    if (!is(ranges1, "list")) ranges1 <- as(ranges1, "list")
    if (!is(ranges2, "list")) ranges2 <- as(ranges2, "list")

    if (is.null(peak.fun)) 
        peak.fun <- function(x)
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
    comb.peaks <- lapply(comb.cov, peak.fun)

    cov1 <- laneCoverage(ranges1, chrom.lens)
    cov2 <- laneCoverage(ranges2, chrom.lens)
    peaks1 <- copyIRangesbyChr(comb.peaks, cov1)
    peaks2 <- copyIRangesbyChr(comb.peaks, cov2)

    peakSummary <-
        do.call(rbind,
                lapply(names(peaks1),
                       function(chr) {
                           if (length(peaks1[[chr]]) == 0) return(NULL)
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


                    
