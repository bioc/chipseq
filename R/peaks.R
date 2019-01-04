
## FIXME: should this go into IRanges? (rename to viewSummary?)
setGeneric("peakSummary", function(x, ...) standardGeneric("peakSummary"))

setMethod("peakSummary", "RleViews", function(x) {
    ir <- ranges(x)
    mcols(ir) <- DataFrame(max = viewMaxs(x),
                           maxpos = mid(viewRangeMaxs(x)),
                           sum = viewSums(x))
    ir
})
setMethod("peakSummary", "RleViewsList", function(x) {
  summaries <- unname(lapply(x, peakSummary))
  ir <- do.call(c, summaries)
  GRanges(rep(names(x), lengths(summaries)), ir)
})
