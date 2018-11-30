
## FIXME: should this go into IRanges? (rename to viewSummary?)
setGeneric("peakSummary", function(x, ...) standardGeneric("peakSummary"))

setMethod("peakSummary", "RleViews", function(x) {
  GRanges(ranges(x), max = viewMaxs(x), maxpos = mid(viewRangeMaxs(x)),
          sum = viewSums(x))
})
setMethod("peakSummary", "RleViewsList", function(x) {
  summaries <- unname(lapply(x, peakSummary))
  rd <- do.call(c, summaries)
  names(rd) <- names(x)[elementNROWS(summaries) > 0]
  rd
})
