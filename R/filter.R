## A srFilter for use with e.g. readAligned()
chipseqFilter <- function(exclude = "[_MXY]",
                          uniqueness = c("location", "sequence",
                            "location*sequence", "none"),
                          hasStrand = TRUE)
{
  if (!is.null(exclude) && (!is.character(exclude) || any(is.na(exclude))))
    stop("'exclude' must be character without NA's")
  uniqueness <- match.arg(uniqueness)
  if (!isTRUEorFALSE(hasStrand))
    stop("'hasStrand' must be TRUE or FALSE")
  filt <- srFilter()
  if (hasStrand)
    filt <- compose(filt, strandFilter(strandLevels=c("-", "+")))
  if (!is.null(exclude))
    filt <- compose(filt, chromosomeFilter(exclude, exclude = TRUE))
  if (uniqueness != "none") {
    withSread <- switch(uniqueness, location = FALSE, sequence = NA,
                        `location*sequence` = TRUE)
    ofilt <-
      occurrenceFilter(withSread = withSread)
    filt <- compose(filt, ofilt)
  }
  filt
}
