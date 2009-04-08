\name{diffPeakSummary}
\alias{diffPeakSummary}
\title{ A function to identify and produce summary statistics for differentially expressed peaks. }
\description{
  Given two sets of reads this function identifies all peaks in the
  combined data with height larger than \code{lower} and then uses those
  regions to compute summary statistics for each of the sets separately.
}
\usage{
diffPeakSummary(ranges1, ranges2, chrom.lens,
                lower = 10, extend = 0, 
                peak.fun = NULL, merge = 0L, islands = FALSE,
                viewSummary = list(sums = viewSums, maxs = viewMaxs))
}
\arguments{
  \item{ranges1}{ First set of reads (as IRanges).}
  \item{ranges2}{ Second set of reads (as IRanges). }
  \item{chrom.lens}{ The lengths of the chromosomes for the organism. }
  \item{lower}{ The height used to declare a peak in the combined samples. }
  \item{extend}{ Currently unused. The intent is to extend peaks by this
    amount before summarizing. }
  \item{peak.fun}{ Function use to use to find peaks.  The default makes
    use of additional arguments \code{merge} and \code{islands}, which
    are otherwise ignored. }
  \item{merge}{ Integer giving the amount of gaps between peaks that
    shoudl be conidered significant.  Smaller gaps are removed by
    combining neighbouring peaks.  }
  \item{islands}{ Logical indicating whether or not to use islands
    (coverage > 0) as peaks. }
  \item{viewSummary}{ A list of the per peak summaries. }
}
\value{
  A \code{data.frame} with one row for each peak in the combined data.
  The chromosome, start and stop nucleotide positions (+ strand) are
  given as are the summary statistics requested.
}
\author{ D. Sarkar }
\examples{
}
