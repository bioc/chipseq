\name{diffPeakSummary}
\alias{diffPeakSummary}
\alias{diffPeakSummary,RleViewsList,RleViewsList-method}

\title{ A function to identify and produce summary statistics for differentially expressed peaks. }
\description{
  Given two sets of peaks, this function combines them and summarizes
  the individual coverage vectors under the combined peak set.
}
\usage{
diffPeakSummary(ranges1, ranges2,
                viewSummary = list(sums = viewSums, maxs = viewMaxs))
}
\arguments{
  \item{ranges1}{ First set of peaks (typically
  an \code{\link[IRanges:RleViewsList-class]{RleViewsList}}).}
  \item{ranges2}{ Second set of peaks (typically
  an \code{\link[IRanges:RleViewsList-class]{RleViewsList}}).}
  \item{viewSummary}{ A list of the per peak summary functions. }
}
\value{
  A \code{data.frame} with one row for each peak in the combined data.
  The chromosome, start and stop nucleotide positions (+ strand) are
  given as are the summary statistics requested.
}
\author{ D. Sarkar }
\examples{
data(cstest)
library(BSgenome.Mmusculus.UCSC.mm9)
seqlengths(cstest) <- seqlengths(Mmusculus)
## find peaks
findPeaks <- function(reads) {
  reads.ext <- resize(reads, width = 200)
  slice(coverage(reads.ext), lower = 8)
}
peakSummary <- diffPeakSummary(findPeaks(cstest$gfp), findPeaks(cstest$ctcf))
}
