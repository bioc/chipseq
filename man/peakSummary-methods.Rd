\name{peakSummary-methods}
\alias{peakSummary-methods}
\alias{peakSummary}
\alias{peakSummary,RleViews-method}
\alias{peakSummary,RleViewsList-method}
\title{ Summarizing peak sets }
\description{
  Summarizes a set of peaks into a
  \code{\link[GenomicRanges:GRanges-class]{GRanges}} object with columns
  of statistics like the peak maxima and integrals (sums).
}
\usage{
peakSummary(x, ...) 
}
\arguments{
  \item{x}{An object containing peaks, usually
    a \code{\link[IRanges:RleViewsList-class]{RleViewsList}}.}
  \item{...}{Arguments to pass to methods}
}

\value{
  A \code{GRanges} object of the peaks, with columns named
  \code{max}, \code{maxpos} (position of the maximum, centered),
  and \code{sum}.
}

\seealso{
  \link[IRanges]{view-summarization-methods} in the IRanges package
  for view summarization methods like \code{viewMaxs} and \code{viewSums}.
}

\keyword{methods}
