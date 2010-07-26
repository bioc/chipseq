\name{peakCutoff}
\alias{peakCutoff}
\title{
  Calculate a peak cutoff
}
\description{
  Calculates a peak cutoff value given an FDR, assuming a Poisson noise
  distribution estimated from the frequency of singleton and doubleton
  islands.
}
\usage{
peakCutoff(cov, fdr.cutoff = 0.001, k = 2:20)
}
\arguments{
  \item{cov}{
    The coverage object, e.g.,
    an \code{\link[IRanges:RleList-class]{RleList}} object.
  }
  \item{fdr.cutoff}{
    The maximum-allowed FDR for calculating the cutoff.
  }
  \item{k}{
    The coverage levels at which to estimate an FDR value. The maximal
    value that is less than \code{fdr.cutoff} is chosen for calculating
    the cutoff. Usually best left to the default.
  }
}
\value{
  A numeric value to use for calling peaks
}
\author{
  D. Sarkar and M. Lawrence
}

\seealso{
  \code{\link{islandDepthPlot}} for the graphical equivalent; the
  vignette for a bit more explanation.
}

\examples{
data(cstest)
cov <- coverage(resize(cstest$ctcf, width=200))
peakCutoff(cov)
}
