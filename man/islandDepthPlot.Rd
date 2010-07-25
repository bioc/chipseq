\name{islandDepthPlot}
\alias{islandDepthPlot}
\title{
Plot island depth distribution
}
\description{
  Plots the distribution of island depths using points for the observed
  islands and a line for the Poisson estimate of the noise. Useful for
  choosing a depth corresponding to a desired FDR.
}
\usage{
islandDepthPlot(x, maxDepth = 20L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    A coverage object, e.g., \code{\link[IRanges:RleList-class]{RleList}}.
  }
  \item{maxDepth}{
    The maximum depth to plot (there are usually some outliers).
  }
}
\author{
D. Sarkar, M. Lawrence
}

\seealso{
\code{\link{peakCutoff}} for calculating a cutoff value for an FDR.
}
\examples{
data(cstest)
cov <- coverage(resize(cstest$ctcf, width=200))
islandDepthPlot(cov)
}
