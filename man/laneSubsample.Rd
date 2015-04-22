\name{laneSubsample}
\Rdversion{1.1}
\alias{laneSubsample}
\title{
  Subsample short read alignment locations
}
\description{
  Subsamples data from multiple lanes on a per-chromosome
  basis.
}
\usage{
laneSubsample(lane1, lane2, fudge = 0.05)
}
\arguments{
  \item{lane1, lane2}{ Two lanes of data, each of class
    \code{"GRanges"}.  }  
  \item{fudge}{ A numeric fudge factor.  For each chromosome, if the
    difference in the sizes relative to the size of the first dataset is
    less than \code{fudge}, no subsampling is done.
  }
}
\value{
  \code{laneSubsample} returns a list similar to its input, but with the
  larger dataset subsampled to be similar to the smaller one.
}
\author{
  D. Sarkar
}
\examples{
data(cstest)
## subsample to compare lanes
cstest.sub <- laneSubsample(cstest[[1]], cstest[[2]])
unlist(cstest.sub)
}
\keyword{manip}
\keyword{utilities} 

