\name{combineLanes}
\Rdversion{1.1}
\alias{combineLanes}
\alias{laneSubsample}
\title{
  Combine or subsample short read alignment locations
}
\description{
  Combines or subsamples data from multiple lanes on a per-chromosome
  basis.
  \strong{THIS FUNCTION IS DEPRECATED}: no longer needed
  given current infrastructure. With a
  \code{\link[GenomicRanges:GRanges-class]{GRanges}}, simply combine
  with \code{c} and use \code{unique} to mimic
  setting \code{keep.unique} to \code{TRUE}.
}
\usage{
combineLanes(x, chromList, keep.unique = FALSE)
laneSubsample(lane1, lane2, fudge = 0.05)
}
\arguments{
  \item{x}{
    Typically a \code{"GenomeDataList"} object representing multiple
    lanes of aligned locations or ranges.  The result will combine the
    locations across lanes on a per-chromosome basis.
  }
  \item{chromList}{
    Character vector specifying Which chromosomes to combine.  Defaults
    to all chromosomes in the first lane.
  }
  \item{keep.unique}{
    logical flag.  If \code{TRUE}, only unique locations will be
    retained.
  }
  \item{lane1, lane2}{ Two lanes of data, each of class
    \code{"GenomeData"}.  }
  
  \item{fudge}{ A numeric fudge factor.  For each chromosome, if the
    difference in the sizes relative to the size of the first dataset is
    less than \code{fudge}, no subsampling is done.
  }
}
\value{
  \code{combineLanes} returns an object of class \code{"GenomeData"}.

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
## two lanes of chr10 become one
##combineLanes(cstest, "chr10") # DEPRECATED
unlist(cstest.sub) # instead (for all chromosomes)
}
\keyword{manip}
\keyword{utilities} 

