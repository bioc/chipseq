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
    logical flag.  If \code{TRUE}, only unique locations/reads will be
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
\keyword{manip}
\keyword{utilities} 

