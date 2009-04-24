\name{combineLanes}
\Rdversion{1.1}
\alias{combineLanes}
\title{
  Combine lanes of sequence fragment locations
}
\description{
  Combines data from multiple lanes on a per-chromosome basis.
}
\usage{
combineLanes(x, chromList, keep.unique = FALSE)
}
\arguments{
  \item{x}{
    Typically a code{"GenomeDataList"} object representing multiple
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
}
\value{
  An object of class \code{"GenomeData"}.
}
\author{
  D. Sarkar
}
\keyword{manip}

