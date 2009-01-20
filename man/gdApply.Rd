\name{gdApply}
\Rdversion{1.1}
\alias{gdApply}
\title{
  Applies a function to elements of a GenomeData
}
\usage{
gdApply(X, FUN, ...)
}
\arguments{
  \item{X}{
    An object of class \code{"GenomeData"} or \code{"GenomeDataList"}
  }
  \item{FUN}{
    A function to be applied to each chromosome-level sub-element of
    \code{X}.
  }
  \item{\dots}{
    Further arguments; passed to \code{FUN}
  }
}
\value{
  Typically an object of the same class as \code{X}.
}
\author{
  Deepayan Sarkar
}
\keyword{manip}

