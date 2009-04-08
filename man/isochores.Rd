\name{isochores.mm8}
\Rdversion{1.1}
\alias{isochores.mm8}
\docType{data}
\title{
  Isochore boundaries for Mus musculus (build mm9).
}
\description{
  Isochore boundaries for Mus musculus (build mm9).  Isochores are large
  segments of the genome such that within-segment variability in GC
  content is substantially lower than between-segment variability.
  These isochores are computationally predicted by IsoFinder (see below).
}
\usage{
data(isochores.mm8)
}
\format{
  A data frame with 32894 observations on the following 4 variables.
  \describe{
    \item{\code{Begin}}{ isochore starts. }
    \item{\code{End}}{ isochore ends. }
    \item{\code{GC}}{ GC content in isochore. }
    \item{\code{chromosome}}{ chromosome identifier. }
  }
}
\source{
  \url{http://bioinfo2.ugr.es/isochores/}
}
\keyword{datasets}
