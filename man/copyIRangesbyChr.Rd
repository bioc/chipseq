\name{copyIRangesbyChr}
\Rdversion{1.1}
\alias{copyIRanges}
\alias{copyIRangesbyChr}
\title{
  Associate ranges to coverage. 
}
\description{

  Associate a set of ranges, typically derived using an independent
  computation, to a coverage as produced by \code{coverage}.  This then
  allows one to compute various summaries such as maximum coverage in
  each range.  \code{copyIRangesbyChr} does this over lists of ranges
  and coverage objects.

}
\usage{
copyIRanges(IR1, newX)
copyIRangesbyChr(IR1, newX)
}

\arguments{
  \item{IR1}{ The set of ranges (an \code{"IRanges"} object) or a list of
    such objects (usually one for each chromosome of interest). }
  \item{newX}{ An \code{"Rle"} object, usually the result of
    \code{link[IRanges:coverage]{coverage}}, or a list of such objects.  }
}
\value{
  A \code{"View"} object, or a list of such objects.
}
\author{
  Deepayan Sarkar
}
\keyword{utilities}

