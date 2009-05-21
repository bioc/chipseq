\name{cstest}
\Rdversion{1.1}
\alias{cstest}
\docType{data}
\title{
  A test Chip-Seq dataset
}
\description{
  A small subset of a Chip-Seq dataset downloaded from the Short-Read
  Archive.
}
\usage{data(cstest)}
\format{

  The dataset is on object of class \code{GenomeDataList} with data from
  three chromosomes in two lanes representing CTCF and GFP pull-down in
  mouse.

  The per-chromosome data is represented as a list of positive and
  negative strand alignment locations.  The recorded locations represent
  the aligned position at the first cycle.
}
\source{
  Short Read Archive, \url{FIXME}
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(cstest)
names(cstest)
cstest$gfp
str(cstest$ctcf$chr10)
}
\keyword{datasets}
