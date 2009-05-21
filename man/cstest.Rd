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
  Short Read Archive, GEO accession number GSM288351
  \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM288351}
}
\references{
  
  Chen X., Xu H., Yuan P., Fang F., Huss M., Vega V.B., Wong E., Orlov
  Y.L., Zhang W., Jiang J., Loh Y.H., Yeo H.C., Yeo Z.X., Narang V.,
  Govindarajan K.R., Leong B., Shahab A.S., Ruan Y., Bourque G., Sung
  W.K., Clarke N.D., Wei C.L., Ng H.H. (2008), \dQuote{Integration of
  External Signaling Pathways with the Core Transcriptional Network in
  Embryonic Stem Cells}. \emph{Cell}, 133:1106-1117.
  
}
\examples{
data(cstest)
names(cstest)
cstest$gfp
str(cstest$ctcf$chr10)
}
\keyword{datasets}
