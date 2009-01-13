\name{geneMouse}
\Rdversion{1.1}
\alias{geneMouse}
\docType{data}
\title{
Genomic coordinates for various gene related quantities in the mouse
mm9 build.
}
\description{
A set of genomic coordinate data for genes, and related quantites for 
M. musculus. Coordinates are relative to the mm9 build.  All coordinates
are in nucleotides from the 5' end of the positive \code{"+"} strand.
}
\usage{data(geneMouse)}
\format{
  A data frame with 49409 observations on the following 12 variables.
  \describe{
    \item{\code{name}}{The name of the gene.}
    \item{\code{chrom}}{The name of the chromosome the gene is located on.}
    \item{\code{strand}}{The strand the gene is coded on, \code{"+"}, or 
         \code{"-"}.}
    \item{\code{txStart}}{Transcription start site.}
    \item{\code{txEnd}}{Transcription stop site.}
    \item{\code{cdsStart}}{Start position of the coding sequence.}
    \item{\code{cdsEnd}}{End position of the coding sequence.}
    \item{\code{exonCount}}{The number of exons.}
    \item{\code{exonStarts}}{A comma separated list of the exon start positions.}
    \item{\code{exonEnds}}{A comma separated list of exon stop positions.}
    \item{\code{proteinID}}{An ID for the protein produced, missing values
      are coded as the empty string, \code{""}}.
    \item{\code{alignID}}{Not quite sure.}
  }
}
\details{
For genes coded on the negative strand the \code{txStart} is really the end,
and similarly for the coding regions.
}

\source{
This table was taken directly from the knownGene table in the UCSC database for mm9.
}
\examples{
data(geneMouse)
str(geneMouse)
}
\keyword{datasets}
