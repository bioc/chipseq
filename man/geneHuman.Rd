\name{geneHuman}
\Rdversion{1.1}
\alias{geneHuman}
\docType{data}
\title{UCSC Gene Predictions for hg18}
\description{
A set of genomic coordinate data for gene transcripts, and related
annotations, for H. sapiens. Coordinates are relative to the hg18 build.
All coordinates are in nucleotides from the 5' end of the positive
\code{"+"} strand.In the below, \dQuote{gene} is synonymous with
\dQuote{spliced transcript}.
}
\usage{data(geneHuman)}
\format{
  A data frame with 56722 observations on the following 12
  variables.
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
      are coded as NA}.
    \item{\code{alignID}}{Unique identifier of each gene and RNA
      alignment pair, apparently redundant with \code{name}.}
  }
}
\details{
For genes coded on the negative strand the \code{txStart} is really the end,
and similarly for the coding regions.
}

\source{
This table was taken directly from the knownGene table in the UCSC
database for hg18, see \url{http://genome.ucsc.edu/cgi-bin/hgTables} and
Hsu F, Kent WJ, Clawson H, Kuhn RM, Diekhans M, Haussler D. The UCSC
Known Genes. Bioinformatics. 2006 May 1;22(9):1036-46.
} 
\examples{
data(geneHuman)
str(geneHuman)
}
\keyword{datasets}
