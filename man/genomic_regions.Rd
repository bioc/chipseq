\name{genomic_regions}
\Rdversion{1.1}
\alias{genomic_regions}
\alias{genomic_exons}
\title{
Functions that compute genomic regions of interest. 
}
\description{
Functions that compute genomic regions of interest such as promotor, 
upstream regions etc, from the genomic locations provided in data like
\code{geneMouse}.
}
\usage{
genomic_regions(genes, proximal = 500, distal = 10000, offset = 1)
genomic_exons(genes)
}
\arguments{
  \item{genes}{
     A data.frame like that provided by \code{geneMouse}.
}
  \item{proximal}{
  The number of bases on either side of TSS and 3'-end for
  the promoter and end region, respectively.
}
  \item{distal}{
The number of bases on either side for upstream/downstream,
i.e. enhancer/silencer regions.
}
  \item{offset}{ An offset to allow for a gap (positive or negative) 
 between the region defined by \code{proximal} and that defined by 
\code{distal}.
}
\details{
Fairly simply addtions/subtractions are made.
}
\value{
For \code{genomic_regions} a data.frame with all components computed.
For \code{genomic_exons} a data.frame with one row per exon.
}
\author{
M. Lawrence.
}

\examples{
##none yet
}
