\name{contextDistribution}
\Rdversion{1.1}
\alias{contextDistribution}
\title{
  Tabulate peak locations according to genomic context
}
\description{
  Given two sets of intervals defined on a genome, tabulates overlap of
  one set with the other.  The first set typically represents
  \dQuote{peak} locations, and the second represents types of genomic
  regions such as promoters, downstream regions, genes, etc. 
}
\usage{
contextDistribution(peaks, gregions, chroms, ...)
}
\arguments{
  \item{peaks}{
    A data frame with one row for each \dQuote{peak}; the location of
    peaks must be defined by the columns \code{chromosome},
    \code{start}, and \code{end}.  Columns \code{up} and \code{down}, if
    present, must be logical, and should indicate peaks that were down
    or upregulated by some definition.  If present, the result will
    include tabulations for the up and down subsets thus defined.
  }
  \item{gregions}{
    Locations of genomic regions of interest.  Currently, this must be
    of the form produced by the function
    \code{\link[GenomicFeatures]{transcripts}}.
  }
  \item{chroms}{
    Which chromosomes to use.  By default, all are used.
  }
  \item{\dots}{
    Further arguments, currently ignored.
  }
}
\value{
  A data frame with overlap counts.
}
\author{
  Deepayan Sarkar
}
\examples{
data(cstest)
library(BSgenome.Mmusculus.UCSC.mm9)
mouse.chromlens <- seqlengths(Mmusculus)
## extend reads, generate peak summary
extRanges <- gdapply(cstest, extendReads, seqLen = 200)
peakSummary <-
    diffPeakSummary(extRanges$gfp, extRanges$ctcf,
                    chrom.lens = mouse.chromlens, lower = 10)
## generate transcripts using GenomicFeatures package
library(GenomicFeatures)
data(geneMouse)
gregions <- transcripts(genes = geneMouse, proximal = 2000)
## finally, calculate context distribution for chr10
contextDistribution(peakSummary, gregions, "chr10")
}
\keyword{manip}

