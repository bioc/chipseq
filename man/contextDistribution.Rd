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
  
  \strong{DEPRECATED}: The user is likely better off using the
  GenomicFeatures package to explore the relationship between peaks and
  genomic annotations. See the vignette for an example.
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
    \code{\link[GenomicFeatures:regions]{transcripts}}.
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
\keyword{manip}

