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
    of the form produced by the function \code{\link{genomic_regions}}.
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

  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
}
\author{
  Deepayan Sarkar
}
\keyword{manip}

