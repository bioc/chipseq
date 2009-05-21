\name{subsetSummary}
\Rdversion{1.1}
\alias{subsetSummary}
\title{
  Compute summaries for cumulative subsets of a short-read data set. 
}
\description{
  Divides a short-read dataset into several subsets, and computes
  various summaries cumulatively.  The goal is to study the
  characteristics of the data as a function of sample size.
}
\usage{
subsetSummary(x, chr, nstep, props = seq(0.1, 1, 0.1),
              chromlens, fg.cutoff = 6, seqLen = 200,
              resample = TRUE, islands = TRUE,
              verbose = getOption("verbose"))
}

\arguments{
  \item{x}{ A \code{"GenomeDataList"} object representing alignment
    locations at the sample level. }
  \item{chr}{ The chromosome for which the summaries are to be obtained.
    Must specify a valid element of \code{x} }
  \item{nstep}{ The number of maps in each increment for the full
    dataset (not per-chromosome).  This will be translated to a
    per-chromosome number proportionally.  }
  \item{props}{ Alternatively, an increasing sequence of proportions
    determining the size of each subset.  Overrides \code{nstep}. }
  \item{chromlens}{ A named vector of per-chromosome lengths, typically
    the result of \code{\link[BSgenome:seqlengths]{seqlengths}}. }
  \item{fg.cutoff}{ The coverage depth above which a region would be
    considered foreground. }
  \item{seqLen}{ The number of bases to which to extend each read before
  computing coverage. }
  \item{resample}{ Logical; whether to randomly reorder the reads before
    dividing them up into subsets.  Useful to remove potential order
    effects (for example, if data from two lanes were combined to
    produce \code{x}). }
  \item{islands}{ Logical.  If \code{TRUE}, the whole island would be
    considered foreground if the maximum depth equals or exceeds
    \code{fg.cutoff}.  If \code{FALSE}, only the region above the cutoff
    would be considered foreground. }
  \item{verbose}{ logical controlling whether progress information will
    be shown during computation (which is potentially long-running).  }
}
\value{
  A data frame with various per-subset summaries.
}
\author{
  Deepayan Sarkar
}
\keyword{ ~kwd1 }
