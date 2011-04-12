\name{estimate.mean.fraglen}
\Rdversion{1.1}

\alias{estimate.mean.fraglen}
\alias{estimate.mean.fraglen,AlignedRead-method}
\alias{estimate.mean.fraglen,GenomeData-method}
\alias{estimate.mean.fraglen,GRanges-method}

\alias{basesCovered}
\alias{densityCorr}
\alias{sparse.density}

\title{
  Estimate summaries of the distribution of fragment lengths in a
  short-read experiment.  The methods are designed for ChIP-Seq
  experiments and may not work well in data without peaks.
}
\description{
  \code{estimate.mean.fraglen} implements three methods for estimating
  mean fragment length.  The other functions are related helper
  functions implementing various methods, but may be useful by
  themselves for diagnostic purposes.  Many of these operations are
  potentially slow.

  \code{sparse.density} is intended to be similar to
  \code{\link{density}}, but returns the results in a run-length encoded
  form.  This is useful when long stretches of the range of the data
  have zero density.
}
\usage{
estimate.mean.fraglen(x, method = c("SISSR", "coverage", "correlation"),
                      \dots)

basesCovered(x, shift = seq(5, 300, 5), seqLen = 35, verbose = FALSE)

densityCorr(x, shift = seq(0, 500, 5), center = FALSE, width = 50, \dots)

sparse.density(x, width = 50, kernel = "epanechnikov",
               experimental = TRUE, from, to)
}
\arguments{
  \item{x}{
    For \code{estimate.mean.fraglen}, typically an
    \link[ShortRead:AlignedRead-class]{AlignedRead} or a
    \link[GenomicRanges:GRanges-class]{GRanges} object. Also supported
    but deprecated, as they do not have formal strand information:
    \link[IRanges:RangedData-class]{RangedData} (with a "strand"
    column), or a list-like object with elements \code{"+"} and \code{"-"}
    representing locations of reads aligned to positive and negative
    strands (the values should be integers denoting the location where
    the first sequenced base matched.) Supported (but again, deprecated)
    list types include: \link[IRanges:RangesList-class]{RangesList},
    \link[IRanges:IntegerList-class]{IntegerList} or an ordinary R
    list. Also a \link[BSgenome:GenomeData-class]{GenomeData} where the
    per-chromosome elements are one of the above list types.

    For \code{basesCovered}, and \code{densityCorr}, a list with elements
    \code{"+"} and \code{"-"} representing locations of reads aligned to
    positive and negative strands (the values should be integers denoting
    the location where the first sequenced base matched.)

    For \code{sparse.density}, a numeric or integer vector for which
    density is to be computed.  
  }
  
  \item{method}{ Character string giving method to be used.
    \code{method = "SISSR"} implements the method described in Jothi et
    al (see References below).  \code{method = "correlation"} implements
    the method described in Kharchenko et al (see References below),
    where the idea is to compute the density of tag start positions
    separately for each strand, and then determine the amount of shift
    that maximizes the correlation between these two densities.
    \code{method = "coverage"} computes the optimal shift for which the
    number of bases covered by any read is minimized.
  }
  
  \item{shift}{ Integer vector giving amount of shifts to be tried when
    optimizing.  The current algorithm simply evaluates all supplied
    values and reports the one giving minimum coverage or maximum
    correlation. }
  
  \item{seqLen}{ For the \code{"coverage"} method, the assumed length of
    each read for computing the coverage.  Typically the read
    length. This is added to the shift estimated by \code{"coverage"}
    and \code{"correlation"} to come up with the actual fragment
    length.}

  \item{verbose}{ Logical specifying whether progress information should
    be printed during execution. }

  \item{center}{ For the \code{"correlation"} method, whether the
    calculations should incorporate centering by the mean density.  The
    default is not to do so; as the density is zero over most of the
    genome, this slightly improves efficiency at negligible loss in
    accuracy. }

  \item{width}{ half-bandwidth used in the computation.  This needs to
    be specified as an integer, data-driven rules are not supported. }
  
  \item{kernel}{ A character string giving the density kernel. }

  \item{experimental}{ logical.  If \code{TRUE} }

  \item{from, to}{ specifies range over which the density is to be
    computed. }

  \item{\dots}{ Extra arguments, passed on as appropriate to other
    functions. }
}
\details{
  These functions are typically used in conjunction with
  \code{\link[BSgenome:gdapply]{gdapply}}.
  
  For the correlation method, the range over which densities are
  computed only cover the range of reads; that is, the beginning and end
  of chromosomes are excluded.
}
\value{
  \code{estimate.mean.fraglen} gives an estimate of the mean fragment
  length.  
  
  \code{basesCovered} and \code{densityCorr} give a vector of the
  corresponding objective function evaluated at the supplied values of
  \code{shift}.
  
  \code{sparse.density} returns an object of class \code{"Rle"}.
}
\references{
  R. Jothi, S. Cuddapah, A. Barski, K. Cui, and K. Zhao. Genome-wide
  identification of in vivo protein-DNA binding sites from ChIP-Seq
  data.  \emph{Nucleic Acids Research}, 36:5221--31, 2008.
  
  P. V. Kharchenko, M. Y. Tolstorukov, and P. J. Park. Design and
  analysis of ChIP experiments for DNA-binding proteins. \emph{Nature
  Biotechnology}, 26:1351--1359, 2008.
}
\author{
  Deepayan Sarkar, Michael Lawrence
}
\examples{
data(cstest)
estimate.mean.fraglen(cstest[["ctcf"]], method = "coverage")
}
\seealso{
  \code{\link[BSgenome:gdapply]{gdapply}}
}
\keyword{univar}
