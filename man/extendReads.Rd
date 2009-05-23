\name{extendReads}
\alias{extendReads}
\title{ A function to extend short reads. }
\description{
  Since the short read is typically represents one end of a longer
fragment there are situations where extending it to the approximate length
of the fragment can be useful.
}
\usage{
extendReads(reads, seqLen = 200, strand = c("+", "-"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reads}{ Either an \code{AlignedReads} object or a \code{list} of
       \code{AlignedReads} objects (or a \code{list} with aligned reads
      for each strand.) }
  \item{seqLen}{ The desired length of the final sequence, assumed to be
   the same for all reads. }
  \item{strand}{ Which strand \code{+} or \code{-} the read is aligned to. }
}
\details{
Read locations are presumed to be the 5' end (relative to the + strand of
the chromosome).  Thus reads on the plus strand are simply extended. Those
that align to the minus strand, we must subtract the read length,
then grow the read towards the 5' end of the + strand (3' end of the 
minus strand).
}
\value{
An \code{IRanges} object with the new ranges, or a list of IRanges objects,
depending on the input.
}
\author{ R. Gentleman}
\examples{
data(cstest)
extRanges <- gdapply(cstest, extendReads, seqLen = 200)
}
