\name{chipseqFilter}
\alias{chipseqFilter}
\title{
Filtering ChIP-seq reads
}
\description{
  Convenience for creating an
  \code{\link[ShortRead:SRFilter-class]{SRFilter}} object appropriate
  for ChIP-seq data. Typically, the result is passed
  to \code{\link[ShortRead]{readAligned}} when loading reads.
}
\usage{
chipseqFilter(exclude = "[_MXY]", uniqueness = c("location", "sequence", "location*sequence", "none"), hasStrand = TRUE)
}
\arguments{
  \item{exclude}{
    A regular expression for excluding chromosomes by name. Just
    like the parameter to \code{\link[BSgenome]{bsapply}}. 
  }
  \item{uniqueness}{
    The criteria used to determine whether a read is unique. A read may
    be unique if it maps to a unique \code{location}, has a
    unique \code{sequence} or both. Specifying \code{none} avoids this
    test entirely.
  }
  \item{hasStrand}{
    Whether to require that the read is mapped to a strand, which usually
    translates to whether the read was mapped at all.
  }
}
\value{
An \code{SRFilter} object
}
\author{
  M. Lawrence
}

\examples{
sp <- SolexaPath(system.file("extdata", package="ShortRead"))

filter <- chipseqFilter()
aln <- readAligned(sp, "s_2_export.txt", filter=filter)
## allow mapping to the same location (but only if sequence is different)
filter <- chipseqFilter(uniqueness = "sequence")
aln <- readAligned(sp, "s_2_export.txt", filter=filter)
## allow sex chromosomes
filter <- chipseqFilter(exclude = "[M_]")
aln <- readAligned(sp, "s_2_export.txt", filter=filter)
}