\name{readReads}
\alias{readReads}
\alias{as.list,AlignedRead-method}
\alias{coerce,AlignedRead,GenomeData-method}

\alias{}
\title{ A function to read in Aligned Short Reads }
\description{
  This is a helper function for reading in aligned reads with a number
  of parameters preset at values we have found useful for analyzing
  ChIP-seq data.
}
\usage{
readReads(srcdir, lane, ...,
          include = "chr[0-9]+$", type = "MAQMapShort",
          simplify = TRUE, minScore = 15)

as.list(x)
as(object, "GenomeData")

}
\arguments{
  \item{srcdir}{The source directory.  }
  \item{lane}{ The name of the file for each lane (logical subset). }
  \item{\dots}{ Additional parameters. }
  \item{include}{ A regular expression indicating which chromosomes to
     retain. }
  \item{type}{ The type of alignment used (MAQ, Bowtie etc). }
  \item{simplify}{ Logical indicating whether the result should be
    reduced to a simpler \code{"GenomeData"} object, which only retains
    the locations of the alignments. }
  \item{minScore}{ A minimum quality score cutoff (possibly MAQ specific). }
  \item{x, object}{ An object of class \code{"AlignedRead"}.
    \code{as.list} converts the object into into a list of
    per-chromosome alignment locations, which can be further coerced
    into a \code{"GenomeData"} object by the coercion method.  These
    functions are used by \code{readReads} when \code{simplify=TRUE}. } 
}
\details{
  This has mainly been used for MAQ alignments.  Our default parameters
  are to include only autosomal chromosomes (there seem to be problems
  with the others that will require details).  We reduce to one read per
  start location and strand.
}

\seealso{
  \code{\link[ShortRead:readAligned]{readAligned}},
  \code{\link[BSgenome:GenomeData-class]{GenomeData}}
}

\value{
  If \code{simplify=FALSE}, a \code{"AlignedRead"} object; otherwise, a
  \code{"GenomeData"} object.
}

\author{ D. Sarkar }
\examples{
##None yet
}
