\name{readReads}
\alias{readReads}
\alias{as.list,AlignedRead-method}
\alias{coerce,AlignedRead,GenomeData-method}
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
}
\details{
  This has mainly been used for MAQ alignments.  Our default parameters
  are to include only autosomal chromosomes (there seem to be problems
  with the others that will require details).  We reduce to one read per
  start location and strand.
}

\section{Coercion}{

  When \code{simplify=TRUE} is specified, the return value is simplified
  to contain only alignemnt locations (and not associated quality
  information, etc.).  This simplification can also be done afterwards
  through coercion methods:
  \describe{
    \item{}{\code{as.list(x)}: 
      where \code{x} is an object of class \code{"AlignedRead"}
    }
    \item{}{\code{as(object, "GenomeData")}: 
      where \code{object} is an object of class \code{"AlignedRead"}
    }
  }
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
\dontrun{
## load reads mapped to chr10 in lane 2 from current working directory
readReads(".", "s_2_export.txt", include = "chr10")
## load all chromosomes in lane 1 from Bowtie output (20 quality cutoff)
readReads(".", "s_1_export.txt", type="Bowtie", minScore=20)
}
}
