\name{readReads}
\alias{readReads}
\title{ A function to read in Aligned Short Reads }
\description{
 This is a helper function for reading in aligned reads with a number of 
parameters preset at values we have found useful for analyzing ChIP-seq data. 
}
\usage{
readReads(srcdir, lane, ..., include = "chr[0-9]+$", type = "MAQMapShort", simplify = TRUE, minScore = 15)
}
\arguments{
  \item{srcdir}{The source directory.  }
  \item{lane}{ The name of the file for each lane (logical subset). }
  \item{\dots}{ Additional parameters. }
  \item{include}{ A regular expression indicating which chromosomes to
     retain. }
  \item{type}{ The type of alignment used (MAQ, Bowtie etc). }
  \item{simplify}{ Indicate whether to return a list (\code{TRUE}). }
  \item{minScore}{ The minimum quality score to include (possibly MAQ specific. }
}
\details{
  This has mainly been used for MAQ data.  Our default parameters are
to include only autosomal chromosomes (there seem to be problems with
the others that will require details).  We reduce to one read per start
location and strand.

}
\seealso{\code{\link[ShortRead:readAligned]{readAligned}}}
\value{
An instance of \code{AlignedRead}.
}
\author{ D. Sarkar }
\examples{
##None yet
}
