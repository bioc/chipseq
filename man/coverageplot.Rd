\name{coverageplot}
\Rdversion{1.1}
\alias{coverageplot}
\title{
  Plot coverage on a small interval.
}
\description{
  A function that plots one or two coverage vectors over a relatively
  small interval in the genome.
}
\usage{
coverageplot(peaks1, peaks2 = NULL, i = 1,
             xlab = "Position", ylab = "Coverage",
             opposite = TRUE, ...)
}

\arguments{
  \item{peaks1, peaks2}{
    A set of peaks as described by ranges over a coverage vector.
  }
  \item{i}{
    Which peak to use.
  }
  \item{xlab, ylab}{
    Axis labels.
  }
  \item{opposite}{
    Logical specifying whether the two peaks should be plotted on
    opposite sides (appropriate for positive and negative strand peaks).
  }
    
  \item{\dots}{
    extra arguments.
  }

}

\author{
  Deepayan Sarkar
}

\keyword{ ~kwd1 }

