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
\examples{
cov <- Rle(c(1:10, seq(10, 1, -2), seq(1,5,2), 4:1), rep(1:2, 11))
peaks <- slice(cov, 3)
peaks.cov <- Views(cov, peaks)
peaks.cov.rev <- rev(peaks.cov)
coverageplot(peaks.cov, peaks.cov.rev, ylab = "Example")
}
\author{
  Deepayan Sarkar
}

\keyword{hplot}

