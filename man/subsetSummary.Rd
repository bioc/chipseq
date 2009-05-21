\name{subsetSummary}
\Rdversion{1.1}
\alias{subsetSummary}
\title{
  Compute summaries for cumulative subsets of a short-read data set. 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
subsetSummary(x, chr, nstep, props = seq(0.1, 1, 0.1),
              chromlens, fg.cutoff = 6, seqLen = 200,
              resample = TRUE, islands = TRUE,
              verbose = getOption("verbose"))
}

\arguments{
  \item{x}{}
  \item{chr}{}
  \item{nstep}{}
  \item{props}{}
  \item{chromlens}{}
  \item{fg.cutoff}{}
  \item{seqLen}{}
  \item{resample}{}
  \item{islands}{}
  \item{verbose}{}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A data frame with various per-subset summaries.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Deepayan Sarkar
}
\note{
%%  ~~further notes~~
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\keyword{ ~kwd1 }
