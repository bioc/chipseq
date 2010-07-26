### FIXME: have this facet by sample after we have sample-level data structure
## This will display the distribution aggregated over all given chromosomes
islandDepthPlot <- function(x, maxDepth = 20L)
{
  tab <- table(unlist(viewMaxs(slice(x, lower = 1))))
  df <- data.frame(depth = as.numeric(names(tab)), count = as.numeric(tab))
  xyplot(log(count) ~ depth, df, 
         subset = (depth <= maxDepth), 
         pch = 16, type = c("p", "g"),
         panel = function(x, y, ...) {
           lambda <- 2 * exp(y[2]) / exp(y[1])
           null.est <- function(xx) {
             xx * log(lambda) - lambda - lgamma(xx + 1)
           }
           log.N.hat <- null.est(1) - y[1]
           panel.lines(1:10, -log.N.hat + null.est(1:10), col = "black")
           panel.xyplot(x, y, ...)
         })
}

peakCutoff <- function(cov, fdr.cutoff = 0.001, k = 2:20)
{
  ## an implementation of the idea in Robertson et al to assess
  ## sufficiency of sampling depth: choose minimum cutoff that gives
  ## an FDR < pre-specified value
  s <- slice(cov, lower = 1)
  y <- table(unlist(viewMaxs(s)))
  lambda <- 2 * y[2] / y[1]
  n <- exp(log(y[1]) - dpois(1, lambda, log = TRUE))
  exp.fd <- n * ppois(k-1, lambda, lower.tail = FALSE)
  obs.d <- integer(length(k))
  for (i in seq_along(k))
    {
      obs.d[i] <- sum(y[as.integer(names(y)) >= k[i]])
    }
  FDR <- ifelse(obs.d == 0, 0, exp.fd / obs.d)
  fdr.ok <- which(FDR < fdr.cutoff)
  if (length(fdr.ok) < 1)
    stop("No cutoff with low enough FDR found")
  fdr.chosen <- fdr.ok[1]
  k[fdr.chosen-1] + (FDR[fdr.chosen-1] - fdr.cutoff) /
    (FDR[fdr.chosen-1] - FDR[fdr.chosen])
}

## Not exported

## Useful summary functions for use with gdApply()
## e.g., as(gdApply(x, countSummary), "data.frame")

## x is a list at the lane->chromosome level, with components "+" and "-"

countSummary <- function(x) 
{
  npos <- length(x$"+")
  nneg <- length(x$"-")
  data.frame(n = npos + nneg, d = npos - nneg, r = npos / nneg)
}    

getSingletons <- function(x, ...)
{
  ## We retain length-400 islands (actually 2 adjoint length 200), but
  ## these are hopefully rare enough not to matter.
  g <- extendReads(x, ...)
  cov <- coverage(g, width = max(end(g) + 400L))
  s <- slice(cov, lower = 1)
  s <- s[viewMaxs(s) == 1]
  0.5 * (start(s) + end(s))
}
