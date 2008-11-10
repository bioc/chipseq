
## FIXME: should generalize to conditional plots if i is a vector

plotPeak <-
    function(peaks1, peaks2 = NULL, i = 1,
             xlab = "Position", ylab = "Coverage",
             ...) 
{
    pos1 <- seq(start(peaks1[i]), end(peaks1[i]))
    cov1 <- as.integer(peaks1[[i]])
    pos1 <- c(head(pos1, 1), pos1, tail(pos1, 1))
    cov1 <- c(0, cov1, 0)
    if (is.null(peaks2)) 
    {
        pos2 <- numeric(0)
        cov2 <- numeric(0)
    }
    else 
    {
        pos2 <- seq(start(peaks2[i]), end(peaks2[i]))
        cov2 <- -as.integer(peaks2[[i]])
        pos2 <- c(head(pos2, 1), pos2, tail(pos2, 1))
        cov2 <- c(0, cov2, 0)
    }
    xyplot(c(cov1, rev(cov2)) ~ c(pos1, rev(pos2)),
           ...,
           panel = function(...) {
               panel.polygon(..., col = "lightgrey")
               panel.abline(h = 0)
           },
           xlab = xlab, ylab = ylab)
}
