
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

