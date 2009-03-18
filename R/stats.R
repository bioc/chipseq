basesCovered <- function(x, shift = seq(0, 500, 5), seqLen = 35, verbose = FALSE)
{
    maxShift <- max(shift)
    rng <- range(unlist(x)) + c(-1, 1) * maxShift
    cov.pos <- coverage(extendReads(x, seqLen = seqLen, strand = "+"), rng[1], rng[2]) > 0
    cov.neg <- coverage(extendReads(x, seqLen = seqLen, strand = "-"), rng[1], rng[2]) > 0
    n <- diff(rng) + 1L
    ans <- shiftApply(shift, cov.pos, cov.neg, function(x, y) sum(x | y),
                      verbose = verbose)
    data.frame(mu = seqLen + shift, covered = ans / ans[1])
}


densityCorr <- function(x, shift = seq(0, 500, 5), ...)
{
    nShift <- length(shift)
    maxShift <- max(shift)
    rng <- range(unlist(x)) + c(-1, 1) * maxShift
    dotArgs <- list(...)
    if ("n" %in% names(dotArgs)) {
        n <- dotArgs[["n"]]
    } else {
        n <- min(65536, round(diff(rng) / 2))
    }
    d1 <- density(x$"+", from = rng[1], to = rng[2], n = n, ...)
    d2 <- density(x$"-", from = rng[1], to = rng[2] + maxShift, n = n, ...)
    d2Shifted <- approx(d2, xout = as.vector(outer(d1$x, shift, FUN = "+")))
    d2Y <- matrix(d2Shifted$y, ncol = nShift)
    data.frame(mu = shift, corr = as.vector(cor(d1$y, d2Y)))
}



## Jothi et al method
## For every tag i in the sense strand, the nearest tag k
## in the antisense strand, downstream of i, is identified. Let j be the
## tag in the sense strand immediately upstream of k. Note that i and j
## could be the same tag. The mean DNA fragment length F is given by

## 2/n \sum_i (d(i,j) + d(j,k)/2),

## where n is the number of sense tags for which there exists a
## k and a j tag, and d(i,k)<=500

jothi.estimate <- function(x, maxDist = 500L)
{
    i <- x$"+"
    message(length(i), " reads")
    pos <- sort(unique(x$"+"))
    neg <- sort(unique(x$"-"))

    whichK <- findInterval(i, neg) + 1L
    whichK[whichK > length(neg)] <- NA_integer_
    k <- neg[whichK]

    whichJ <- findInterval(k - 1L, pos)
    whichJ[whichJ < 1L] <- NA_integer_
    j <- pos[whichJ]

    keep <- !is.na(j) & !is.na(k) & abs(i-k) < maxDist
    i <- i[keep]
    j <- j[keep]
    k <- k[keep]

    2 * sum(abs(i-j) + 0.5 * abs(j-k)) / sum(keep)
}
