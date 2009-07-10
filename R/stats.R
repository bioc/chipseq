## The first set of functions provide a version of density for
## 'sparse' data (large stretches of empty space where density is 0;
## so large that standard density() would not be able to compute it.
## FIXME: this is definitely not the right place to put this, but I'm
## not sure what is (IRanges has the underlying classes, but that does
## not necessarily make it the right choice).

## The use-case of interest is Kharchencko et al's density-correlation
## method for estimating average fragmet length

## FIXME: the "efficient" summary functions are not as important now
## because IRanges has improved.  At some point, should consider
## getting rid of them.


dKernel <-
    function(width = 50,
             kernel = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine"))
{
    kernel <- match.arg(kernel) 
    kords <- seq(-width, width)
    switch(kernel,
           gaussian = dnorm(kords, sd = width / 3),
           rectangular = {
               a <- width
               ifelse(abs(kords) < a, 0.5/a, 0)
           }, triangular = {
               a <- width
               ax <- abs(kords)
               ifelse(ax < a, (1 - ax/a)/a, 0)
           }, epanechnikov = {
               a <- width
               ax <- abs(kords)
               ifelse(ax < a, 3/4 * (1 - (ax/a)^2)/a, 0)
           }, biweight = {
               a <- width
               ax <- abs(kords)
               ifelse(ax < a, 15/16 * (1 - (ax/a)^2)^2/a, 0)
           }, cosine = {
               a <- width
               ifelse(abs(kords) < a, (1 + cos(pi * kords/a))/(2 * a), 0)
           }, optcosine = {
               a <- width
               ifelse(abs(kords) < a, pi/4 * cos(pi * kords/(2 * a))/a, 0)
           })
}


## pad a 0 in the beginning to fill in stretch of 0's before

doDensity1 <- function(x, kernel, width) ## calls R's density() 
{
    n <- length(x)
    from <- x[1] - width
    to <- x[n] + width
    c(0, n * density(x, kernel = kernel, width = 2 * width, from = from, to = to, n = to - from + 1)$y)
}

doDensity2 <- function(x, dk, width) ## explicit looping (but dk not recomputed)
{
    n <- length(x)
    if (n == 1)
    {
        c(0, dk)
    }
    else
    {
        x <- x - (x[1] - width - 1)
        from <- x[1] - width
        to <- x[n] + width
        len <- to - from + 1
        ans <- numeric(len)
        for (xi in x)
        {
            si <- seq(xi - width, xi + width)
            ans[si] <- ans[si] + dk
        }
        c(0, ans)
    }
}

doDensity3 <- function(x, dk, width) ## loop in C
{
    .Call("do_naive_density",
          as.integer(x), as.double(dk), as.integer(width),
          PACKAGE = "chipseq")
}


## experimental=FALSE calls R's density() within each "island".
## experimental=TRUE pre-computes the kernel function for the given
## width and contructs the density by looping in R.  This gives a ~60x
## speedup (most islands have few reads, and the loops are short
## there).  This (the tapply part) can easily be ported to C to be
## even faster.

sparse.density <- function(x, width = 50, kernel = "epanechnikov", experimental = TRUE,
                           from = start(rix)[1] - 10L,
                           to = end(rix)[length(rix)] + 10L)
{
    if (!is.numeric(x))
        stop("'x' must be an integer vector")
    if (!is.integer(x))
        x <- as.integer(x)
    x <- sort(x)
    ix <- IRanges(x-width, x+width)
    rix <- sort(reduce(ix))
    ## we will calculate density on a range containing the data.  If
    ## necessary, we will subset later (FIXME: TODO).
    from0 <- min(from, start(rix)[1] - 1L)
    to0 <- max(to, end(rix)[length(rix)] + 1L)
    if (from > from0 || to < to0) stop("[from, to] smaller than support not implemented yet (but easy to add)")
    ## ox <- overlap(rix, x, multiple = FALSE)
    ox <- findInterval(x, start(rix)) # equivalent, but a little faster
    island.densities <- 
        if (!experimental)
        {
            ## tapply(x, ox, doDensity1, kernel = kernel, width = width, simplify = FALSE) # really slow
            dk <- dKernel(width = width, kernel = kernel)
            tapply(x, ox, doDensity2, dk = dk, width = width, simplify = FALSE)
        }
        else 
        {
            dk <- dKernel(width = width, kernel = kernel)
            tapply(x, ox, doDensity3, dk = dk, width = width, simplify = FALSE)
        }
    ## result will be a "Rle" object.  Need to compute @values and
    ## @lengths.  We have the pieces, each of which will have
    ## values=island.densities[[i]][-1] and lengths=rep(1,width(rix)).
    ## In between, we have stretches of 0.
    zero.lengths <- c(start(rix), to0) - c(from0, end(rix)) - 1L
    all.values <- c(unlist(island.densities, recursive = FALSE, use.names = FALSE), 0)
    all.lengths <- rep(1L, length(all.values))
    all.lengths[cumsum(c(1L, sapply(island.densities, length)))] <- zero.lengths
    ## ans <- Rle(all.values, all.lengths) # has unnecessary fancy checks
    ans <- new("Rle", values = all.values, lengths = as.integer(all.lengths))
    ans
}


basesCovered <- function(x, shift = seq(5, 300, 5), seqLen = 35, verbose = FALSE)
{
    if (!is.list(x))
        stop("'x' must be a list object")
    if (!all(c("+", "-") %in% names(x)))
        stop("x must have named elements '+' and '-'")
    maxShift <- max(shift)
    rng <- range(unlist(x)) + c(-1, 1) * maxShift
    cov.pos <- coverage(extendReads(x, seqLen = seqLen, strand = "+"), shift = 1-rng[1], width = 1+diff(rng)) > 0
    cov.neg <- coverage(extendReads(x, seqLen = seqLen, strand = "-"), shift = 1-rng[1], width = 1+diff(rng)) > 0
    n <- diff(rng) + 1L
    ## ans <- shiftApply(shift, cov.pos, cov.neg, function(x, y) sum(x | y), verbose = verbose)
    ans <- shiftApply(shift, cov.pos, cov.neg, RleSumAny, verbose = verbose)
    data.frame(mu = seqLen + shift, covered = ans / ans[1])
}


## An efficient version of sum(e1 | e2), where e1 and e2 are Rle objects

RleSumAny <- function (e1, e2)
{
    len <- length(e1)
    stopifnot(len == length(e2))
    x1 <- runValue(e1); n1 <- runLength(e1); s1 <- cumsum(n1)
    x2 <- runValue(e2); n2 <- runLength(e2); s2 <- cumsum(n2)
    .Call("rle_sum_any",
          as.integer(x1), as.integer(n1), as.integer(s1),
          as.integer(x2), as.integer(n2), as.integer(s2),
          as.integer(len),
          PACKAGE = "chipseq")
}


## correlation from Rle objects.

similarity.corr <- function(pos, neg, center = FALSE)
{
    ## pos, neg are "Rle" objects
    if (center)
    {
        pos <- pos - mean(pos)
        neg <- neg - mean(neg)
    }
    RleSumProd(pos, neg) / sqrt(RleSumProd(pos, pos) * RleSumProd(neg, neg))
}


## this needs chromosome lengths, and includes all the 0-s on either
## side.  Not exported, as 

correlationProfile <-
    function(x, chrom,
             shift = seq(5, 300, 5),
             chrom.lengths,
             center = FALSE,
             ...)
{
    dl <- lapply(x[[chrom]], sparse.density,
                 from = 1L, to = chrom.lengths[chrom], ...)
    if (center) dl <- lapply(dl, function(x) { x - mean(x) })
    len <- length(dl[[1]])
    wid <- len - max(shift)
    cl <-
        sapply(shift,
               function(s) {
                   sumxy <- 
                       with(dl,
                            RleSumProd(subseq(`+`, start = 1L, width = wid),
                                       subseq(`-`, start = 1L + s, width = wid)))
                   sumxy
               })
    data.frame(mu = shift, corr = cl / with(dl( sqrt( RleSumProd(`+`, `+`) * RleSumProd(`-`, `-`)  ) )))
}


### really really slow
## RleSumProd <- function (x1, n1, x2, n2)
## {
##     ok <- e1 != 0 & e2 != 0
##     sum(e1[ok] * e2[ok])
## }

## An efficient version of sum(e1 * e2), where e1 and e2 are Rle objects

RleSumProd <- function (e1, e2)
{
    len <- length(e1)
    stopifnot(len == length(e2))
    x1 <- runValue(e1); n1 <- runLength(e1); s1 <- cumsum(n1)
    x2 <- runValue(e2); n2 <- runLength(e2); s2 <- cumsum(n2)
    ## rle_sum_prod_prototype(x1, n1, s1, x2, n2, s2, len)
    .Call("rle_sum_prod",
          as.double(x1), as.integer(n1), as.integer(s1),
          as.double(x2), as.integer(n2), as.integer(s2),
          as.integer(len),
          PACKAGE = "chipseq")
}


## a prototype, should go away.

rle_sum_prod_prototype <- function (x1, n1, s1, x2, n2, s2, len)
{
    i1 <- i2 <- k <- 1 #i: RLE pointer, k: virtual pointer of underlying seq
    ans <- 0
    while (k <= len)
    {
        if (x1[i1] == 0 || x2[i2] == 0) # both 0, no contribution
        {
            i1 <- i1 + 1
            i2 <- i2 + 1
            ## move lagging pointer ahead to location of the one ahead
            while (s1[i1-1] < s2[i2-1]) i1 <- i1 + 1
            while (s1[i1-1] > s2[i2-1]) i2 <- i2 + 1
            k <- 1 + max(s1[i1-1], s2[i2-1])
        }
        else # at this point, moving the lagging one forward by one must reach or jump over the other
        {
            next.k <- 1 + min(s1[i1], s2[i2])
            ans <- ans + (next.k - k) * x1[i1] * x2[i2]
            if (s1[i1] == next.k - 1) i1 <- i1 + 1
            if (s2[i2] == next.k - 1) i2 <- i2 + 1
            k <- next.k
        }
    }
    ans
}


## this version only uses the range of the data

densityCorr <- function(x, shift = seq(0, 500, 5), center = FALSE, width = 50, ...)
{
    if (!is.list(x))
        stop("'x' must be a list object")
    if (!all(c("+", "-") %in% names(x)))
        stop("x must have named elements '+' and '-'")
    maxShift <- max(shift)
    rng <- range(unlist(x)) + c(-1, 1) * (maxShift +  width)
    dl <- lapply(x, sparse.density, from = rng[1], to = rng[2], ...)
    if (center) dl <- lapply(dl, function(x) { x - mean(x) })
    len <- length(dl[[1]])
    wid <- len - max(shift)
    ## cl <- shiftApply(shift[1:10], dl$"+", dl$"-", FUN = similarity.corr, simplify = TRUE)
    ## cl <- shiftApply(shift[1:10], dl$"+", dl$"-", FUN = function(x, y) sum(x * y), simplify = TRUE)
    cl <- shiftApply(shift, dl$"+", dl$"-", FUN = RleSumProd, simplify = TRUE)
    ##     cl <-
    ##         sapply(shift,
    ##                function(s) {
    ##                    sumxy <- 
    ##                        with(dl,
    ##                             RleSumProd(subseq(`+`, start = 1L, width = wid),
    ##                                        subseq(`-`, start = 1L + s, width = wid)))
    ##                    sumxy
    ##                })
    data.frame(mu = shift, corr = cl / with(dl, sqrt( RleSumProd(`+`, `+`) * RleSumProd(`-`, `-`)  ) ))
}


## densityCorr <- function(x, shift = seq(0, 500, 5), ...)
## {
##     nShift <- length(shift)
##     maxShift <- max(shift)
##     rng <- range(unlist(x)) + c(-1, 1) * maxShift
##     dotArgs <- list(...)
##     if ("n" %in% names(dotArgs)) {
##         n <- dotArgs[["n"]]
##     } else {
##         n <- min(65536, round(diff(rng) / 2))
##     }
##     d1 <- density(x$"+", from = rng[1], to = rng[2], n = n, ...)
##     d2 <- density(x$"-", from = rng[1], to = rng[2] + maxShift, n = n, ...)
##     d2Shifted <- approx(d2, xout = as.vector(outer(d1$x, shift, FUN = "+")))
##     d2Y <- matrix(d2Shifted$y, ncol = nShift)
##     data.frame(mu = shift, corr = as.vector(cor(d1$y, d2Y)))
## }



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

## point estimates for the other methods

removeIsolated <- function(u, min.close = 500L) 
{
    u <- sort(u)
    du <- diff(u)
    min.dist <- pmin(c(Inf, du), c(du, Inf))
    u[min.dist <= min.close]
}


coverage.estimate <- function(x, maxDist = 500L, ...)
{
    if (!is.null(maxDist))
        x <- lapply(x, removeIsolated)
    d <- basesCovered(x, ...)
    with(d, mu[which.min(covered)])
}

correlation.estimate <- function(x, maxDist = 500L, ...)
{
    if (!is.null(maxDist))
        x <- lapply(x, removeIsolated)
    d <- densityCorr(x, ...)
    with(d, mu[which.max(corr)])
}

setGeneric("estimate.mean.fraglen", signature = "x",
           function(x, method = c("SISSR", "coverage", "correlation"), ...)
               standardGeneric("estimate.mean.fraglen"))

setMethod("estimate.mean.fraglen", "list",
          function(x, method = c("SISSR", "coverage", "correlation"), ...)
          {
              if (!all(c("+", "-") %in% names(x)))
                  stop("x must have named elements '+' and '-'")
              method <- match.arg(method)
              switch(method,
                     SISSR = jothi.estimate(x, ...),
                     coverage = coverage.estimate(x, ...),
                     correlation = correlation.estimate(x, ...))
          })

setMethod("estimate.mean.fraglen", "GenomeData",
          function(x, method = c("SISSR", "coverage", "correlation"), ...)
              unlist(lapply(x, estimate.mean.fraglen, ...)))

setMethod("estimate.mean.fraglen", "RangedData",
          function(x, method = c("SISSR", "coverage", "correlation"), ...) {
              unlist(lapply(x,
                            function(xElt) {
                                estimate.mean.fraglen(split(ifelse(strand(xElt) == "-",
                                                                   end(ranges(xElt)),
                                                                   start(ranges(xElt))),
                                                            strand(xElt)),
                                                      method = method, ...)
                            }))
          })

setMethod("estimate.mean.fraglen", "AlignedRead",
          function(x, method = c("SISSR", "coverage", "correlation"), ...) {
              splitData <-
                split(data.frame(strand = strand(x),
                                 start =
                                 ifelse(strand(x) == "-",
                                        position(x) + width(x) - 1L,
                                        position(x))),
                      chromosome(x))
              unlist(lapply(splitData,
                            function(y) {
                                estimate.mean.fraglen(split(y[["start"]],
                                                            y[["strand"]]),
                                                      method = method, ...)
                            }))
           })
