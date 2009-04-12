

## The first set of functions provide a version of density for
## 'sparse' data (large stretches of empty space where density is 0;
## so large that standard density() would not be able to compute it.
## FIXME: this is definitely not the right place to put this, but I'm
## not sure what is (IRanges has the underlying classes, but that does
## not necessarily make it the right choice).

## The use-case of interest is Kharchencko et al's density-correlation
## method for estimating average fragmet length


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
    x <- sort(x)
    ix <- IRanges(x-width, x+width)
    rix <- sort(reduce(ix))
    ## we will calculate density on a range containing the data.  If
    ## necessary, we will subset later (FIXME: TODO).
    from0 <- min(from, start(rix)[1] - 10L)
    to0 <- max(to, end(rix)[length(rix)] + 10L)
    if (from > from0 || to < to0) stop("[from, to] smaller than support not implemented yet (but easy to add)")
    ## ox <- overlap(rix, x, multiple = FALSE)
    ox <- findInterval(x, start(rix)) # equivalent, but a little faster
    island.densities <- 
        if (!experimental)
        {
            tapply(x, ox, doDensity1, kernel = kernel, width = width, simplify = FALSE)
        }
        else 
        {
            dk <- dKernel(width = width, kernel = kernel)
            tapply(x, ox, doDensity2, dk = dk, width = width, simplify = FALSE)
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
    maxShift <- max(shift)
    rng <- range(unlist(x)) + c(-1, 1) * maxShift
    cov.pos <- coverage(extendReads(x, seqLen = seqLen, strand = "+"), shift = 1-rng[1], width = 1+diff(rng)) > 0
    cov.neg <- coverage(extendReads(x, seqLen = seqLen, strand = "-"), shift = 1-rng[1], width = 1+diff(rng)) > 0
    n <- diff(rng) + 1L
    ans <- shiftApply(shift, cov.pos, cov.neg, function(x, y) sum(x | y),
                      verbose = verbose)
    data.frame(mu = seqLen + shift, covered = ans / ans[1])
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
    sum(pos * neg) / sqrt(sum(pos * pos) * sum(neg * neg))
}



## this needs chromosome lengths, and includes all the 0-s on either
## side 

correlationProfile <-
    function(x, chrom,
             shift = seq(5, 300, 5),
             chrom.lengths,
             center = FALSE,
             ...)
{
    dl <- lapply(x[[chrom]], sparse.density,
                 from = 1L, to = chrom.lengths[chrom], ...)
    len <- length(dl[[1]])
    wid <- len - max(shift)
    cl <-
        sapply(shift,
               function(s) {
                   corr <- 
                       with(dl,
                            similarity.corr(subseq(`+`, start = 1L, width = wid),
                                            subseq(`-`, start = 1L + s, width = wid),
                                            center = FALSE))
                   ## if (interactive()) print(corr)
                   corr
               })
    data.frame(mu = shift, corr = cl)
}



## this version only uses the range of the data

densityCorr <- function(x, shift = seq(0, 500, 5), ...)
{
    maxShift <- max(shift)
    rng <- range(unlist(x)) + c(-1, 1) * maxShift
    dl <- lapply(x, sparse.density, from = rng[1], to = rng[2], ...)
    len <- length(dl[[1]])
    wid <- len - max(shift)
    cl <-
        sapply(shift,
               function(s) {
                   corr <- 
                       with(dl,
                            similarity.corr(subseq(`+`, start = 1L, width = wid),
                                            subseq(`-`, start = 1L + s, width = wid),
                                            center = FALSE))
                   ## if (interactive()) print(corr)
                   corr
               })
    data.frame(mu = shift, corr = cl)
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

coverage.estimate <- function(x, ...)
{
    d <- basesCovered(x, ...)
    with(d, mu[which.min(covered)])
}

correlation.estimate <- function(x, ...)
{
    d <- densityCorr(x, ...)
    with(d, mu[which.max(corr)])
}

estimate.mean.fraglen <-
    function(x, method = c("SISSR", "coverage", "correlation"),
             ...)
{
    method <- match.arg(method)
    switch(method,
           SISSR = jothi.estimate(x, ...),
           coverage = coverage.estimate(x, ...),
           correlation = correlation.estimate(x, ...))
}

    
