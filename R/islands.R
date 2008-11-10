
logscale.components <- function(axis = c("x", "y"), base = 2)
{
    axis <- match.arg(axis)
    switch(axis,
           x = 
           function(...) {
               ans <- xscale.components.default(...)
               ans$bottom$labels$labels <- 
                   base^(ans$bottom$labels$at)
               ans
           },
           y = 
           function(...) {
               ans <- yscale.components.default(...)
               ans$left$labels$labels <- 
                   base^(ans$left$labels$at)
               ans
           })
}


## common infrastructure to summarize reads in the form of nested lists: 
##  reads.list=list("1" = list("chr1" = list("+"=..., "-"=...), 
##                             "chr2"=...), 
##                  "2" = list(...))

## NOTE: summary.fun must return a data frame in the current
## implementation (and shouldn't be very big at that)

summarizeLane <- function(clist, summary.fun, ...)
{
    ## clist is a list at the lane level, with one list("+"=, "-"=) for each chromsome
    ans <- do.call(lattice::make.groups, lapply(clist, summary.fun, ...))
    names(ans)[names(ans) == "which"] <- "chromosome"
    ## cbind(chr = factor(colnames(ans), levels = colnames(ans)), as.data.frame(t(ans)))
    ans
}

summarizeLane2 <- function(clist, summary.fun, ..., seqlen)
{
    ## clist is a list at the lane level, with one list("+"=, "-"=) for each chromsome
    stopifnot(all(names(clist) %in% names(seqlen)))
    seqlen <- seqlen[names(clist)]
    mapply(summary.fun, clist, seqlen, ..., SIMPLIFY = FALSE)
}



summarizeReads <- 
    function(reads.list, lanes = names(reads.list), ..., verbose = FALSE)
{
    if (verbose) cat(paste("Processing lanes", paste(lanes, collapse = ",")), fill = TRUE)
    ans <- do.call(lattice::make.groups, lapply(reads.list[lanes], summarizeLane, ...))
    names(ans)[names(ans) == "which"] <- "lane"
    ans
}

## different summary.fun can give different useful summaries    

countSummary <- function(x) 
{
    ## x is a list at the lane->chromosome level, with components "+" and "-"
    npos <- length(x$"+")
    nneg <- length(x$"-")
    data.frame(n = npos + nneg, d = npos - nneg, r = npos / nneg)
}    

## get islands.  But this needs more care; 
## coverage by chromosome, possibly different species.  
## We'll avoid this by just using the furthest hit

## library("BSgenome.Mmusculus.UCSC.mm9")
## library("BSgenome.Hsapiens.UCSC.hg18")

## mouse.seqlens <- seqlengths(Mmusculus)
## human.seqlens <- seqlengths(Hsapiens)

sliceSummary <- 
    function(x, lower = 1,
             viewSummary = list(sums = viewSums, maxs = viewMaxs))
    ## x is a list at the lane->chromosome level, with components "+" and "-"
{
    g <- growSeqs(x)
    cov <- coverage(g, 1, max(end(g) + 400L))
    s <- slice(cov, lower = lower)
    ans <- data.frame(start = start(s), end = end(s))
    if (is.list(viewSummary)) 
    {
        for (nm in names(viewSummary))
            ans[[nm]] <- viewSummary[[nm]](s)
    }
    else ans[["summary"]] <- viewSummary(s)
    ans
}


## Not really a summary function. FIXME: coverage should eventually be
## generic, dealing with 'AlignedList'-type objects

## coverageSummary <-
##     function(x, max = NULL, strand = c("+", "-"))
##     ## x is a list at the lane->chromosome level, with components "+" and "-"
## {
##     g <- extendReads(x, strand=strand)
##     if( missing(max) ) max = max(end(g)) + 400L
##     coverage(g, 1, max)
## }


getSingletons <-
    function(x, jitter = TRUE)
    ## x is a list at the lane->chromosome level, with components "+" and "-"
{
    g <- growSeqs(x)
    cov <- coverage(g, 1, max(end(g) + 400L))
    s <- slice(cov, lower = 1)
    s <- s[viewMaxs(s) == 1]
    ## We retain length-400 islands (actually 2 adjoint length 200), but
    ## these are hopefully rare enough not to matter.
    if (jitter)
        0.5 * (start(s) + end(s) + runif(length(s)))
    else
        0.5 * (start(s) + end(s))
}

count.singletons <- function(x) ## x <- summarizeReads(.)
{
    ans <-
        do.call(make.groups,
                lapply(x,
                       function(x) {
                           as.data.frame.table(sapply(x, length))
                       }))
    names(ans)[names(ans) == "Freq"] <- "count"
    names(ans)[names(ans) == "Var1"] <- "chromosome"
    names(ans)[names(ans) == "which"] <- "lane"
    ans
}

