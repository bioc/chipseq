
## calculations on successive subsets to assess if enough data has been collected


estimate.bg.rate <- function(s, seqLen)
{
    ## y <- table(viewMaxs(s))
    y1 <- sum(viewMaxs(s) == 1)
    y2 <- sum(viewMaxs(s) == 2)
    p <- y2 / y1
    alpha <- 1 - (1-p)^(1/seqLen)
    alpha
}

summary.subsets <- function(...)
{
    .Deprecated("subsetSummary")
    subsetSummary(...)
}

cutoff.fdr <- function(cov, fdr.cutoff = 0.001, k = 2:20)
{
    ## an implementation of the idea in Robertson et al to assess
    ## sufficiency of sampling depth: choose minimum cutoff that gives
    ## an FDR < pre-specified value, then compute number of peaks at
    ## that cutoff.
    s <- slice(cov, lower = 1)
    y <- table(viewMaxs(s))
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



subsetSummary <- 
    function(x,
             chr,
             nstep, ## number of reads in each increment for full data
             props = seq(.1, 1, .1),
             chromlens,
             fg.cutoff = 6, seqLen = 200,
             fdr.cutoff = 0.001,
             use.fdr = FALSE,
             resample = TRUE, islands = TRUE,
             verbose = getOption("verbose"))
{
    g <- extendReads(x[[chr]], seqLen = seqLen)
    if (resample) g <- g[sample(length(g))]
    ##     if (!missing(nstep) && missing(props))
    ##         props <- seq(nstep / length(g), 1, by = nstep / length(g))
    if (!missing(nstep) && missing(props))
    {
        nreads.total <-
          sum(unlist(lapply(x, function(u) sum(sapply(u, length)))))
        num.steps <- floor(nreads.total / nstep)
        if (verbose) message(num.steps, " steps\n.")
        props <- seq(0, 1, length.out = num.steps)[-1]
    }
    ids <- as.integer(round(props * length(g)))
    if (verbose)
      message(length(g), " reads in ", chr, ". Increments: ",
              paste(ids, collapse = ", "))
    old.peaks <- IRanges()
    start <- 1L
    ans.cols <-
        c("alpha.hat", "bg.rate", "old.bg", "old.fg",
          "new.bg", "new.fg", "old.fg.area", "old.total.area",
          "reads.converted", "npeaks", "npeaks.fdr", "npeaks.fdr.lower",
          "npeaks.fdr.higher", "fdr.count.cutoff")
    ans <- matrix(0, nrow = length(ids), ncol = length(ans.cols))
    colnames(ans) <- ans.cols
    for (i in seq_along(ids))
    {
        old.reads <- sort(head(g, start-1L))
        cum.reads <- head(g, ids[i])
        if (length(old.reads) > 0)
        {
            old.islands <- slice(coverage(old.reads, width = chromlens[chr]),
                                 lower = 1)
            nreads <- as.integer(viewSums(old.islands) / seqLen)
            old.total.area <- sum(width(old.islands))
        }
        else
        {
            old.total.area <- 0
            old.islands <- IRanges()
        }
        ## old fg area
        old.fg.area <- sum(width(old.peaks))
        new.reads <- sort(g[start:(ids[i])]); start <- ids[i] + 1L
        current.cov <- coverage(cum.reads, width = chromlens[chr])
        fdr.interp <- cutoff.fdr(current.cov, fdr.cutoff = fdr.cutoff)
        fdr.floor <- floor(fdr.interp)
        fdr.ceiling <- ceiling(fdr.interp)
        current.islands <- slice(current.cov, lower = 1)
        alpha.hat <- estimate.bg.rate(current.islands, seqLen = seqLen)
        bg.rate <- alpha.hat / ids[i]
        if (islands)
        {
            peaks.fixed <-
              current.islands[viewMaxs(current.islands) >= fg.cutoff]
            peaks.fdr.lower <- current.islands[viewMaxs(current.islands) >=
                                               fdr.ceiling]
            peaks.fdr.higher <- current.islands[viewMaxs(current.islands) >=
                                                fdr.floor]
        }
        else 
        {
            peaks.fixed <- slice(current.cov, lower = fg.cutoff)
            peaks.fdr.lower <- slice(current.cov, lower = fdr.ceiling)
            peaks.fdr.higher <- slice(current.cov, lower = fdr.floor)
        }
        peaks.fdr <- peaks.fdr.lower
        npeaks.fdr <- length(peaks.fdr.lower) + (fdr.ceiling - fdr.interp) *
          (length(peaks.fdr.higher) - length(peaks.fdr.lower))
        current.peaks <- peaks.fixed
        if (use.fdr)
          current.peaks <- peaks.fdr
        ## old reads that hit old peaks
        old.peak.hits.old <- old.reads %in% old.peaks
        ## old reads that hit current peaks
        current.peak.hits.old <- old.reads %in% current.peaks
        ## old background hits
        old.bg <- sum(!old.peak.hits.old)
        old.fg <- sum(old.peak.hits.old)
        ## number of old reads that go from bg to fg
        reads.converted <- sum(current.peak.hits.old & !old.peak.hits.old)
        ## new reads that hit old fg
        old.peak.hits.new <- new.reads %in% old.peaks
        ## new reads that hit old something (that is, not blank)
        old.total.hits.new <- new.reads %in% old.islands
        new.fg <- sum(old.peak.hits.new)
        new.bg <- sum(old.total.hits.new & !old.peak.hits.new)
        ##print(table(old.total.hits.new, old.peak.hits.new))
        ##browser()
        ans[i, ] <- c(alpha.hat, bg.rate, old.bg, old.fg,
                      new.bg, new.fg, old.fg.area, old.total.area,
                      reads.converted, length(peaks.fixed), npeaks.fdr,
                      length(peaks.fdr.lower), length(peaks.fdr.higher),
                      fdr.ceiling)
        old.peaks <- current.peaks
    }
    ans <- cbind(chromosome = chr, proportion = props, size = diff(c(0, ids)),
                 cumsize = ids, as.data.frame(ans))
    ans
}


