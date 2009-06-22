
## calculations on successive subsets to assess if enough data has been collected


estimate.bg.rate <- function(s, seqLen)
{
    y <- table(viewMaxs(s))
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

subsetSummary <- 
    function(x,
             chr,
             nstep, ## number of reads in each increment for full data
             props = seq(.1, 1, .1),
             chromlens,
             fg.cutoff = 6, seqLen = 200,
             resample = TRUE, islands = TRUE,
             verbose = getOption("verbose"))
{
    g <- extendReads(x[[chr]], seqLen = seqLen)
    if (resample) g <- g[sample(length(g))]
    ##     if (!missing(nstep) && missing(props))
    ##         props <- seq(nstep / length(g), 1, by = nstep / length(g))
    if (!missing(nstep) && missing(props))
    {
        nreads.total <- sum(unlist(lapply(x, function(u) sum(sapply(u, length)))))
        num.steps <- floor(nreads.total / nstep)
        if (verbose) message(num.steps, " steps\n.")
        props <- seq(0, 1, length.out = num.steps)[-1]
    }
    ids <- as.integer(round(props * length(g)))
    if (verbose) message(length(g), " reads in ", chr, ". Increments: ", paste(ids, collapse = ", "))
    old.peaks <- IRanges()
    start <- 1L
    ans.cols <-
        c("alpha.hat", "bg.rate", "old.bg", "old.fg",
          "new.bg", "new.fg", "old.fg.area", "old.total.area",
          "reads.converted", "npeaks")
    ans <- matrix(0, nrow = length(ids), ncol = length(ans.cols))
    colnames(ans) <- ans.cols
    for (i in seq_along(ids))
    {
        old.reads <- sort(head(g, start-1L))
        cum.reads <- head(g, ids[i])
        if (length(old.reads) > 0)
        {
            old.islands <- slice(coverage(old.reads, width = chromlens[chr]), lower = 1)
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
        if (islands)
        {
            current.peaks <- slice(coverage(cum.reads, width = chromlens[chr]), lower = 1)
            alpha.hat <- estimate.bg.rate(current.peaks, seqLen = seqLen)
            bg.rate <- alpha.hat / ids[i]
            current.peaks <- current.peaks[viewMaxs(current.peaks) >= fg.cutoff]
            current.peaks <- as(current.peaks, "IRanges")
        }
        else 
        {
            alpha.hat <- estimate.bg.rate(slice(coverage(cum.reads, width = chromlens[chr]), lower = 1),
                                          seqLen = seqLen)
            bg.rate <- alpha.hat / ids[i]
            current.peaks <- as(slice(coverage(head(g, ids[i]), width = chromlens[chr]),
                                      lower = fg.cutoff), "IRanges")
        }
        ## old reads that hit old peaks
        old.peak.hits.old <- !is.na(overlap(old.peaks, old.reads, multiple = FALSE))
        ## old reads that hit current peaks
        current.peak.hits.old <- !is.na(overlap(current.peaks, old.reads, multiple = FALSE))
        ## old background hits
        old.bg <- sum(!old.peak.hits.old)
        old.fg <- sum(old.peak.hits.old)
        ## number of old reads that go from bg to fg
        reads.converted <- sum(current.peak.hits.old & !old.peak.hits.old)
        ## new reads that hit old fg
        old.peak.hits.new <- !is.na(overlap(old.peaks, new.reads, multiple = FALSE))
        ## new reads that hit old something (that is, not blank)
        old.total.hits.new <- !is.na(overlap(old.islands, new.reads, multiple = FALSE))
        new.fg <- sum(old.peak.hits.new)
        new.bg <- sum(old.total.hits.new & !old.peak.hits.new)
        ##print(table(old.total.hits.new, old.peak.hits.new))
        ##browser()
        ans[i, ] <- c(alpha.hat, bg.rate, old.bg, old.fg,
                      new.bg, new.fg, old.fg.area, old.total.area,
                      reads.converted, length(current.peaks))
        old.peaks <- current.peaks
    }
    ans <- cbind(chromosome = chr, proportion = props, size = diff(c(0, ids)),
                 cumsize = ids, as.data.frame(ans))
    ans
}


