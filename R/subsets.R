
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

### FIXME: simplify/update this function:
## 'chr' is unncessary (why always limit to one chromosome?)
## 'chromlens' comes from our GRanges
## 'seqLen' is unnecessary (just resize first)
##
### OR: just remove it -- have we ever seen saturation?
subsetSummary <- 
    function(x,
             chr,
             nstep, ## number of reads in each increment for full data
             props = seq(.1, 1, .1),
             chromlens = seqlengths(x),
             fg.cutoff = 6, seqLen = 200,
             fdr.cutoff = 0.001,
             use.fdr = FALSE,
             resample = TRUE, islands = TRUE,
             verbose = getOption("verbose"))
{
    .Defunct()
    g <- GRanges(chr, ranges(resize(x[seqnames(x) == chr], width=seqLen)))
    seqlengths(g) <- chromlens[chr]
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
    old.peaks <- RangesList()
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
        old.reads <- head(g, start-1L)
        old.reads <- old.reads[order(start(old.reads))]
        cum.reads <- head(g, ids[i])
        if (length(old.reads) > 0) {
            old.islands <- slice(coverage(old.reads), lower = 1)
            nreads <- as.integer(unlist(viewSums(old.islands) / seqLen,
                                        use.names=FALSE))
            old.islands <- ranges(old.islands)
            old.total.area <- sum(width(old.islands))
            old.fg.area <- sum(width(old.peaks))
        } else {
            old.total.area <- 0
            old.fg.area <- 0
            old.islands <- RangesList()
        }
        new.reads <- g[start:(ids[i])]
        new.reads <- new.reads[order(start(new.reads))]
        start <- ids[i] + 1L
        current.cov <- coverage(cum.reads)
        fdr.interp <- peakCutoff(current.cov, fdr.cutoff = fdr.cutoff)
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
        npeaks.fdr <- sum(elementLengths(peaks.fdr.lower)) +
          (fdr.ceiling - fdr.interp) *
          (sum(elementLengths(peaks.fdr.higher)) -
           sum(elementLengths(peaks.fdr.lower)))
        current.peaks <- peaks.fixed
        if (use.fdr)
          current.peaks <- peaks.fdr
        ## old reads that hit old peaks
        old.peak.hits.old <- !is.na(findOverlaps(old.reads, old.peaks,
                                                 select="first"))
        ## old reads that hit current peaks
        current.peak.hits.old <- !is.na(findOverlaps(old.reads,
                                                     ranges(current.peaks),
                                                     select="first"))
        ## old background hits
        old.bg <- sum(!old.peak.hits.old)
        old.fg <- sum(old.peak.hits.old)
        ## number of old reads that go from bg to fg
        reads.converted <- sum(current.peak.hits.old & !old.peak.hits.old)
        ## new reads that hit old fg
        old.peak.hits.new <- !is.na(findOverlaps(new.reads, old.peaks,
                                                 select="first"))
        ## new reads that hit old something (that is, not blank)
        old.total.hits.new <- !is.na(findOverlaps(new.reads, old.islands,
                                                  select="first"))
        new.fg <- sum(old.peak.hits.new)
        new.bg <- sum(old.total.hits.new & !old.peak.hits.new)
        ##print(table(old.total.hits.new, old.peak.hits.new))
        ##browser()
        ans[i, ] <- c(alpha.hat, bg.rate, old.bg, old.fg,
                      new.bg, new.fg, old.fg.area, old.total.area,
                      reads.converted, sum(elementLengths(peaks.fixed)),
                      npeaks.fdr,
                      sum(elementLengths(peaks.fdr.lower)),
                      sum(elementLengths(peaks.fdr.higher)),
                      fdr.ceiling)
        old.peaks <- ranges(current.peaks)
    }
    ans <- cbind(chromosome = chr, proportion = props, size = diff(c(0, ids)),
                 cumsize = ids, as.data.frame(ans))
    ans
}


