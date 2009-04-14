

## rdafiles <- list.files(".", pattern = "locs.*rda")

## for (file in rdafiles)
## {
##     local({
##         cat(file, fill = TRUE)
##         load(file)
##         assign(gsub(".rda", "", file, fixed = TRUE),
##                locs,
##                envir = .GlobalEnv)
##     })
## }


if (FALSE)
{

library("org.Mm.eg.db")
eg.sym <- 
    toTable(org.Mm.egSYMBOL2EG[c("Tnnc2", "Des", "Myh3", "Myog", "Mylpf",
                                 "Myl1", "Acta1", "Ckm", "Cdh15", "Myl4",
                                 "Itga7")])
eg.chrloc <-
    toTable(org.Mm.egCHRLOC[ eg.sym$gene_id ])
merge(eg.sym, eg.chrloc)

}


library(lattice)
library(simplehmm)


allLocs <- function(chr, lane, jitter = TRUE)
    ## jitter to make the locations more or less unique (otherwise we
    ## have problems downstream)
{
    load(sprintf("locs_%s_%s.rda", lane, chr))
    locs <- with(locs, c(forward, reverse))
    if (jitter) locs + runif(length(locs))
    else locs
}


base.lane <- 8
data.lane <- 2

target.hits <- 30
initial.states <- c(5, 20, 40)

## create list suitable for subsequent work.  We use sapply instead of
## lapply because it uses the first argument as names.

chrs.to.use <-
    c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", ## "chr8",
      "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
      "chr16", "chr17", "chr18", "chr19")


base.locs <-
    sapply(chrs.to.use,
           allLocs, lane = base.lane,
           simplify = FALSE)

data.locs <-
    sapply(chrs.to.use,
           allLocs, lane = data.lane,
           simplify = FALSE)

names(data.locs)
chromosomes <- names(data.locs)

if (FALSE)
{
    histogram(~ data/1e6 | which,
              data =
              make.groups(base = base.locs$chr6,
                          data = data.locs$chr6),
              nint = 1000,
              layout = c(1, 2), type = "count",
              col = "magenta", border = "transparent")



    QAhist <- function(file = "locs_3_chr6.rda", xlab = "Location (mb)", ...)
    {
        load(file)
        foo <-
            histogram(~ I(data/1e6) | which,
                      data = do.call(make.groups, locs),
                      nint = 1000,
                      main = file,
                      xlab = xlab,
                      layout = c(1, 2), type = "count",
                      col = "magenta", border = "transparent", ...)
        plot(foo)
    }

    pdf("QA-hist.pdf")
    QAhist("locs_1_chr6.rda")
    QAhist("locs_2_chr6.rda")
    QAhist("locs_3_chr6.rda")
    QAhist("locs_4_chr6.rda")
    QAhist("locs_5_chr6.rda")
    QAhist("locs_6_chr6.rda")
    QAhist("locs_7_chr6.rda")
    QAhist("locs_8_chr6.rda")
    dev.off()

}



base.hits <- sapply(base.locs, length)
data.hits <- sapply(data.locs, length)

rel <- data.hits / base.hits


getIntervals <-
    function(reflocs, obslocs,
             target = 15,
             prop = length(obslocs) / length(reflocs))
{
    nbins <- ceiling(length(reflocs) * prop / target)
    quantile(reflocs, prob = ppoints(nbins, a = 1), names = FALSE)
}


bin.list <-
    sapply(chromosomes,
           function(chr, ...)
           getIntervals(base.locs[[chr]], data.locs[[chr]], ...),
           target = target.hits,
           prop = median(rel),
           simplify = FALSE)

## str(bin.list)

datahits.list <-
    sapply(chromosomes,
           function(chr)
           as.numeric(table(cut(x = data.locs[[chr]],
                                breaks = bin.list[[chr]], labels = NULL))),
           simplify = FALSE)

str(datahits.list)


if (FALSE)
{
    ## just a histogram of the raw counts, with no smoothing

    countList <- 
        mapply(FUN =
               function(counts, bins) {
                   stopifnot(length(bins) == length(counts) + 1)
                   data.frame(bin.center = 0.5 * (bins[-length(bins)] + bins[-1]) / 1e6,
                              counts = counts)
               }, datahits.list, bin.list, SIMPLIFY = FALSE)
    countDF <- do.call(make.groups, countList)

    pdf("spikes.pdf", width = 20, height = 24)

##     xyplot(counts ~ bin.center | which, countDF, type = c("l", "g"),
##            layout = c(1, length(chromosomes)), lwd = 0.5,
##            main = sprintf("Lane %g normalized by lane %g", data.lane, base.lane),
##            scales =
##            list(x = list(tick.number = 20),
##                 y = list(relation = "free", rot = 0)),
##            strip = FALSE, strip.left = TRUE)

    xyplot(counts ~ bin.center | cut(bin.center, 10),
           countDF, subset = (which == "chr2"),
           type = c("l", "g"),
           layout = c(1, 10),
           lwd = 0.5, col = "black",
           main = sprintf("Chr 2, Lane %g normalized by lane %g", data.lane, base.lane),
           scales =
           list(x = list(relation = "sliced", tick.number = 20),
                y = list(relation = "same", rot = 0)),
           strip = FALSE)

    
    dev.off()


}



fm3 <-
    coverageHmm(datahits.list,
                family =
                hmm.family("nbinom",
                           mu = 1,
                           states.scale = initial.states,
                           size = 20,
                           states.free = TRUE))

fm3 <- update(fm3, iterations = 150, verbose = TRUE)
fm3 <- update(fm3, iterations = 2, verbose = TRUE)

summary(fm3)

pdf(sprintf("comparison_%g_%g.pdf", data.lane, base.lane),
    width = 10, height = 14)

dotplot(sort(rel), xlab = "Relative coverage per chromosome")


xyplot(fm3, decode = "viterbi",
       abscissa =
       lapply(bin.list,
             function(x) (x[-1] + x[-length(x)]) / 2e6 ),
       col = c('darkgrey', 'black'),
       ylim = c(0, 100),
       scales = list(x = list(tick.number = 20, axs = "r")),
       xlab = "Location (Mb)", main = "Decoded path (Viterbi)",
       lty = 1,
       strip.left = TRUE)

dev.off()



## print(rootogram(fm3, xlim = c(0, 50), strip = TRUE,
##                 main = "Marginal rootogram"))


## print(rootogram(fm3, marginal = FALSE,
##                 main = "Conditional rootogram",
##                 xlim = c(0, 50), strip = TRUE))



## long.output <- paste(prefix, "hmm-results-long.csv", sep = "-")
## short.output <- paste(prefix, "hmm-results-short.csv", sep = "-")


## write.hmm(fm3,
##           bins = bin.list,
##           decode = "viterbi",
##           file = long.output)

## write.hmm(fm3,
##           bins = bin.list,
##           decode = "viterbi",
##           combine = TRUE,
##           file = short.output)





## densityplot(~data, do.call(make.groups, locs),
##             groups = which,
##             n = 1000, bw = 5000,
##             plot.points = FALSE)

