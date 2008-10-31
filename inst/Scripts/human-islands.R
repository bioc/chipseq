

library(chipseq)
load("pairedReads.rda")


summarizeLane <- function(clist, summary.fun, ..., seqlen)
{
    ## clist is a list at the lane level, with one list("+"=, "-"=) for each chromsome
    stopifnot(all(names(clist) %in% names(seqlen)))
    seqlen <- seqlen[names(clist)]
    mapply(summary.fun, clist, seqlen, ..., SIMPLIFY = FALSE)
}


summarizeReads <- 
    function(reads.list, lanes = c("1", "2", "3", "4", "6", "7", "8"), ...,
             verbose = TRUE)
{
    if (verbose) cat(paste("Processing lanes", paste(lanes, collapse = ",")), fill = TRUE)
    lapply(reads.list[lanes], summarizeLane, ...)
}


coverageSummary <- 
    function(x, max = max(end(g)) + 400L)
    ## x is a list at the lane->chromosome level, with components "+" and "-"
{
    g <- growSeqs(x)
    coverage(g, 1, max)
}



## paired end reads

## lane1: mouse fibroblasts expressing Myod, 3 antibodies combined
## lane2: C2C12 myotube, 3 antibodies combined
## lane3: human fibroblast expressing Myod, antibody 7311
## lane4: human fibroblast expressing Myod, antibody 6975b

## lane6: human fibroblast expressing Myod, antibody 6196
## lane7: human fibroblast controls, 3 antibodies combined
## lane8: human fibroblast expressing Myod, 3 antibodies combined

library("BSgenome.Hsapiens.UCSC.hg18")
human.seqlens <- seqlengths(Hsapiens)

pairedEndHumanCoverage <- 
    summarizeReads(pairedReads,
                   lanes = c("3", "4", "6", "7", "8"),
                   summary.fun = coverageSummary,
                   seqlen = human.seqlens)


save(pairedEndHumanCoverage, file = "pairedEndHumanCoverage.rda")

isldf <- 
    do.call(make.groups, 
            lapply(pairedEndHumanCoverage,
                   function(x) viewMaxs(slice(x[["chr2"]], lower = 1))))

dotplot(xtabs(~data + which, isldf),
        horizontal = FALSE, groups = FALSE,
        scales = list(y = list(log = 2), x = list(rot = 90)),
        aspect = "xy")


dotplot(xtabs(~data + which, isldf),
        horizontal = FALSE, type = "o", 
        par.settings = simpleTheme(pch = 16),
        auto.key = list(columns = 5),
        aspect = "xy",
        scales = list(y = list(log = 2), x = list(rot = 90)))

## conclusion: 10 seems to be a good coverage cutoff

