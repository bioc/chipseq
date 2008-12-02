
library(ShortRead)
library(chipseq)

simulatedReads <-
    readAligned("/home/dsarkar/chipseq-simulation/",
                pattern = "allChr.map",
                type = "MAQMapShort",
                filter = alignQualityFilter(15))

simulatedReads <- as.list(simulatedReads)

save(simulatedReads, file = "simulatedReads.rda")
sessionInfo()

## for across chromosome comparison, we should have started with
## number of reads proportional to chromosome lengths.  Instead, we
## had the same number from all chromosomes.

## sampling the results proportional to chromosome lengths will not
## work because the effective lengths of chromosomes are not the same
## as official lengths (due to gaps, repeats, etc.).  However, we can
## estimate this "thinning" effect by the proportion of original reads
## for which alignments were found.

library(BSgenome.Mmusculus.UCSC.mm9)

ref.chromlens <- seqlengths(Mmusculus)[names(simulatedReads)]

aligned.reads <- sapply(simulatedReads, function(x) length(x[["+"]]))
## aligned.props <- aligned.reads / 1e6

keep.props <- (aligned.reads/1e6) * (ref.chromlens/1e6)
keep.props <- keep.props / max(keep.props)

## for each chromosome, choose a suitable subset

for (chr in names(simulatedReads))
{
    cat(chr, fill = TRUE)
    n <- aligned.reads[chr]
    x <- rbinom(1, size = n, prob = keep.props[chr])
    simulatedReads[[chr]][["+"]] <-
        sample(simulatedReads[[chr]][["+"]], x)
}

save(simulatedReads, file = "simulatedReadsSampled.rda")

