

source("/home/dsarkar/svn/Projects/eboxes/context.R")

if (file.exists("geneMouse.rda")) load("geneMouse.rda") else
{
    geneMouse <-
        read.table("geneMouse.txt", header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
    save(geneMouse, file = "geneMouse.rda")
}

gregions <- genomic_regions(genes = geneMouse)
gregions$gene <- as.character(gregions$gene)

library(lattice) # for make.groups
library(Matrix)
library(IRanges)

load("peakSummary.rda")


countHits <- function(subject, query)
{
##     tree <- IntervalTree(subject)
##     result.nm <- overlap(tree, query, multiple = FALSE)
##     ## result <- overlap(tree, query)
##     str(result)
    sum(!is.na(overlap(subject, query, multiple = FALSE)))
}


doPeakSet <- function(peaks, gregions)
{
    query <- with(peaks, IRanges(start, end))
    irangeByType <-
        function(type = c("promoter", "threeprime",
                          "upstream", "downstream", "gene"))
        {
            type <- match.arg(type)
            istarts <- sprintf("%s.start", type)
            iends <- sprintf("%s.start", type)
            IRanges(start = unique(gregions[[istarts]]),
                    end = unique(gregions[[iends]]))
        }
    subject <-
        list(promoter = irangeByType("promoter"),
             threeprime = irangeByType("threeprime"),
             upstream = irangeByType("upstream"),
             downstream = irangeByType("downstream"),
             gene = irangeByType("gene"))
    c(total = length(query), sapply(subject, countHits, query = query))
}

doChromosome <- function(chr)
{
    gregions.sub <- subset(gregions, chrom == chr)
    tube.peaks <-
        subset(peakSummary.blasts.wrt.tubes,
               chromosome == chr)
    blast.peaks <-
        subset(peakSummary.tubes.wrt.blasts,
               chromosome == chr)
    ans <-
        as.data.frame(rbind(tube = doPeakSet(tube.peaks, gregions.sub),
                            blast = doPeakSet(blast.peaks, gregions.sub),
                            tube.only = doPeakSet(subset(tube.peaks, resids < -3), gregions.sub),
                            blast.only = doPeakSet(subset(blast.peaks, resids < -3), gregions.sub)))
    ans <- cbind(type = rownames(ans), ans)
    ans
}

## all.chroms <- levels(peakSummary.blasts.wrt.tubes$chromosome)

doAll <- function(chroms = paste("chr", 1:19, sep = ""))
{
    ans <- do.call(make.groups, sapply(all.chroms, doChromosome, simplify = FALSE))
    names(ans)[names(ans) == "which"] <- "chromosome"
    rownames(ans) <- NULL
    ans
}

doAll()
       
