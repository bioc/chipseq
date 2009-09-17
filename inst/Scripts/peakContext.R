
source("/home/dsarkar/svn/Projects/eboxes/context.R")

data(geneMouse)

gregions.500 <- genomic_regions(genes = geneMouse, proximal = 500)
gregions.2000 <- genomic_regions(genes = geneMouse, proximal = 2000)

gregions <- gregions.2000
gregions$gene <- as.character(gregions$gene)

library(lattice) # for make.groups
library(Matrix)
library(IRanges)

load("peakSummary.rda")


countHits <- function(subject, query)
{
##     tree <- IntervalTree(subject)
##     result.nm <- findOverlaps(query, tree, multiple = FALSE)
##     ## result <- findOverlaps(query, tree)
##     str(result)
    sum(!is.na(findOverlaps(query, subject, multiple = FALSE)))
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
            iends <- sprintf("%s.end", type)
            keep <- !duplicated(gregions[[istarts]]) ## what's the right thing to do here???
            IRanges(start = gregions[[istarts]][keep],
                    end = gregions[[iends]][keep])
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
                            tube.only = doPeakSet(subset(tube.peaks, resids < -2), gregions.sub),
                            blast.only = doPeakSet(subset(blast.peaks, resids < -2), gregions.sub)))
    ans <- cbind(type = factor(rownames(ans), levels = unique(rownames(ans))), ans)
    ans
}

## all.chroms <- levels(peakSummary.blasts.wrt.tubes$chromosome)

doAll <- function(chroms = paste("chr", 1:19, sep = ""))
{
    ans <- do.call(make.groups, sapply(chroms, doChromosome, simplify = FALSE))
    names(ans)[names(ans) == "which"] <- "chromosome"
    rownames(ans) <- NULL
    ans
}

ans <- doAll()
       
sumtab <- 
    as.data.frame(rbind(total = xtabs(total ~ type, ans),
                        promoter = xtabs(promoter ~ type, ans),
                        threeprime = xtabs(threeprime ~ type, ans),
                        upstream = xtabs(upstream ~ type, ans),
                        downstream = xtabs(downstream ~ type, ans),
                        gene = xtabs(gene ~ type, ans)))


sumtab <-
    within(sumtab,
       {
           `tube/blast` <- tube / blast
           `tube.only/blast.only` <- tube.only / blast.only
       })

sumtab[c("tube", "blast", "tube/blast", "tube.only", "blast.only", "tube.only/blast.only")]

