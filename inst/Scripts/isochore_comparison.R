
##

library(Matrix)
library(IRanges)
library(chipseq)


if (FALSE)
{
    
    data(isochores.mm8)
    library(BSgenome.Mmusculus.UCSC.mm8)
    isochores.df <- split(isochores.mm8, isochores.mm8$chromosome)

    gcContent <- function(isochore, chrom = "chr1")
    {
        v <- with(isochore[[chrom]],
                  Views(Mmusculus[[chrom]], start = Begin, end = End))
        rowSums(alphabetFrequency(v, baseOnly = TRUE, freq = TRUE)[, c("G", "C")])
    }

    for (chr in names(isochores.df))
    {
        print(chr)
        isochores.df[[chr]]$newGC <-
            gcContent(isochore = isochores.df, chrom = chr)
    }
    
    usedf <- do.call(rbind, isochores.df)
    
    with(usedf, 
     {
         xyplot(newGC ~ GC,
                panel = function(...) {
                    panel.smoothScatter(...)
                    panel.loess(..., col.line = "black")
                })
     })
   

}


data(isochores.mm9)

## str(isochore.mm9) # error because @annotation slot missing

dfparams <- RDApplyParams(isochores.mm9, as.data.frame)

isochore.df <- rdapply(dfparams) ## per-chromosome list
isochore.df <- isochore.df[paste("chr", 1:19, sep = "")]



load("myodMyo.rda")
load("myodFibro.rda")

load("simulatedReadsSampled.rda")


combineLaneReads <- function(laneList, chromList = names(laneList[[1]])) {
    names(chromList) = chromList ##to get the return value named
    lapply(chromList,
           function(chr) { ## sample() to make order random
               list("+" = sample(unlist(lapply(laneList, function(x) x[[chr]][["+"]]), use.names = FALSE)),
                    "-" = sample(unlist(lapply(laneList, function(x) x[[chr]][["-"]]), use.names = FALSE)))
           })
}

combinedMyo <- 
    list(cblasts = combineLaneReads(myodMyo[c("1","3","6")]),
         ctubes = combineLaneReads(myodMyo[c("2","4","7")]))

combinedFibro <- 
    list(cfibro = combineLaneReads(myodFibro[c("1","3","6")]),
         cfibromyod = combineLaneReads(myodFibro[c("2","4","7")]))


## recompute GC content

library(BSgenome.Mmusculus.UCSC.mm9)

gcContent <- function(isochore = isochore.df,
                      chrom = "chr1")
{
    v <- with(isochore[[chrom]],
              Views(Mmusculus[[chrom]], start = start, end = end))
    rowSums(alphabetFrequency(v, baseOnly = TRUE, freq = TRUE)[, c("G", "C")])
}


countOverlap <- function(isochore = isochore.df,
                         reads,
                         chrom = "chr1",
                         singletons = TRUE)
{
    isoranges <- with(isochore[[chrom]], IRanges(start, end))
    readranges <- extendReads(reads[[chrom]])
    if (singletons)
    {
        s <- slice(coverage(readranges, width = max(end(readranges)) + 300L), lower = 1)
        readranges <- s[viewMaxs(s) == 1]
    }
    readranges <- readranges[order(start(readranges))]
    ans <- overlap(isoranges, readranges, multiple = TRUE)
    rowSums(ans@matchMatrix)
    ## almost the same as 
    ##   table(overlap(isoranges, readranges, multiple = FALSE))
    ## except that 0-count ones will be missing
}


for (chrom in names(isochore.df))
{
    message("Processing ", chrom)
    isochore.df[[chrom]]$newGC <-
        gcContent(isochore = isochore.df, chrom = chrom)
    isochore.df[[chrom]]$cblast <- 
        countOverlap(isochore = isochore.df,
                     reads = combinedMyo$cblasts,
                     chrom = chrom)
    isochore.df[[chrom]]$ctube <- 
        countOverlap(isochore = isochore.df,
                     reads = combinedMyo$ctubes,
                     chrom = chrom)
    isochore.df[[chrom]]$sim <- 
        countOverlap(isochore = isochore.df,
                     reads = simulatedReads,
                     chrom = chrom, singletons = FALSE)
}

## usedf <- isochore.df$chr2

usedf <- do.call(rbind, isochore.df)

pdf("isochore_comparison.pdf")

with(usedf, 
 {
     xyplot(newGC ~ vals,
            panel = function(...) {
                panel.smoothScatter(...)
                panel.loess(..., col.line = "black")
            })
 })


with(usedf, 
 {
     logWidth <- equal.count(log(width), 9)
     xyplot(asinh(cblast + ctube) ~ newGC | logWidth,
            ## type = c("p", "g", "smooth"),
            panel = function(...) {
                panel.smoothScatter(...)
                panel.loess(..., col.line = "black")
            })
 })

## with(usedf,
##  {
##      logWidth <- equal.count(log(width), 9)
##      xyplot(asinh(sim) + asinh(cblast + ctube) ~ newGC | logWidth,
##             type = c("p", "g"))
##  })


with(usedf,
 {
     logWidth <- equal.count(log(width), 9)
     xyplot(asinh((cblast + ctube) / sim) ~ newGC | logWidth, subset = sim > 2,
            ## type = c("p", "g", "smooth"),
            panel = function(...) {
                panel.smoothScatter(...)
                panel.loess(..., col.line = "black")
            })
 })




with(usedf,
 {
     logWidth <- equal.count(log(width), 9)
     xyplot(asinh(cblast + ctube) ~ asinh(sim) | logWidth, subset = sim > 2,
            ## type = c("p", "g", "smooth"),
            panel = function(...) {
                panel.smoothScatter(...)
                panel.loess(..., col.line = "black")
            })
 })



with(usedf,
 {
     xyplot(asinh(cblast) + asinh(ctube) + asinh(sim) ~ log(width),
            outer = TRUE,
            ## type = c("p", "g", "smooth"),
            panel = function(...) {
                panel.smoothScatter(...)
                panel.loess(..., col.line = "black")
            })
 })


dev.off()


