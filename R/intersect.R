

effective.glength <- function(IR1, IR2, width = 200L)
{
    
    ## Motivation: think of IR1 as ranges of singleton islands in a
    ## sample.  IR2 represents similar ranges in a different sample.
    ## Assuming random distribution, the chance that IR2[i] intersects
    ## IR1 is, up to the assumption that elements of IR1 are well
    ## separated, is p = sum(width(IR1) + width) / G, where width is
    ## the width of IR2[i].  Since all our widths are the same, a
    ## simple approximation is

    ## p = 2 * width * length(IR1) / G

    ## Then, the number of elements in IR2 that intersects, X ~ Bin(m, p)

    ## where m = length(IR2).  This may be used to estimate G.

    if (!is(IR1, "IRanges")) IR1 <- IRanges(start = as.integer(IR1 - width/2L), width = width)
    if (!is(IR2, "IRanges")) IR2 <- IRanges(start = as.integer(IR2 - width/2L), width = width)
    m <- overlap(IR1, IR2, multiple = FALSE)
    phat <- sum(!is.na(m)) / length(m)
    Ghat <- 2 * width * length(IR1) / phat
    Ghat
}


effective.glength.byChr <-
    function(IRL1, IRL2, chroms = intersect(names(IRL1), names(IRL2)))
{
    sapply(chroms,
           function(chr) {
               effective.glength(IRL1[[chr]],
                                 IRL2[[chr]])
           })
}


