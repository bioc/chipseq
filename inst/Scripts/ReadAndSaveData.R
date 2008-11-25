
## FIXME: update to use filtering in ShortRead

library(ShortRead)
library(chipseq)

readReads <-
    function(srcdir, lane, exclude = "[MXY]|rand", type = "MAQMapShort", ...)
{
    message(sprintf("reading data from lane %s [%s]", lane, srcdir))
    ans <- 
        readAligned(srcdir, lane, type = type, 
                    filter = compose(strandFilter(strandLevels=c("-", "+")),
                                     chromosomeFilter(regex = "chr[0-9]+"),
                                     uniqueFilter(withSread = FALSE),
                                     alignQualityFilter(15)))
    as.list(ans)
    ## readAndClean(srcdir, lane, exclude = exclude, type = type, ...))
}



readFirstRead <-
    function(srcdir = "/home/jdavison/ycao/01-09-2008/text", lane,
             exclude = "[MXY]|rand", minScore = 15, dropDups = TRUE)
{
    message(sprintf("reading data from lane %s [%s]", lane, srcdir))
    aln <- readAligned(srcdir, lane, type = "MAQMapview")
    ## trimmed to the first read of the pair
    aln <- aln[grep("/1", as.character(id(aln)), fixed=TRUE)]
    exChr <- grep(exclude, aln@chromosome)
    aln <- aln[-exChr]
    keep <- (aln@alignQuality@quality >= minScore)
    s2 <- aln[keep]
    if (dropDups) 
        as.list(s2[!srduplicated(sread(s2))])
    else
        as.list(s2)
}


pat.lanes <- sprintf("s_%g", 1:8)
names(pat.lanes) <- as.character(1:8)
pat.lanes <- pat.lanes[-5]

## paired end reads

## lane1: mouse fibroblasts expressing Myod, 3 antibodies combined
## lane2: C2C12 myotube, 3 antibodies combined
## lane3: human fibroblast expressing Myod, antibody 7311
## lane4: human fibroblast expressing Myod, antibody 6975b

## lane6: human fibroblast expressing Myod, antibody 6196
## lane7: human fibroblast controls, 3 antibodies combined
## lane8: human fibroblast expressing Myod, 3 antibodies combined


pairedReads <- lapply(pat.lanes, function(s) readFirstRead(lane = s))
save(pairedReads, file = "pairedReads.rda")
rm(pairedReads)
gc()


## myoD_myo
## lanes 1, 3, 6 are Myoblasts
## lanes 2, 4, 7 are Myotubes
## lane 8 is a reference lane

myodMyo <- lapply(pat.lanes, function(s) readReads(srcdir = "/home/jdavison/ycao/26-06-2008/binary", lane = paste(s, ".map", sep = "")))
save(myodMyo, file = "myodMyo.rda")
rm(myodMyo)
gc()

## myoD_myo
## lanes 1, 3, 6 are Fibroblasts
## lanes 2, 4, 7 are Fibroblasts + MyoD
## lane 8 is a reference lane

myodFibro <- lapply(pat.lanes, function(s) readReads(srcdir = "/home/jdavison/ycao/25-04-2008/binary", lane = paste(s, ".map", sep = "")))
save(myodFibro, file = "myodFibro.rda")
rm(myodFibro)
gc()


sessionInfo()


