
## FIXME: update to use filtering in ShortRead

library(ShortRead)
library(chipseq)

readReads <-
    function(srcdir, lane, ...,
             exclude = "[M]|rand", type = "MAQMapShort",
             simplify = TRUE)
{
    filt <-
        compose(strandFilter(strandLevels=c("-", "+")),
                chromosomeFilter(regex = "chr[0-9]+$"),
                uniqueFilter(withSread = FALSE),
                alignQualityFilter(15),
                ...)
    message(sprintf("reading data from lane %s [%s], using filter %s", lane, srcdir, name(filt)))
    ans <- readAligned(srcdir, lane, type = type, filter = filt)
    if (simplify) as.list(ans)
    else ans
    ## readAndClean(srcdir, lane, exclude = exclude, type = type, ...))
}

## In ShortRead now

## idFilter <- 
##     function (regex = character(0), fixed = FALSE, .name = "IdFilter") 
## {
##     stopifnot(length(regex) == 1)
##     srFilter(function(x) {
##         .idx <- logical(length(x))
##         .idx[grep(regex, as.character(id(x)), fixed=fixed)] <- TRUE
##         .idx
##     }, name = .name)
## }

## readFirstRead <-
##     function(srcdir = "/home/jdavison/ycao/01-09-2008/text", lane,
##              exclude = "[MXY]|rand", minScore = 15, dropDups = TRUE)
## {
##     message(sprintf("reading data from lane %s [%s]", lane, srcdir))
##     aln <- readAligned(srcdir, lane, type = "MAQMapview")
##     ## trimmed to the first read of the pair
##     aln <- aln[grep("/1", as.character(id(aln)), fixed=TRUE)]
##     exChr <- grep(exclude, aln@chromosome)
##     aln <- aln[-exChr]
##     keep <- (aln@alignQuality@quality >= minScore)
##     s2 <- aln[keep]
##     if (dropDups) 
##         as.list(s2[!srduplicated(sread(s2))])
##     else
##         as.list(s2)
## }

## paired end reads

## lane1: mouse fibroblasts expressing Myod, 3 antibodies combined
## lane2: C2C12 myotube, 3 antibodies combined
## lane3: human fibroblast expressing Myod, antibody 7311
## lane4: human fibroblast expressing Myod, antibody 6975b

## lane6: human fibroblast expressing Myod, antibody 6196
## lane7: human fibroblast controls, 3 antibodies combined
## lane8: human fibroblast expressing Myod, 3 antibodies combined


pat.lanes <- sprintf("s_%g", 1:8)
names(pat.lanes) <- as.character(1:8)
pat.lanes <- pat.lanes[-5]

pairedReads <-
    lapply(pat.lanes,
           function(s) {
               readReads(srcdir = "/home/jdavison/ycao/01-09-2008/text", lane = s,
                         type = "MAQMapview",
                         idFilter(regex = "/1", fixed = TRUE))
           })

save(pairedReads, file = "pairedReads.rda")
rm(pairedReads)
gc()



pat.lanes <- sprintf("s_%g", 1:8)
names(pat.lanes) <- as.character(1:8)
pat.lanes <- pat.lanes[-5]

pairedFragmentSummary <-
    function(srcdir, lane)
{
    filt <-
        compose(strandFilter(strandLevels=c("-", "+")),
                chromosomeFilter(regex = "chr[0-9]+$"),
                uniqueFilter(withSread = FALSE),
                alignQualityFilter(15))
    message(sprintf("reading data from lane %s [%s], using filter %s", lane, srcdir, name(filt)))
    ans <- readAligned(srcdir, lane, type = "MAQMapview", filter = filt)

    ids <- as.character(id(ans))
    match2 <- grep("/2", ids)
    match.names <- gsub("/2", "/1", ids[match2])
    match1 <- match(match.names, ids)
    smry <- 
        data.frame(strand1 = strand(ans)[match1],
                   chrom1 = chromosome(ans)[match1],
                   position1 = position(ans)[match1],
                   strand2 = strand(ans)[match2],
                   chrom2 = chromosome(ans)[match2],
                   position2 = position(ans)[match2])
    smry <- subset(smry, (chrom1 == chrom2) & (strand1 != strand2))
    smry$length <-
        with(smry, 
             35L + ifelse(strand1 == "-", position1 - position2, position2 - position1))
                   
##     ## add 35?
##     bwplot(chrom1 ~ length, smry, xlim = c(-500, 500))
##     summary(smry$length)

    smry
}

pairedFragmentLengths <-
    lapply(pat.lanes,
           function(s) {
               pairedFragmentSummary(srcdir = "/home/jdavison/ycao/01-09-2008/text", lane = s)
           })

save(pairedFragmentLengths, file = "pairedFragmentLengths.rda")
rm(pairedFragmentLengths)
gc()


paired.frag.size <- 
    do.call(lattice::make.groups, 
            lapply(pairedFragmentLengths,
                   function(x) {
                       m <- with(x, tapply(length, chrom1, median))
                       data.frame(chrom = names(m),
                                  mu.est = as.numeric(m))
                   }))

stripplot(reorder(which, mu.est) ~ mu.est, paired.frag.size,
          jitter = TRUE)






## myoD_myo
## lanes 1, 3, 6 are Myoblasts
## lanes 2, 4, 7 are Myotubes
## lane 8 is a reference lane

myodMyo <-
    lapply(pat.lanes,
           function(s) {
               readReads(srcdir = "/home/jdavison/ycao/26-06-2008/binary", lane = s,
                         type = "MAQMapShort")
           })

save(myodMyo, file = "myodMyo.rda")
rm(myodMyo)
gc()

## myoD_myo
## lanes 1, 3, 6 are Fibroblasts
## lanes 2, 4, 7 are Fibroblasts + MyoD
## lane 8 is a reference lane

myodFibro <-
    lapply(pat.lanes,
           function(s) {
               readReads(srcdir = "/home/jdavison/ycao/25-04-2008/binary", lane = s,
                         type = "MAQMapShort")
           })

## myodFibro <- lapply(pat.lanes, function(s) readReads(srcdir =
## "/home/jdavison/ycao/25-04-2008/binary", lane = paste(s, ".map", sep = "")))

save(myodFibro, file = "myodFibro.rda")
rm(myodFibro)
gc()


sessionInfo()


