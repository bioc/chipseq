
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
                alignQualityFilter(9),
                ...)
    message(sprintf("reading data from lane %s [%s], using filter %s", lane, srcdir, name(filt)))
    ans <- readAligned(srcdir, lane, type = type, filter = filt)
    if (simplify) as.list(ans)
    else ans
}





## First 54-cycle read

## lane1: C2C12 0h DNA methylation ChIP (rep1)
## lane2: C2C12 0h DNA methylation ChIP (rep 2)
## lane3: C2C12 96h DNA methylation ChIP (rep 1)
## lane4: C2C12 96h DNA methylation ChIP (rep 2)
## lane5: phiX
## lane6: Rhabdomyosarcoma(RD) 24h Myod ChIP -- RD is a human muscle tumor with defects in muscle differentiation program
## lane7: primary mouse myotubes 72h Myod ChIP (Antibody: 6975)
## lane8: primary mouse myotubes 72h Myod ChIP (Antibody: 6196)

pat.lanes <- sprintf("s_%g.map", 1:8)
names(pat.lanes) <- as.character(1:8)
pat.lanes <- pat.lanes[-c(5)] 

solexa54 <-
    lapply(pat.lanes,
           function(s) {
               readReads(srcdir = "/home/jdavison/ycao/29-12-2008/binary", lane = s,
                         type = "MAQMapShort")
           })


sapply(solexa54, function(x) sum(sapply(x, function(u) sum(sapply(u, length)))))

save(solexa54, file = "solexa54.rda")
rm(solexa54)
gc()


## also read CTCF data

pat.lanes <- c("SRR001985.map", "SRR001986.map", "SRR001987.map")
names(pat.lanes) <- as.character(1:3)

ctcf <-
    lapply(pat.lanes,
           function(s) {
               readReads(srcdir = "/home/jdavison/externalData/ES_CTCF/maps", lane = s,
                         type = "MAQMapShort")
           })

sapply(ctcf, function(x) sum(sapply(x, function(u) sum(sapply(u, length)))))

save(ctcf, file = "ctcf.rda")




sessionInfo()




readOne <-
    function(srcdir, lane, ...,
             exclude = "[M]|rand", type = "MAQMapShort",
             simplify = TRUE)
{
    filt <-
        compose(strandFilter(strandLevels=c("-", "+")),
                chromosomeFilter(regex = "chr[0-9]+$"),
                uniqueFilter(withSread = FALSE),
                ##  alignQualityFilter(10),
                ...)
    message(sprintf("reading data from lane %s [%s], using filter %s", lane, srcdir, name(filt)))
    ans <- readAligned(srcdir, lane, type = type, filter = filt)
    if (simplify) as.list(ans)
    else ans
}

lane2 <- readOne(srcdir = "/home/jdavison/ycao/29-12-2008/binary",
                 lane = "s_2",
                 type = "MAQMapShort", simplify = FALSE)
lane6 <- readOne(srcdir = "/home/jdavison/ycao/29-12-2008/binary",
                 lane = "s_6",
                 type = "MAQMapShort", simplify = FALSE)
lane8 <- readOne(srcdir = "/home/jdavison/ycao/29-12-2008/binary",
                 lane = "s_8",
                 type = "MAQMapShort", simplify = FALSE)

barchart(table(quality(alignQuality(lane2))))
barchart(table(quality(alignQuality(lane6))))
barchart(table(quality(alignQuality(lane8))))

