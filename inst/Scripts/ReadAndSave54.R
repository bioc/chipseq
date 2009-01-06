
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

pat.lanes <- sprintf("s_%g", 1:8)
names(pat.lanes) <- as.character(1:8)
pat.lanes <- pat.lanes[-c(5, 6, 7, 8)] # just 5 eventually, but 6 and 7 haven't finished yet

solexa54 <-
    lapply(pat.lanes,
           function(s) {
               readReads(srcdir = "/home/jdavison/ycao/29-12-2008/binary", lane = s,
                         type = "MAQMapShort")
           })

save(solexa54, file = "solexa54.rda")
rm(solexa54)
gc()



sessionInfo()




readOne <-
    function(srcdir, lane, ...,
             exclude = "[M]|rand", type = "MAQMapShort",
             simplify = TRUE)
{
    filt <-
        compose(strandFilter(strandLevels=c("-", "+")),
                chromosomeFilter(regex = "chr[0-9]+$"),
                ## uniqueFilter(withSread = FALSE),
                ##  alignQualityFilter(10),
                ...)
    message(sprintf("reading data from lane %s [%s], using filter %s", lane, srcdir, name(filt)))
    ans <- readAligned(srcdir, lane, type = type, filter = filt)
    if (simplify) as.list(ans)
    else ans
}

foo <- readOne(srcdir = "/home/jdavison/ycao/29-12-2008/binary",
               lane = "s_8",
               type = "MAQMapShort", simplify = FALSE)
