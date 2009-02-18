library(chipseq)

## also read 24-mer CTCF & GFP data
readUniqueMappings <-
function (srcdir, lane, ..., include = "chr[0-9]+$", type = "MAQMapShort", 
          simplify = TRUE)
{
    filters <-
      compose(strandFilter(strandLevels = c("-", "+")),
              chromosomeFilter(regex = include),
              ...)
    message(sprintf("reading data from lane %s [%s], using filter %s", 
                    lane, srcdir, name(filters)))
    ans <- readAligned(srcdir, lane, type = type, filter = filters)
    if (!all(is.na(quality(alignQuality(ans)))))
        ans <- ans[order(chromosome(ans), strand(ans), position(ans),
                         - quality(alignQuality(ans)))]
    ans <-
      ans[!duplicated(data.frame(chromosome(ans), strand(ans), position(ans)))]
    if (simplify) {
        ans <- XDataFrame(chromosome = Rle(chromosome(ans)),
                          strand = Rle(strand(ans)),
                          start = ifelse(strand(ans) == "-",
                                         position(ans) + width(ans) - 1L,
                                         position(ans)),
                          quality = quality(alignQuality(ans)))
    } else {
        ans@sread <- DNAStringSet(as.character(ans@sread))
        ans@quality@quality <- BStringSet(as.character(ans@quality@quality))
    }
    ans
}

maqMaps <- c("SRR001985.map", "SRR001986.map", "SRR001987.map")
names(maqMaps) <- as.character(seq_len(length(maqMaps)))

bowtieMaps <- c("SRR001985.bowtie_map", "SRR001986.bowtie_map", "SRR001987.bowtie_map")
names(bowtieMaps) <- as.character(seq_len(length(bowtieMaps)))

ctcf <-
  list("maq" =
       do.call(XDataFrameList, lapply(maqMaps, function(s)
           readUniqueMappings(srcdir = "/home/jdavison/externalData/ES_CTCF/bases3-26/maq/maps",
                             lane = s, type = "MAQMapShort"))),
       "bowtie"=
       do.call(XDataFrameList, lapply(bowtieMaps, function(s)
           readUniqueMappings(srcdir = "/home/paboyoun/externalData/ES_CTCF/bases3-26/bowtie/maps",
                              lane = s, type = "Bowtie"))))

sapply(ctcf, nrow)

save(ctcf, file = "ctcf.rda")
rm(ctcf)
gc()


maqMaps <- c("SRR001996.map", "SRR001997.map", "SRR001998.map", "SRR001999.map")
names(maqMaps) <- as.character(seq_len(length(maqMaps)))

bowtieMaps <- c("SRR001996.bowtie_map", "SRR001997.bowtie_map", "SRR001998.bowtie_map", "SRR001999.bowtie_map")
names(bowtieMaps) <- as.character(seq_len(length(bowtieMaps)))

gfp <-
  list("maq" =
       do.call(XDataFrameList, lapply(maqMaps, function(s)
           readUniqueMappings(srcdir = "/home/jdavison/externalData/ES_CTCF/bases3-26/maq/GFP_background/maps",
                              lane = s, type = "MAQMapShort"))),
       "bowtie" =
       do.call(XDataFrameList, lapply(bowtieMaps, function(s)
           readUniqueMappings(srcdir = "/home/paboyoun/externalData/ES_CTCF/bases3-26/bowtie/GFP_background/maps",
                              lane = s, type = "Bowtie"))))

sapply(gfp, nrow)

save(gfp, file = "gfp.rda")
rm(gfp)
gc()


sessionInfo()
