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
    else
        ans <- ans[order(chromosome(ans), strand(ans), position(ans))]
    ans <-
      ans[!duplicated(data.frame(chromosome(ans), strand(ans), position(ans)))]
    if (simplify) {
        ans <- DataFrame(chromosome = Rle(chromosome(ans)),
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

maqMapFiles <- c("SRR001985.map", "SRR001986.map", "SRR001987.map")
names(maqMapFiles) <- sprintf("CTCF_%g", seq_len(length(maqMapFiles)))

bowtieMapFiles <- c("SRR001985.bowtie_map", "SRR001986.bowtie_map", "SRR001987.bowtie_map")
names(bowtieMapFiles) <- sprintf("CTCF_%g", seq_len(length(bowtieMapFiles)))

ctcfMaps <-
  list("maq" =
       do.call(DataFrameList, lapply(maqMapFiles, function(s)
           readUniqueMappings(srcdir = "/home/jdavison/externalData/ES_CTCF/bases3-26/maq/maps",
                             lane = s, type = "MAQMapShort"))),
       "bowtie"=
       do.call(DataFrameList, lapply(bowtieMapFiles, function(s)
           readUniqueMappings(srcdir = "/home/paboyoun/externalData/ES_CTCF/bases3-26/bowtie/maps",
                              lane = s, type = "Bowtie"))))

sapply(ctcfMaps, nrow)

save(ctcfMaps, file = "ctcfMaps.rda")
rm(ctcfMaps)
gc()


maqMapFiles <- c("SRR001996.map", "SRR001997.map", "SRR001998.map", "SRR001999.map")
names(maqMapFiles) <- sprintf("GFP_%g", seq_len(length(maqMapFiles)))

bowtieMapFiles <- c("SRR001996.bowtie_map", "SRR001997.bowtie_map", "SRR001998.bowtie_map", "SRR001999.bowtie_map")
names(bowtieMapFiles) <- sprintf("GFP_%g", seq_len(length(bowtieMapFiles)))

gfpMaps <-
  list("maq" =
       do.call(DataFrameList, lapply(maqMapFiles, function(s)
           readUniqueMappings(srcdir = "/home/jdavison/externalData/ES_CTCF/bases3-26/maq/GFP_background/maps",
                              lane = s, type = "MAQMapShort"))),
       "bowtie" =
       do.call(DataFrameList, lapply(bowtieMapFiles, function(s)
           readUniqueMappings(srcdir = "/home/paboyoun/externalData/ES_CTCF/bases3-26/bowtie/GFP_background/maps",
                              lane = s, type = "Bowtie"))))

sapply(gfpMaps, nrow)

save(gfpMaps, file = "gfpMaps.rda")
rm(gfpMaps)
gc()


sessionInfo()
