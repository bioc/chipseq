
 library("chipseq")

 library("BSgenome.Mmusculus.UCSC.mm9")

  ##lanes 1, 3, 6 are Myoblasts
  ##lanes 2, 4, 7 are Myotubes
  ##lane 8 is a reference lane

  lanes <- c(1, 2, 3, 4, 6, 7, 8)
 
  reads = vector("list", length = length(lanes))
  names(reads) = as.character(lanes)

  for (i in seq_along(lanes) ) {
    lane <- lanes[i]
    message("Starting Lane ", lane)
    pat = paste("s_", lane, ".map", sep="")
    reads[[i]] = readAndClean("/home/jdavison/ycao/26-06-2008/binary",
       pattern = pat, exclude = "[MXY]|rand")
  }

 ##we drop the sex chromosomes and mitochondria.
 ## mouse has 19 chromosomes

 chrom.list <- paste("chr", c(1:19), sep = "")
 nchrom = length(chrom.list)
 chromLens = rep(NA, nchrom)
 names(chromLens) = chrom.list
 for( i in 1:nchrom) 
    chromLens[i] = nchar(unmasked(Mmusculus[[chrom.list[i]]])) 

 seqRanges = lapply(reads, growSeqs)


 cblasts = combineLanes(seqRanges[c(1,3,5)])
 ctubes = combineLanes(seqRanges[c(2,4,6)])

 all = combineLanes(list(cblasts, ctubes))

 ss1 = laneSubsample(cblasts, ctubes)

 covblasts = laneCoverage(cblasts, chromLens)

 blastIslands = islands(covblasts)

 blastIcts = readsPerIsland(blastIslands)

