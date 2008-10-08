
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
    pat <- paste("s_", lane, ".map", sep="")
    reads[[i]] <- readAndClean("/home/jdavison/ycao/26-06-2008/binary",
                               pattern = pat, exclude = "[MXY]|rand")
}

## we drop the sex chromosomes and mitochondria.
## mouse has 19 chromosomes

chrom.list <- paste("chr", c(1:19), sep = "")

## nchrom <- length(chrom.list)
## chromLens = rep(NA, nchrom)
## names(chromLens) = chrom.list
## for( i in 1:nchrom) 
##     chromLens[i] = nchar(unmasked(Mmusculus[[chrom.list[i]]])) 

chromLens <-
    sapply(chrom.list,
           function(chr) {
               nchar(unmasked(Mmusculus[[chr]]))
               ## same as 'length(Mmusculus[[chr]])' ?
           },
           simplify = TRUE)

system.time(seqRanges.old <- lapply(reads, growSeqs), gcFirst=TRUE)

## basically same, but retains order of chromosomes
system.time(seqRanges <- lapply(lapply(reads, as.list), growSeqs), gcFirst=TRUE)



cblasts = combineLanes(seqRanges[c(1,3,5)])
ctubes = combineLanes(seqRanges[c(2,4,6)])

covblasts = laneCoverage(cblasts[chrom.list], chromLens)
covtubes = laneCoverage(ctubes[chrom.list], chromLens)
covctrl = laneCoverage(seqRanges[["8"]], chromLens)


tubes.islands <- islandSummary(islands(covtubes))
save(tubes.islands, file = "tubes.islands.rda")


 all = combineLanes(list(cblasts, ctubes))

 ss1 = laneSubsample(cblasts, ctubes)

 blastIslands = islands(covblasts)

 blastIcts = readsPerIsland(blastIslands)

 blastSummaries = islandSummary(blastIslands)



 ##now lets see how to do the subtraction that Zizhen is doing
 ## basic work flow: given two IRanges objects, find peaks in
 ## one,  merge those that are close, and then look in the
 ## other to see how many reads there are in that one in the same
 ## location -  would this be more easily done by some form of density
 ## estimation? Does that remove the bias we see due to differences in
 ## counts for each lane?


 blastp12 = lapply(covblasts, slice, lower=12)
 mblastp12 = lapply(blastp12, merge, maxgap=100)

 blastinctrl = copyIRangesbyChr(mblastp12, covctrl)

 maxBlast = lapply(mblastp12, viewMaxs)
 maxCtrl =  lapply(blastinctrl, viewMaxs)

 ratios = vector("list", length=length(maxBlast))
 names(ratios) = names(maxBlast)
 for(i in names(maxBlast)) 
     ratios[[i]] = maxCtrl[[i]]/maxBlast[[i]]

sapply(ratios, function(x) sum(x > 1/6))  ##these are the bad ones

 tubesp12 = lapply(covtubes, slice, lower = 12)
 mtubep12 = lapply(tubesp12, merge, maxgap=100)


##read in Zizhen's coverage vector - in some compressed form
if(FALSE) {

 if(FALSE) {
 ##read in only those records for Chr2
 my = read.delim("/home/yzizhen/chip_seq/myoD_myo/s_1_3_6.coverage",
 as.is=TRUE, skip=763857,nrows=815599, head=F)

 Z1 = rep(0, 181748087)

 for(i in 1:nrow(my)) {
     Z1[(my[i,2]):(my[i,3])] = my[i,5]
     if(  (i %% 10000000) == 0 ) print(i)
 }

} else  load("Z1.rda")
 
 my2 = as(covblasts[["chr2"]], "integer")
 Z1 = Z1[1:chromLens["chr2"]] ##since Zizhen wrote off the end

 diff1 = my2-Z1

 nonz = which(diff1 != 0 )

 mstr = toupper("cttgtgggcacagctcgtgggcacagcagccctgt")
 

}
