
library(chipseq)
library("BSgenome.Mmusculus.UCSC.mm9")

data(abc)
qs = seqQScores(abc)

schr =  names(Mmusculus)[1:21]

simv = as.integer(0.02 * seqlengths(Mmusculus)[schr])

names(simv) = schr

set.seed(123)
sims24 = simulateReads(simv, Mmusculus, 24, qs[,3:26])

printSim(sims24, file="sim24.fastq")
