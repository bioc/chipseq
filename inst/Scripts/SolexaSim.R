library(chipseq)
library("BSgenome.Mmusculus.UCSC.mm9")

data(abc)
qs = seqQScores(abc)

schr =  names(Mmusculus)[1:21]

simv = as.integer(0.02 * seqlengths(Mmusculus)[schr])

names(simv) = schr

set.seed(123)

N = 100
for (i in seq_len(length(N))) {
    cat("Iteration", i, "\n")
    sims24 = simulateReads(simv/N, Mmusculus, 24, qs[,3:26])
    writeFastq(sims24, file = paste("sim24_", i, ".fastq", sep = ""))
}
