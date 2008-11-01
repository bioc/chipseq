
 library(chipseq)
            
 library("BSgenome.Mmusculus.UCSC.mm9")

 data(abc)
 qs = seqQScores(abc)

 schr =  names(Mmusculus)[1:21]

 simv = rep(1000, length(schr))

 names(simv) = schr

 sims = simulateReads(simv, Mmusculus, 35, qs)

 printSim(sims, file="sim.out")

