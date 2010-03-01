## Load required libraries
library("chipseq")
library("BSgenome.Mmusculus.UCSC.mm9")

## create Position Weight Matrix for CTCF binding
ctcfPWM <-
  t(rbind(c(0.22, 0.24, 0.33, 0.22),
          c(0.06, 0.35, 0.21, 0.37),
          c(0.20, 0.19, 0.41, 0.19),
          c(0.28, 0.06, 0.57, 0.09),
          c(0.05, 0.87, 0.05, 0.03),
          c(0.02, 0.96, 0.01, 0.01),
          c(0.67, 0.02, 0.14, 0.16),
          c(0.05, 0.53, 0.41, 0.01),
          c(0.13, 0.54, 0.12, 0.21),
          c(0.90, 0.02, 0.04, 0.05),
          c(0.01, 0.01, 0.97, 0.01),
          c(0.39, 0.01, 0.59, 0.01),
          c(0.03, 0.02, 0.62, 0.33),
          c(0.02, 0.01, 0.97, 0.01),
          c(0.02, 0.03, 0.91, 0.04),
          c(0.17, 0.72, 0.02, 0.10),
          c(0.39, 0.02, 0.51, 0.08),
          c(0.04, 0.52, 0.41, 0.03),
          c(0.15, 0.38, 0.14, 0.34),
          c(0.31, 0.28, 0.33, 0.08)))
rownames(ctcfPWM) <- c("A", "C", "G", "T")
ctcfPWM <- 100 * ctcfPWM

## match CTCF binding to mouse genome (takes 8 - 10 minutes)
matchPWMList <- function(pwm, subject, min.score = "80%") {
    revPWM <- pwm[c("T", "G", "C", "A"), ncol(pwm):1]
    rownames(revPWM) <- c("A", "C", "G", "T")

    subject <- unmasked(subject)
    list("+" = start(matchPWM(pwm = pwm, subject = subject, min.score = min.score)),
         "-" = end(matchPWM(pwm = revPWM, subject = subject, min.score = min.score)))
}
bsParams <-
  new("BSParams", X = Mmusculus, FUN = matchPWMList,
      exclude = c("X","Y","M","random"), simplify = FALSE)
ctcfPWMMatches <- GenomeData(bsapply(bsParams, pwm = ctcfPWM, min.score = "85%"))
save(ctcfPWMMatches, file = "ctcfPWMMatches.rda")

ctcfPWMRanges <- extendReads(ctcfPWMMatches, seqLen = 20)
ctcfPWMScores <-
  lapply(structure(names(ctcfPWMMatches), names = names(ctcfPWMMatches)),
         function(i) {
             revPWM <- ctcfPWM[c("T", "G", "C", "A"), ncol(ctcfPWM):1]
             rownames(revPWM) <- c("A", "C", "G", "T")
             subject <- unmasked(Mmusculus[[i]])
             scores <-
               c(PWMscore(pwm = ctcfPWM, subject = subject,
                          start = ctcfPWMMatches[[i]][["+"]]),
                 PWMscore(pwm =  revPWM, subject = subject,
                          start = ctcfPWMMatches[[i]][["-"]] - ncol(ctcfPWM) + 1L))
             scores[order(c(ctcfPWMMatches[[i]][["+"]],
                            ctcfPWMMatches[[i]][["-"]] - ncol(ctcfPWM) + 1L))]
         })
do.call(rbind, lapply(ctcfPWMMatches, sapply, length))
do.call(c, lapply(ctcfPWMRanges, length))

ctcfPWMStringSet <-
  lapply(structure(names(ctcfPWMMatches), names = names(ctcfPWMMatches)),
         function(i) {
             subject <- unmasked(Mmusculus[[i]])
             strings <-
               c(as.character(Views(subject,
                   IRanges(start = ctcfPWMMatches[[i]][["+"]], width = 20))),
                 as.character(reverseComplement(DNAStringSet(as.character(Views(subject,
                   IRanges(end = ctcfPWMMatches[[i]][["-"]], width = 20)))))))
             DNAStringSet(strings[order(c(ctcfPWMMatches[[i]][["+"]],
                                          ctcfPWMMatches[[i]][["-"]] - ncol(ctcfPWM) + 1L))])
         })
fragmentLength <- 140L
readLength <- 24L
ctcfPWMPotentialRanges <-
  lapply(structure(names(ctcfPWMMatches), names = names(ctcfPWMMatches)),
         function(i) {
             list("+" =
                  IRanges(start = unlist(lapply(ctcfPWMMatches[[i]][["+"]],
                                                function(j) {
                                                (j - (fragmentLength - 1L)):(j + (20L - readLength))
                                                #(j - (fragmentLength - ncol(ctcfPWM))):(j + (ncol(ctcfPWM) - readLength))
                                                })),
                          width = readLength),
                  "-" =
                  IRanges(end = unlist(lapply(ctcfPWMMatches[[i]][["-"]],
                                              function(j) {
                                              (j + (20L - readLength) + (fragmentLength - 1L)):j
                                              #(j + (fragmentLength - ncol(ctcfPWM))):(j - (ncol(ctcfPWM) - readLength))
                                              })),
                          width = readLength))
         })
ctcfPWMPotentialReducedRanges <-
  lapply(structure(names(ctcfPWMMatches), names = names(ctcfPWMMatches)),
         function(i) {
             list("+" = reduce(ctcfPWMPotentialRanges[[i]][["+"]]),
                  "-" = reduce(ctcfPWMPotentialRanges[[i]][["-"]]))
                })
ctcfPWMPotentialReads <-
  DNAStringSet(sort(unique(unlist(
    lapply(structure(names(ctcfPWMMatches), names = names(ctcfPWMMatches)),
           function(i) {
               subject <- unmasked(Mmusculus[[i]])
               c(as.character(Views(subject, ctcfPWMPotentialRanges[[i]][["+"]])),
                 as.character(reverseComplement(DNAStringSet(as.character(Views(subject,
                                                ctcfPWMPotentialRanges[[i]][["-"]]))))))
           })))))

data(abc)
qualityScores <- seqQScores(abc)[,3:26]
ctcfPWMPotentialReadQuality <-
  SFastqQuality(BStringSet(do.call(paste,
    c(lapply(seq_len(readLength), function(i)
             qualityScores[as.character(narrow(ctcfPWMPotentialReads, i, i)), i]),
      sep=""))))

writeFastq(ShortReadQ(sread = ctcfPWMPotentialReads,
                      quality = ctcfPWMPotentialReadQuality),
           file = "ctcfPWMPotentialReads140upstream.fastq")

## Using the MAQ mappings from CTCFBindingSites.R
data(ctcfPWMMatches)
ctcfSegMaqMaps <-
  readUniqueMappings(srcdir = "/home/jdavison/simulations/maq/paboyoun/ctcf/maps",
                     lane = "ctcfPWMPotentialReads.map", type = "MAQMapShort")

ctcfSegMaqCoverage <-
  lapply(structure(as.character(unique(ctcfSegMaqMaps[["chromosome"]])),
                   names = as.character(unique(ctcfSegMaqMaps[["chromosome"]]))),
         function(i, seqLen = 140) {
             ctcfSubset <-
               ctcfSegMaqMaps[as.vector(ctcfSegMaqMaps[["chromosome"]] == i),]
             ctcfSplit <- split(ctcfSubset[["start"]], as.vector(ctcfSubset[["strand"]]))
             ctcfReads <-
               IRanges(start = c(ctcfSplit[["+"]], ctcfSplit[["-"]] - seqLen + 1L),
                       end = c(ctcfSplit[["+"]] + seqLen - 1L, ctcfSplit[["-"]]))
             coverage(ctcfReads, start = 1, end = seqlengths(Mmusculus)[[i]])
         })

ctcfMappingTables <-
list("+" =
     lapply(as.character(unique(ctcfSegMaqMaps[["chromosome"]])), function(chr) {
                mapping <-
                  IRanges(start =
                          ctcfSegMaqMaps[as.vector(ctcfSegMaqMaps[["chromosome"]] == chr) &
                                         as.vector(ctcfSegMaqMaps[["strand"]] == "+"),"start"],
                          width = readLength)
                y <- table(findOverlaps(mapping, ctcfPWMPotentialReducedRanges[[chr]][["+"]], select = "first"))
                cbind(count = unname(y),
                      total = width(ctcfPWMPotentialReducedRanges[[chr]][["+"]])[as.integer(names(y))] - readLength + 1L)
            }),
     "-" =
     lapply(as.character(unique(ctcfSegMaqMaps[["chromosome"]])), function(chr) {
                mapping <-
                  IRanges(start =
                          ctcfSegMaqMaps[as.vector(ctcfSegMaqMaps[["chromosome"]] == chr) &
                                         as.vector(ctcfSegMaqMaps[["strand"]] == "-"),"start"],
                          width = readLength)
                 y <- table(findOverlaps(mapping, ctcfPWMPotentialReducedRanges[[chr]][["-"]], select = "first"))
                 cbind(count = unname(y),
                       total = width(ctcfPWMPotentialReducedRanges[[chr]][["-"]])[as.integer(names(y))] - readLength + 1L)
             }))
