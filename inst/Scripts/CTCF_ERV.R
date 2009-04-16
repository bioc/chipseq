## Load libraries
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg18))
suppressMessages(library(rtracklayer))
suppressMessages(library(chipseq))

## Create Position Weight Matrix for CTCF binding
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

revPWM <- ctcfPWM[c("T", "G", "C", "A"), ncol(ctcfPWM):1]
rownames(revPWM) <- c("A", "C", "G", "T")

## Read in human LTR RepeatMasker regions
ervRangedData <- import("hgu_rmsk_erv.bed")
if (any(width(ranges(ervRangedData)) < 20)) {
    ervRangedData <- ervRangedData[width(ranges(ervRangedData)) >= 20, ]
}

## Read in human LINE RepeatMasker regions
lineRangedData <- import("hgu_rmsk_line.bed")
if (any(width(ranges(lineRangedData)) < 20)) {
    lineRangedData <- lineRangedData[width(ranges(lineRangedData)) >= 20, ]
}

## Generate CTCF PWM scores
scoreByPWM <-
function(chr, data, pwmLength = 20L) {
    cat("\r", chr)
    widths <- width(ranges(data)[[chr]]) - (pwmLength - 1L)
    interval <- rep.int(seq_len(length(widths)), widths)
    starts <- start(ranges(data)[[chr]])
    ends <- end(ranges(data)[[chr]])
    ats <- eval(parse(text = paste("c(", paste(starts, ":",
                                               ends - (pwmLength - 1L),
                                               sep = "", collapse = ","),
                                   ")")))
    chrString <- unmasked(Hsapiens[[chr]])
    posPWMScores <- PWMscoreStartingAt(ctcfPWM, chrString, ats)
    negPWMScores <- PWMscoreStartingAt(revPWM, chrString, ats)
    RangedData(ranges(data)[[chr]],
               Chromosome = Rle(chr, length(starts)),
               Name = data[chr][["name"]],
               PosScore = unlist(lapply(split(posPWMScores, interval), max)),
               NegScore = unlist(lapply(split(negPWMScores, interval), max)),
               Width = width(ranges(data)[[chr]]),
               GCpercent =
               rowSums(alphabetFrequency(Views(chrString, ranges(data)[[chr]]),
                                         baseOnly=TRUE, freq=TRUE)[,c("C", "G")]))
}

ervPWMTable <-
  do.call(RangedDataList,
          lapply(structure(paste("chr", 1:22, sep = ""),
                           names = paste("chr", 1:22, sep = "")),
                 scoreByPWM, data = ervRangedData))
save(ervPWMTable, file = "ervPWMTable.rda")

ervPWMScores <-
  do.call(rbind,
          lapply(ervPWMTable, function(x)
                 data.frame(name = x[["Name"]],
                            score = pmax(x[["PosScore"]], x[["NegScore"]]),
                            width = x[["Width"]],
                            gc = x[["GCpercent"]])))
ervPWMScores$name <- ervPWMScores$name[drop = TRUE]
tail(sort(tapply(ervPWMScores$score, ervPWMScores$name, function(x) mean(x > 1076))))


linePWMTable <-
  do.call(RangedDataList,
          lapply(structure(paste("chr", 1:22, sep = ""),
                           names = paste("chr", 1:22, sep = "")),
                 scoreByPWM, data = lineRangedData))
save(linePWMTable, file = "linePWMTable.rda")
linePWMScores <-
  do.call(rbind,
          lapply(linePWMTable, function(x)
                 data.frame(name = x[["Name"]],
                            score = pmax(x[["PosScore"]], x[["NegScore"]]),
                            width = x[["Width"]],
                            gc = x[["GCpercent"]])))
linePWMScores$name <- linePWMScores$name[drop = TRUE]
tail(sort(tapply(linePWMScores$score, linePWMScores$name, function(x) mean(x > 1076))))
