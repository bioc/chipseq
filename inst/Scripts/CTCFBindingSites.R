# Load required libraries
library("Biostrings")
library("BSgenome.Mmusculus.UCSC.mm9")

# Create forward and reverse Position Weight Matrices for CTCF binding
ctcfForwardPWM <-
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
rownames(ctcfForwardPWM) <- c("A", "C", "G", "T")
ctcfForwardPWM <- 100 * ctcfForwardPWM

ctcfReversePWM <- ctcfForwardPWM[c("T", "G", "C", "A"), ncol(ctcfForwardPWM):1]
rownames(ctcfReversePWM) <- c("A", "C", "G", "T")

# Match CTCF binding to mouse genome (takes 8 - 10 minutes)
unmaskedMatchPWM <- function(pwm, subject, min.score = "80%") {
    active(masks(subject)) <- FALSE
    matchPWM(pwm = pwm, subject = subject, min.score = min.score)
}
bsParams <- new("BSParams", X = Mmusculus, FUN = unmaskedMatchPWM, simplify = FALSE)
ctcfForwardMatches <- bsapply(bsParams, pwm = ctcfForwardPWM, min.score = "80%")
ctcfReverseMatches <- bsapply(bsParams, pwm = ctcfReversePWM, min.score = "80%")

# Summarize the output
# Counts
sapply(ctcfForwardMatches, length)
sapply(ctcfReverseMatches, length)

# Consensus matrices from the matches
lapply(ctcfForwardMatches, consensusMatrix, baseOnly = TRUE)
lapply(ctcfReverseMatches, consensusMatrix, baseOnly = TRUE)
