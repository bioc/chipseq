
setClass("AlignedList", contains = "list")

setClass("GenomeList",
         contains = "list",
         representation(genome = "character"))


## Long-term goal: use BSgenome classes "GenomeData", "GenomeDataList"

## setAs("list", "ChromosomeData",
##       function(from) {
##           new("ChromosomeData",
##               elements = from)
##       })


## setAs("list", "SampleChromosomeData",
##       function(from) {
##           el <- lapply(from, function(x) as(x, "ChromosomeData"))
##           new("SampleChromosomeData", elements = el)
##       })

## 


setGeneric("gdApply",
           function(X, FUN, ...) {
               standardGeneric("gdApply")
           })

setMethod("gdApply",
          signature(X = "GenomeDataList", FUN = "function"),
          function(X, FUN, ...) {
              NX <- length(X)
              new.elements <- vector(mode = "list", length = NX)
              names(new.elements) <- names(X)
              for (i in seq_len(NX))
              {
                  new.elements[[i]] <- gdApply(X[[i]], FUN, ...)
              }
              cls <- lapply(new.elements, class)
              ucl <- unique(unlist(cls))
              if (identical(ucl, "GenomeData"))
                  GenomeDataList(new.elements)
              else new.elements
          })

setMethod("gdApply",
          signature(X = "GenomeData", FUN = "function"),
          function(X, FUN, ...) {
              NX <- length(X)
              new.elements <- vector(mode = "list", length = NX)
              names(new.elements) <- names(X)
              for (i in seq_len(NX))
              {
                  new.elements[[i]] <- FUN(X[[i]], ...)
              }
              cls <- lapply(new.elements, class)
              ucl <- unique(unlist(cls))
              if (length(ucl) == 1)
                  GenomeData(new.elements, organism = X@organism)
              else new.elements
          })

setAs("GenomeData", "data.frame",
      function(from) {
          ans <- 
              do.call(rbind, 
                      sapply(names(from),
                             function(chr) {
                                 cbind(as(from[[chr]], "data.frame"), chromosome = chr)
                             }, simplify = FALSE))
          row.names(ans) <- NULL
          ans
      })

setAs("GenomeDataList", "data.frame",
      function(from) {
          ans <- 
              do.call(rbind, 
                      sapply(names(from),
                             function(sample) {
                                 cbind(as(from[[sample]], "data.frame"), sample = sample)
                             }, simplify = FALSE))
          row.names(ans) <- NULL
          ans
      })






## summarizeLane <- function(clist, summary.fun, ...)
## {
##     ## clist is a list at the lane level, with one list("+"=, "-"=) for each chromsome
##     ans <- do.call(lattice::make.groups, lapply(clist, summary.fun, ...))
##     names(ans)[names(ans) == "which"] <- "chromosome"
##     ## cbind(chr = factor(colnames(ans), levels = colnames(ans)), as.data.frame(t(ans)))
##     ans
## }

## summarizeLane2 <- function(clist, summary.fun, ..., seqlen)
## {
##     ## clist is a list at the lane level, with one list("+"=, "-"=) for each chromsome
##     stopifnot(all(names(clist) %in% names(seqlen)))
##     seqlen <- seqlen[names(clist)]
##     mapply(summary.fun, clist, seqlen, ..., SIMPLIFY = FALSE)
## }



## summarizeReads <- 
##     function(reads.list, lanes = names(reads.list), ..., verbose = FALSE)
## {
##     if (verbose) cat(paste("Processing lanes", paste(lanes, collapse = ",")), fill = TRUE)
##     ans <- do.call(lattice::make.groups, lapply(reads.list[lanes], summarizeLane, ...))
##     names(ans)[names(ans) == "which"] <- "lane"
##     ans
## }

