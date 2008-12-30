
## get isochore tracks from http://bioinfo2.ugr.es/isochores/

## http://bioinfo2.ugr.es/isochores/UpdatedAssemblies/hg18/iso_chr1.coo


chromosomes <- paste("chr", c(1:19, "X", "Y"), sep = "")

## chr <- chromosomes[1]

isochores <-
    lapply(chromosomes, 
           function(chr) {
               message("Processing ", chr)
               srcurl <- sprintf("http://bioinfo2.ugr.es/isochores/UpdatedAssemblies/mm8/iso_%s.coo", chr)
               srcdata <- readLines(url(srcurl))
               srcdata <- srcdata[-c(1:4, length(srcdata))]
               isoch <- read.table(textConnection(paste(srcdata, collapse = "\n")))
               names(isoch) <- c("Begin", "End", "Size", "GC")
               isoch$chromosome <- chr
               isoch
           })

isochores.mm8 <- do.call(rbind, isochores)
isochores.mm8$Size <- NULL
rownames(isochores.mm8) <- NULL

save(isochores.mm8, file = "../../data/isochores.mm8.rda")


