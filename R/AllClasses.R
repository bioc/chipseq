
setClass("AlignedList", contains = "list")

setClass("GenomeList",
         contains = "list",
         representation(genome = "character"))

