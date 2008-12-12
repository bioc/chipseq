
setClass("AlignedList", contains = "list")

setClass("GenomeList",
         contains = "list",
         representation(genome = "character"))


## Long-term goal: use BSgenome classes "ChromosomeData", "SampleChromosomeData"

setAs("list", "ChromosomeData",
      function(from) {
          new("ChromosomeData",
              elements = from)
      })


setAs("list", "SampleChromosomeData",
      function(from) {
          el <- lapply(from, function(x) as(x, "ChromosomeData"))
          new("SampleChromosomeData", elements = el)
      })


