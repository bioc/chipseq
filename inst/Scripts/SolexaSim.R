 library(ShortRead)

# sp <- SolexaPath("/mnt/fred/solexa/ycao/080623_HWI-EAS88_0001")
# rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")
# abc <- alphabetByCycle(rfq)

# save(abc, file="abc.rda")

 load("abc.rda")

 q <- rawToChar(as.raw(32+0:93)) # these are valid quality encodings
 intQ = as(SFastqQuality(q), "numeric") + 33 # actually, 'integer' :( 

 x1 = apply(abc, c(1,3), 
           function(x) if( sum(x) > 0 ) round(sum(x*intQ)/sum(x)) else 0)

 x1 = x1[ rowSums(x1) != 0, ]
 x1 = x1[-5,]  #drop the N's

 x2 = ss[x1]
 dim(x2) = dim(x1)
 dimnames(x2) = dimnames(x1)

 library(BSgenome.Mmusculus.UCSC.mm9)

 chr1 = unmasked(Mmusculus[[X]])

 myreads = sample(1:nchar(chr1), 500000)

 myV = Views(chr1, start=myreads, end=myreads+34)

 SimChars = as.character(myV)

 rm(chr1)

 Char2Row = sapply(SimChars, function(x) chartr("ACGTN", "12345", x))

 SimQ = sapply(Char2Row, function(x) {
            r = as.integer(unlist(strsplit(x, "")))
	    return(x1[matrix(c(r,1:35), nc=2)]) })


            

