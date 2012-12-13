# Time-stamp: < modified by qinqianhappy :2012-12-13 20:04:03 >
## For analysis TCGA RNA-seq Data

require(DESeq)


ratioHist <- function(){
   # get an example count data set -- or use your data:
   cds <- makeExampleCountDataSet()
   cds <- newCountDataSet( countsTable, conds )
   # estimate the size factors:
   cds <- estimateSizeFactors( cds )

   # calculate the gene-wise geometric means
   geomeans <- exp( rowMeans( log( counts(cds) ) ) )

   # choose a sample whose size factor estimate we want to check
   j <- 1

   # just for clarity: the size factor was estimated as described above,
   # so the following two lines give the same result
   print( sizeFactors(cds)[j] )
   print( median( ( counts(cds)[,j]/geomeans )[ geomeans>0 ] )  )

   # Plot a histogram of the ratios of which we have just taken
   # the median:
   hist( log2( counts(cds)[,j] / geomeans ), breaks=100 )

   # This histogram should be unimodal, with a clear peak at the value
   # of the size factor. Mark it in red:
   abline( v=log2( sizeFactors(cds)[ j ] ), col="red" )
}
