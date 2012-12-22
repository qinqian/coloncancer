###################
## main for TCGA###
###################

source("Array.R")
source("RNAseq.R")

## data.orig <- IO(args[1], args[2], args[3])
data.orig <- IO("../data/",
               "colon_cancer_TCGA_agilent_expression.xls",
               "colon_cancer_mutation_all.maf",
               "RNAseq_all_expression_Trim.txt")

################################################
## preprocessing
SeqArraymatch <- function(datalist=""){
  ## format data tables columns in uniform pattern
  ## match mutation, RNAseq and array patients ID
  # patients mutation
  mutation.patient <- substr(unique(datalist$data2$Tumor_Sample_Barcode), 1, 15)
  exp.patient <- colnames(datalist$data1)
  mutationexp.match <- match(mutation.patient, exp.patient)
  mutationexp.patient <- exp.patient[mutationexp.match]
  ## write table of matched expression data
  expseq.patient <- intersect(exp.patient, colnames(datalist$data3))
  RNAseq.matchmut <- mutation.patient %in% colnames(datalist$data3)
  RNAseq.matchmut <- match(mutation.patient[RNAseq.matchmut], colnames(datalist$data3))
  Seqmut.patient <- colnames(datalist$data3)[RNAseq.matchmut]
  SeqArrayMut <- Seqmut.patient %in% mutationexp.patient
  SeqArrayMut.patient <- Seqmut.patient[SeqArrayMut]
  ## match for expression between RNAseq and arrays
  SeqExpinter.gene <- intersect(rownames(datalist$data3), as.vector(datalist$data1[,1]))
  print(head(SeqExpinter.gene))
  ## exceptions for RNASeq annotated gene SLC35E2
  ## exp.match["SLC35E2",]
  ## Overall statistics of mutations genes in different patient's samples
  Mutgenes.freq <- sort(table(datalist$data2[,1]), decreasing = T)
  mut_class <- data.frame(genes=names(Mutgenes.freq),freq=as.numeric(Mutgenes.freq))
  ## mut filter by mutation frequency
  write.table(mut_class, file="allpatientColonMutationTop.txt", quote=F, sep="\t", col.names = T)
  mutation.table <- Mut.Table(datalist$data2, mutation.patient, mut_class)
  ## match with expression, 51 patients
  rownames(mutation.table) <- mutation.table$genes
  write.table(mutation.table, file="allpatient_colon_somatic_mutations.txt", quote=F, sep="\t", col.names = T, row.names=T)
  ## match mutation genes symbol with expression genes symbol
  ## For downstream genes analysis, use all expressed genes
  ## mutation.table denotes patient with mutation and array data
  ## mutation.tableSeqArray denotes patient with mutation,array and RNASeq
  ## Overlap RNASeq data with array and mutation
  mutation.tableSeqArray <- mutation.table[, SeqArrayMut.patient]
  match_patient_mutseq <- sort(rowSums(mutation.tableSeqArray), decreasing = T)
  match_mutclass <- data.frame(match_patient_mutseq)
  write.table(match_mutclass, file="matchpatientColonMutationTop.txt", quote=F, sep="\t", col.names = T)
  write.table(mutation.tableSeqArray, file="match_patient_colon_somatic_mutations.txt", quote=F, sep="\t", col.names = T, row.names=T)
  result <- list(mutclass=match_mutclass, mutable = mutation.tableSeqArray,
                 expseq.gene=SeqExpinter.gene, expseq.patient=expseq.patient)
  result
}

## get correlation of RNAseq and array overlapped genes
SeqExp.correlation <- function(expall="", cutoff=0.5, method="pearson"){
  ## get log2 for RNAseq RPKM part,
  SeqExp.cor <- apply(expall, 1, function(x) cor(as.numeric(x[1:(length(expall[1,])/2)]),
                                                 as.numeric(x[(length(expall[1,])/2+1):length(expall[1,])]),
                                                 method=method))
  ## cutoff 0.5
  SeqExp.cor
}

normalization <- function(x ,transpose=TRUE, log2t=TRUE, pseudo = 0.0000001){
  library(preprocessCore) ##
  if (is.vector(x))
    scale(x)
  else if (is.matrix(x) || is.data.frame(x)){
    if (log2t){
      x.norm <- normalize.quantiles(log2(as.matrix(x + pseudo)))
    }
    else x.norm <- normalize.quantiles(x + pseudo)
  }
  colnames(x.norm) <- colnames(x)
  rownames(x.norm) <- rownames(x)
  return(data.frame(x.norm))
}

################################################
## step1. match data among mutation, Seq and array, normalization
# 51  for all three data
# 53 for mutation and array
matchdata <- SeqArraymatch(datalist=data.orig)
attach(matchdata)

## array Combat, Seq quantile normalization of RPKM
## mutation type analysis
## SNP: point mutation
## separate by two alleles, normal and tumor and validated ones
## use all the validated genes
table(data.orig$data2$Variant_Type)
table(data.orig$data2["Variant_Classification"])

## use the genes with mutation above 5
## top: filtered by mutation frequency in matched 51 patients
filtered <- rownames(subset(mutclass, match_patient_mutseq >= 5))
mutfiltertype <- subset(data.orig$data2, Hugo_Symbol %in% filtered)
(table(as.vector(mutfiltertype$Hugo_Symbol)))
table(mutfiltertype$Variant_Type)
write.table(
  table(mutfiltertype$Variant_Classification), quote=F,
  file="filtered_gene_mut_type.txt"
  )

MutType <- function(genes, mut){
  ## get all types of mutations for each high frequent mutated genes
  # genes for all interesting genes
  # mut for filtered table of somatic mutation
  # mutl for
  mut.type <- list()
  mut.class <- list()
  for (g in genes){
    mut.type[[g]] = table(as.vector(subset(mut, mut$Hugo_Symbol %in% g)$Variant_Type))
    mut.class[[g]] = table(as.vector(subset(mut, mut$Hugo_Symbol %in% g)$Variant_Classification))
  }
  out <- list(tp=mut.type, cl=mut.class)
  out
}

mutct <- MutType(genes=filtered, mutfiltertype)

################################################
## step2 QC measurement
## Seq and array
## normalization of expression
## 151 patients and 16145 genes
Seq.genematch <- data.orig$data3[expseq.gene, expseq.patient]
Seq.genematchnorm <- normalization(Seq.genematch)
write.table(Seq.genematchnorm, file="../results/colon_seq_normaftermatchexp.txt", sep="\t", quote=F)

rownames(data.orig$data1) <- data.orig$data1$ProbeID
Array.genematch <- data.orig$data1[expseq.gene,expseq.patient]
write.table(Array.genematch, file="../results/colon_array_normaftermatchexp.txt", sep="\t", quote=F)

## all matched 151 patients' expression data seq and arrays
expseq.match <- cbind(data.frame(Seq.genematchnorm), Array.genematch+0.0000001)

## pearson
SeqExp.cor <- SeqExp.correlation(expseq.match, cutoff=0.5) ## cancell cutoff to get distribution of correlation
## spearman
cor.spearman <- SeqExp.correlation(expseq.match, cutoff=0.5, method="spearman") ## cancell cutoff to get distribution of correlation
## evaluate correlation
## get rid of outlier effects
pdf("../results/Spearman_Pearson.pdf", width = 25, height=10)
par(mfrow = c(1,2))
plot(SeqExp.cor, cor.spearman, main="Outlier control by Spearman correlation vs Pearson correlation",
     xlab="pearson correlation", ylab="spearman correlation", sub="16145 overlapped genes between RNAseq and Agilent array",
     pch='.', col="red")
Lab.palette <-
  colorRampPalette(c("blue", "orange", "red"), space = "Lab")
smoothScatter(SeqExp.cor, cor.spearman, main="Outlier control by Spearman correlation vs Pearson correlation",
     xlab="pearson correlation", ylab="spearman correlation", sub="16145 overlapped genes between RNAseq and Agilent array",
     pch='.')#colramp = Lab.palette)
dev.off()
pdf("../results/colon_correlation.pdf", width = 25, height=10)
par(mfrow = c(1,2))
plot.ecdf(SeqExp.cor, main="normalized RNAseq and Agilent arrays pearson correlation distribution",
          xlab="Correlation Value", ylab="Accumulative distribution", sub="Colon Cancer 51 matched patients with 16145 genes")
text(0.3,0.8, "Cutoff>=0.5, 6906 genes", cex=1)
abline(v=0.5,lty=2, col="red")
plot.ecdf(cor.spearman, main="normalized RNAseq and Agilent arrays spearman correlation distribution",
          xlab="Correlation Value", ylab="Accumulative distribution", sub="Colon Cancer 51 matched patients with 16145 genes")
text(0.3,0.8, "Cutoff>=0.5, 10235 genes", cex=1)
abline(v=0.5,lty=2, col="red")
dev.off()

## pearson
indexall <- names(SeqExp.cor)
index0.5 <- names(which(SeqExp.cor>=0.5))
index0.2 <- names(which(SeqExp.cor>=0.2))
index0.3 <- names(which(SeqExp.cor>=0.3))

## spearman
indexsall <- names(cor.spearman)
indexs0.5 <- names(which(cor.spearman>=0.5))
indexs0.3 <- names(which(cor.spearman>=0.3))
indexs0.2 <- names(which(cor.spearman>=0.2))


match3.pt <- colnames(mutable)
## non normalizaed seq with combat arrays
seq_array_non <- cbind(data.orig$data3[expseq.gene, match3.pt], data.orig$data1[expseq.gene, match3.pt])

## normalized without correlation filtering genes 16145
seq_m_norm<- normalization(data.orig$data3[expseq.gene, match3.pt])
names(seq_m_norm) <- gsub('.','-', names(seq_m_norm),fixed=T)
seq_array_mut <- cbind(seq_m_norm, data.orig$data1[expseq.gene, match3.pt])
write.table(seq_array_mut, file="../results/matchexp_seq_arraybeforefilter.txt", sep="\t", quote=F)

corcutoff.out <- function(index="", cutoff=0.5, all){
  ## normalized match data
  pdf(paste("../results/QCseqarray", cutoff, ".pdf", sep="", collapse = ""))
  par(mfrow = c(2,1))
  ## expseq.match is the all matched RNAseq and arrays
  boxplot(expseq.match[index, 1: 151], main="RNAseq expression over patients QC", col=c(4:7), outline=F)
  boxplot(expseq.match[index, 152:302], main="array expression over patients QC", col=c(4:7), outline=F)
  dev.off()
  pdf(paste("../results/QCseqarraymut", cutoff, ".pdf",sep="", collapse = ""))
  par(mfrow = c(3,1))
  boxplot(data.orig$data3[index, match3.pt], main="RNAseq expression over patients", col=c(4:7), outline=F)
  boxplot(all[index, 1:51], main="RNAseq normalization expression over patients QC", col=c(4:7), outline=F)
  boxplot(all[index, 52:102], main="array expression over patients QC", col=c(4:7), outline=F)
  dev.off()
  write.table(all[index, ], file=paste("../results/matchexp_seq_array",cutoff, '.txt', sep="_", collapse=""), sep="\t", quote=F)
  pdf(paste("../results/scatter_mean_seq_array", cutoff,".pdf", sep=""))
  scatterplot(rowMeans(all[index, 1:51]), rowMeans(all[index, 52:102]), xlab="RNAseq", ylab="Array",
              main=paste("correlation of RNAseq and Array, cutoff=", cutoff))
  dev.off()
  return(all[index,])
}

corcutoff.out(indexall, cutoff="all", seq_array_mut)
corcutoff.out(index0.5, cutoff=0.5, seq_array_mut)
corcutoff.out(index0.3, cutoff=0.3, seq_array_mut)
corcutoff.out(index0.2, cutoff=0.2, seq_array_mut)
corcutoff.out(indexs0.5, cutoff="spearman0.5", seq_array_mut)
corcutoff.out(indexs0.2, cutoff="spearman0.2", seq_array_mut)
## after spearman correlation filter and normalization
all_norm <- corcutoff.out(indexs0.3, cutoff="spearman0.3", seq_array_mut)

## check important driver genes
## Cancer Driver
## Archilles BRAF
braf <- read.delim("../../archilles/archilles_BRAF_han.txt", header=F)
names(braf) <- "genes"
Driver <- read.xls("../../archilles/Drivers, Oct 16 2012.xlsx")
driver <- Driver$HUGO
matchDriver <- function(genes,..){
  cat(length(intersect(genes, indexall)), "pearson all \n")
  cat(length(intersect(genes, index0.2)), "pearson 0.2 \n")
  cat(length(intersect(genes, index0.5)), "pearson 0.5 \n")
  cat(length(intersect(genes, index0.3)), "pearson 0.3 \n")
  cat(length(intersect(genes, indexsall)), "spearman all \n")
  cat(length(intersect(genes, indexs0.2)), "spearman 0.2 \n")
  cat(length(intersect(genes, indexs0.5)), "spearman 0.5 \n")
  cat(length(intersect(genes, indexs0.3)), "spearman 0.3 \n")
}

## correlation levels
matchDriver(braf$genes)
matchDriver(driver)

## QC of outlier point
source("../code/nscore.R")
library(copa)
library(DriverNet) ## very similar to what we are doing
seq_outlier <- getPatientOutlierMatrix(all_norm[, 1:51]) ## DriverNet
length(which(seq_outlier == TRUE))  # outlier 18116
length(which(seq_outlier == FALSE)) # 637642

test1 <- nscore(as.vector(all_norm[,2]))
bak <- backtr(test1$nscore, test1)# tails="separate")
cor(all_norm[,2], bak)

## for mutation analysis
## driver genes detected
library(seqinr)
library(CorMut)
## optional steps paradigm shift, kegggraph, cytoscape
library(iCluster)
library(som) # use som to replace iCluster knn
library(nnet)
require(CancerMutationAnalysis)

################################################
## step3, differential analysis of overlap
## mutation and arrays
metaRank <- function(rank1, rank2){
  ## meta analysis between different data type,
  ## e.g. RNASeq and Array
  sort(rank1*rank2)
}

## mutation, arrays and Seq
## mutclass, data.orig$data2, mutable, filtered 156 high frequent mutated genes
seq_array_mut.index <- apply(mutable[filtered,], 1, function(x) {which(as.numeric(x) >= 1)})
## call exp.patientclass to draw plot in batch
## classify patients' expression by mutation frequent genes and draw QQ-plot, heatmap
########################################
## downstream  matchseq analysis
## rank only by mutation frequency
################################
## TODO mutation types classification involve
################################
## topmatch is the patient mutation data (matched with seq and array)
## take 53 patients mutation frequency now,
## TODO: adjust mutation to 51 matched patients, done
## all_norm, with correlation filtered 0.3 spearman correlation
## 1:51 seq, 52: 102 array all_norm
pdf("../results/seq_colon_seq_array_match_normfilter_downstream0.01.pdf")
seq_match_common <- exp.patientclass(seq_array_mut.index,
                                     all_norm[, 1:(length(all_norm[1,])/2)],
                                     head(filtered, 2), cutoff=0.01, "downstream","matchseq")
dev.off()

pdf("../results/array_colon_seq_array_match_normfilter_downstream0.01.pdf")
array_match_common <- exp.patientclass(seq_array_mut.index,
                                       all_norm[, (length(all_norm[1,])/2+1):length(all_norm[1,])],
                                       head(filtered,2),cutoff= 0.01, "downstream","matcharray")
dev.off()

## seq not normalization without correlation filter seq_array_mut differential expressed
pdf("../results/seq_colon_seq_array_match_normfilter_downstream0.01origbefore.pdf")
seq_match_common_before <- exp.patientclass(seq_array_mut.index,
                                            seq_array_non[, c(1:51)],
                                            head(filtered, 2), cutoff=0.01, "downstream","matchseqnonfiler")
dev.off()

## array not normalization without correlation filter seq_array_mut differential expressed
pdf("../results/array_colon_seq_array_match_normfilter_downstream0.01origbefore.pdf")
array_match_common <- exp.patientclass(seq_array_mut.index,
                                       seq_array_non[, (length(all_norm[1,])/2+1):length(all_norm[1,])],
                                       head(filtered,2),cutoff= 0.01, "downstream","matcharray")
dev.off()

## take APC as examples, to calculate array and seq significant genes overlap
collective_gene <- function(gene=""){
  ## using t.test value
  array <- c()
  array[1:length(all_norm[,1])] <- "none"
  array[array_match_common[[gene]]$tu] <- "up"
  array[array_match_common[[gene]]$td] <- "down"
  seq <- c()
  seq[1:length(all_norm[,1])] <- "none"
  seq[seq_match_common[[gene]]$tu] <- "up"
  seq[seq_match_common[[gene]]$td] <- "down"

  write.table(rownames(all_norm)[intersect(array_match_common[[gene]]$tu,seq_match_common[[gene]]$tu)], file=paste("collective_up", gene, ".xls", sep=""), quote=F, sep = "\t")
  write.table(rownames(all_norm)[intersect(array_match_common[[gene]]$td,seq_match_common[[gene]]$td)], file=paste("collective_down", gene, ".xls",sep=""), quote=F, sep = "\t")
  data.frame(table(array, seq))
  write.table(data.frame(table(array, seq)), file=paste(gene,"0.01", sep = "", collapse = ""), quote=F, sep="\t")
  return(table(array, seq))
}

# filtered
APC <- collective_gene("APC")
SYNE1 <- collective_gene("SYNE1")
TTN <- collective_gene("TTN")
TP53 <- collective_gene("TP53")
APOB <- collective_gene("APOB")

## limmma
seq_match_limma <- limma.patientclass(patient.match.mutatedSeqArray,
                                      log2(all[SeqExp.cor, 1:(length(all[1,])/2)]),
                                      head(topmatch$genes), 1)


arraylimmamatchSeq.result <- limma.patientclass(patient.match.mutatedSeqArray,
                                                log2(all[SeqExp.cor, 1:(length(all[1,])/2)]),
                                                head(topmatch$genes), 1)
dev.off()
source("../code/GSEA-P-R/GSEA.1.0.R")
## exp.match is arrays matched with mutation table

##   ## array match with RNAseq 52:102
 #   pdf("../results/Colon_Cancer_Scatter_downstreammatchSeq.pdf")
  #  dev.off()
##   ## ## For feedback genes regulation analsysis, use overlapped genes from mutation and array
##   ## exp.match.withgene <- cbind(genes=rownames(exp.match), exp.match)
##   ## exp.match.overlap <- merge(Colon.Mutclass, exp.match.withgene, sort=F, all=F)
##   ## exp.match.overlaped <- exp.match.overlap[,-c(1,2)]
##   ## rownames(exp.match.overlaped) <- exp.match.overlap[,1]
##   ## pdf("../results/Colon_Cancer_Scatter_feedbackarray.pdf")
##   ## feedback_array <- exp.patientclass(patient.match.mutated, exp.match.overlaped, toparraymatch, toparraymatch$genes, 0.01, "feedback")
##   ## dev.off()
##   ## call rna.diff
##   ## limma differential expression results list
##   ## rank products
##   result <- list(seq=arraylimmamatchSeq.result, array=arraylimmamatcharray.result)
##   print(toparraymatch$genes==topgenelistSeqArray$genes)
##   for (gene in toparraymatch$genes){
## ##  mtscaled <- as.matrix(scale(mtcars))
## ##  heatmap(mtscaled, Colv=F, scale='none')
## ##  hc.rows <- hclust(dist(mtscaled))
## ##  plot(hc.rows)
## ##  hc.cols <- hclust(dist(t(mtscaled)))
## ##  heatmap(mtscaled, Colv=as.dendrogram(hc.cols), scale='none')
## ## # draw heatmap for first cluster
## ##  heatmap(mtscaled[cutree(hc.rows,k=2)==1,], Colv=as.dendrogram(hc.cols), scale='none')
## ## # draw heatmap for second cluster
## ## heatmap(mtscaled[cutree(hc.rows,k=2)==2,], Colv=as.dendrogram(hc.cols), scale='none')
## ## cutoff p.value 0.01, adjust p.value have no significant genes
##     seqmatch.gene <- (arraylimmamatchSeq.result[[gene]][arraylimmamatchSeq.result[[gene]]$P.Value<=0.01,])
##     arraymatch.gene <- (arraylimmamatcharray.result[[gene]][arraylimmamatcharray.result[[gene]]$P.Value<=0.01,])
##     print(intersect(seqmatch.gene, arraymatch.gene))
##     ## pdf(paste("limma", gene, "heatmap.pdf", sep = "", collapse = ""))
##     ## heatmap.2(as.matrix(seqexp), symm=TRUE, trace="none", symbreaks=TRUE, density.info="none",
##     ##           Rowv = FALSE, cexRow = 0.5,
##     ##           ## redgreen(75)
##     ##           col = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)),
##     ##           lmat=rbind(c(4, 4, 4), c(2,1,0), c(0,3,0)), lwid=c(1, 6, 1), lhei=c(2, 10, 1),
##     ##           scale="row")
##     ## dev.off()
match.difftable <- SeqArray.limma(expall=all, toparraymatch=head(toparray), topRNAseqmatch=head(topmatch),
                                  exp.match=exp.match, corgene=SeqExp.cor)

## cutoff 0.05 for non, upregulated and downregulated genes
## no need to normailized by exon lengths of genes
## topgenelist.norm <- topN.class()

## optionally, outlier of exon
pdf("../results/exonvsmutation.pdf")
mutation.norm <- exon.explore("../temp/hg19_exon.ref", Colon.Mutclass, 1e5)
dev.off()
write.table(mutation.norm$outlier, file="../results/outlier_exon_1e5.txt", quote=F, sep="\t", col.names = F)

## args <- commandArgs(trailingOnly = TRUE)
## print(args)
## main(args)
main <- function(){
  print();
}
