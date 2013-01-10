###################
## main for TCGA###
###################
################################################
## for mutation analysis
## driver genes detected
library(seqinr)
library(CorMut)
## limma has zscore calculations
## optional steps paradigm shift, kegggraph, cytoscape
library(iCluster)
library(som) # use som to replace iCluster knn
library(nnet)
require(CancerMutationAnalysis)
################################################
source("Array.R", verbose=F)
source("RNAseq.R")
################################################

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
  ## som has normalize methods for mean centering
  ## affy has C implementation of quantilte normalization
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
################################################
# 51  for all three data
# 53 for mutation and array
################################################
matchdata <- SeqArraymatch(datalist=data.orig)
attach(matchdata)
## array Combat, Seq quantile normalization of RPKM
## mutation type analysis
## SNP: point mutation
## separate by two alleles, normal and tumor and validated ones
## use all the validated genes
barplot(table(data.orig$data2$Variant_Type), ylim=c(0,15000), xlim=c(0,4), col=c(0,1,2))
barplot(table(data.orig$data2["Variant_Classification"]), ylim=c(0,8000), xlim=c(0,10), col=c(0:7))

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
par(mfrow = c(1,2))
barplot(mutct$cl$APC, main="APC mutation classification")
barplot(mutct$tp$APC, main="APC mutation type")
write.table(mutct$cl$APC, file="apc_type.txt", quote=F)

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
  all[index,]
}

corcutoff.out(indexall, cutoff="all", seq_array_mut)
corcutoff.out(index0.5, cutoff=0.5, seq_array_mut)
corcutoff.out(index0.3, cutoff=0.3, seq_array_mut)
corcutoff.out(index0.2, cutoff=0.2, seq_array_mut)
corcutoff.out(indexs0.5, cutoff="spearman0.5", seq_array_mut)
corcutoff.out(indexs0.2, cutoff="spearman0.2", seq_array_mut)
## after spearman correlation filter and normalization
all_norm <- corcutoff.out(indexs0.3, cutoff="spearman0.3", seq_array_mut)

## QC of outlier point
source("../code/nscore.R")
library(copa)
library(DriverNet) ## very similar to what we are doing
seq_outlier <- getPatientOutlierMatrix(all_norm[, 1:51]) ## DriverNet
length(which(seq_outlier == TRUE))  # outlier 18116
length(which(seq_outlier == FALSE)) # 637642
## test for NST
test1 <- nscore(as.vector(all_norm[,5]))
bak <- backtr(test1$nscore, test1)# tails="separate")
cor(all_norm[,5], bak)

## use NST for all
nst <- function(x){
  y <- nscore(as.vector(x))
  back <- backtr(y$nscore, y) ## back transform the nscored transform
  cor=cor(x, back)
}

all.brcor <- apply(all_norm, 2, function(x) nst(x)) ## back transform has high correlation
all.nst <- apply(all_norm, 2, function(x) nscore(x)$nscore) ## transformed matrix
rownames(all.nst) <- rownames(all_norm)
write.table(all.nst, file="nstransformed_norm_filterSpear0.3.txt", quote=F,sep="\t")
write.table(all_norm, file="nottransformed_norm_filterSpear0.3.txt", quote=F,sep="\t")
all.transcor <- apply(all_norm, 2, function(x) cor(nscore(x)$nscore,x)) ## back transform has high correlation

################################################
## step3, differential analysis of overlap
## mutation, arrays and Seq
## mutclass, data.orig$data2, mutable, filtered 156 high frequent mutated genes
seq_array_mut.index <- apply(mutable[filtered,], 1, function(x) {which(as.numeric(x) >= 1)})
## call exp.patientclass to draw plot in batch
## classify patients' expression by mutation frequent genes and draw QQ-plot, heatmap
########################################
## downstream  matchseq analysis
################################
## TODO mutation types classification involve
################################
## topmatch is the patient mutation data (matched with seq and array)
## adjust mutation to 51 matched patients
## all_norm, with correlation filtered 0.3 spearman correlation
## 1:51 seq, 52: 102 array all_norm
## cutoff only for P.value, No fold change cutoff yet
pdf("../results/seq_colon_seq_array_match_normfilter_downstream0.01.pdf")
seq_match_common <- exp.patientclass(seq_array_mut.index,
                                     all_norm[, 1:(length(all_norm[1,])/2)],
                                     filtered[1:5], cutoff=0.05, "downstream","matchseq")
dev.off()

pdf("../results/seq_colon_seq_array_match_normfilter_downstream0.01nst.pdf")
seq_match_commonnst <- exp.patientclass(seq_array_mut.index,
                                     all.nst[, 1:(length(all.nst[1,])/2)],
                                     filtered, cutoff=0.05, "downstream","matchseqnst")
dev.off()

pdf("../results/array_colon_seq_array_match_normfilter_downstream0.01.pdf")
array_match_common <- exp.patientclass(seq_array_mut.index,
                                       all_norm[, (length(all_norm[1,])/2+1):length(all_norm[1,])],
                                       filtered[1:5],cutoff= 0.05, "downstream","matcharray")
dev.off()

pdf("../results/array_colon_seq_array_match_normfilter_downstream0.01nst.pdf")
array_match_commonnst <- exp.patientclass(seq_array_mut.index,
                                       all.nst[, (length(all.nst[1,])/2+1):length(all.nst[1,])],
                                       filtered,cutoff= 0.05, "downstream","matcharraynst")
dev.off()

##############################
## To compare normalization and correlation filter effectiveness
## seq not normalization without correlation filter seq_array_mut differential expressed
## remove constant point
seqorig <- log2(seq_array_non[, c(1:51)] + 0.0000001)
remove <- apply(seqorig, 1, function(x) length(unique(x)) == 1) ## remove 51
seqremove <- seqorig[!remove,]
pdf("../results/seq_colon_seq_array_match_normfilter_downstream0.01origbefore.pdf")
seq_match_common_before <- exp.patientclass(seq_array_mut.index,
                                            seqremove,
                                            head(filtered, 2), cutoff=0.05, "downstream","matchseqnonfiler")
dev.off()

## array without seq normalization without correlation filter seq_array_mut differential expressed
pdf("../results/array_colon_seq_array_match_normfilter_downstream0.01origbefore.pdf")
## remove the array genes which overlaps with constant RNAseq points, to match in collective_gene
array_match_common_before<- exp.patientclass(seq_array_mut.index,
                                             seq_array_non[!remove, (length(all_norm[1,])/2+1):length(all_norm[1,])],
                                             head(filtered,2),cutoff= 0.05, "downstream","matcharray")
dev.off()
################################

collective_gene <- function(gene="", expr="", seqe="", arraye="", cutoff=0.05){
  ## using t.test value
  arrayt <- c()
  arrayt[1:length(expr[,1])] <- "none"
  print(length(arrayt))
  arrayt[arraye[[gene]]$tu] <- "up"
  arrayt[arraye[[gene]]$td] <- "down"
  print(length(arrayt))
  seqt <- c()
  print(names(head(seqe[[gene]]$tu)))
  seqt[1:length(expr[,1])] <- "none"
  seqt[seqe[[gene]]$tu] <- "up"
  seqt[seqe[[gene]]$td] <- "down"
  tu <- rownames(expr)[intersect(arraye[[gene]]$tu,seqe[[gene]]$tu)]
  write.table(tu, file=paste("t_collective_up", gene, ".xls", sep=""), quote=F, sep = "\t")
  td <- intersect(names(arraye[[gene]]$td),names(seqe[[gene]]$td))
  write.table(td, file=paste("t_collective_down", gene, ".xls",sep=""), quote=F, sep = "\t")
  write.table(data.frame(table(arrayt, seqt)), file=paste(gene,"0.05ttest", sep = "", collapse = ""), quote=F, sep="\t")
  arrayw <- c()
  arrayw[1:length(expr[,1])] <- "none"
  arrayw[arraye[[gene]]$wu] <- "up"
  arrayw[arraye[[gene]]$wd] <- "down"
  seqw <- c()
  seqw[1:length(expr[,1])] <- "none"
  seqw[seqe[[gene]]$wu] <- "up"
  seqw[seqe[[gene]]$wd] <- "down"
  wu <- rownames(expr)[intersect(arraye[[gene]]$wu,seqe[[gene]]$wu)]
  write.table(wu, file=paste("w_collective_up", gene, ".xls", sep=""), quote=F, sep = "\t")
  wd <- intersect(names(arraye[[gene]]$wd),names(seqe[[gene]]$wd))
  write.table(wd, file=paste("w_collective_down", gene, ".xls",sep=""), quote=F, sep = "\t")
  write.table(data.frame(table(arrayw, seqw)), file=paste(gene,"0.05wtest", sep = "", collapse = ""), quote=F, sep="\t")
  output <- list(t=table(arrayt, seqt), w=table(arrayw, seqw), tl=list(up=tu, down=td), wl=list(up=wu, down=wd))
  return(output)
}

# filtered
filter_overlap <- lapply(filtered, collective_gene, expr=all_norm, seqe=seq_match_common, arraye=array_match_common)
names(filter_overlap) <- filtered
## APC.non <- collective_gene("APC", expr=seqremove, seqe=seq_match_common_before, arraye = array_match_common_before)## seqe=seq_match_common, arraye=array_match_common)
save.image("all_seq_array.RData")  ## save workspace

setwd("../data/")
load("all_seq_array.RData")  ## recover



######################################
## limmma
## seq_match_limma <- limma.patientclass(patient.match.mutatedSeqArray,
##                                       log2(all[SeqExp.cor, 1:(length(all[1,])/2)]),
##                                       head(topmatch$genes), 1)
## arraylimmamatchSeq.result <- limma.patientclass(patient.match.mutatedSeqArray,
##                                                 log2(all[SeqExp.cor, 1:(length(all[1,])/2)]),
##                                                 head(topmatch$genes), 1)

########################################
## GSEA using expression diff
## For APC testing
##
source("../code/GOHyperGAll.R", verbose = F)
test.up <- filter_overlap$APC$wl$up
test.down <- filter_overlap$APC$wl$down
write.table(test.up, file="TEST_LIST_UP", quote=F)
write.table(test.down, file="TEST_LIST_DOWN", quote=F)
GOhyper2GSEA(myfile=c("", ""), type="all") ## input file names, output gmt genesets db


#################

GSEA.batch <- function(gene, data.type, analyze.type){
  source("../code/GSEA-P-R/GSEA.1.0.R", verbose=F, max.deparse.length=9999)
  o <- paste(gene, data.type, analyze.type,"_GSEA" sep = ".", collapse = "")
  system("mkdir ", o)
  GSEA(                                                  # Input/Output Files :-------------------------------------------
       input.ds =  "./Datasets/P53.gct",                 # Input gene expression Affy dataset file in RES or GCT format
       input.cls = "./Datasets/P53.cls",                 # Input class vector (phenotype) file in CLS format
       gs.db =     "../GSEA.1.0.R/GeneSetDatabases/C2.gmt",          # Gene set database in GMT format, Use curated and oncogenic datasets
       output.directory      = o,                        # Directory where to store output and results (default: "")
       ##  Program parameters :-------------------------------------------------------------------------------------------------------------------------
       doc.string            = paste(o,"analysis", sep = ""),        # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
       non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
       reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels"
       nperm                 = 1000,            # Number of random permutations (default: 1000)
       weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
       nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
       fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
       fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
       topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
       adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
 gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
       gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
       reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
       preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
       random.seed           = 760435,          # Random number generator seed. (default: 123456)
       perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
       fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
       replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
       save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
       OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
       use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
       )
## Overlap and leading gene subset assignment analysis of the GSEA results
  GSEA.Analyze.Sets(
    directory  = o,                                # Directory where to store output and results (default: "")
    topgs = 20,                                    # number of top scoring gene sets used for analysis
    height = 16,
    width = 16
    )
}


#### Regulator Potential
## TFRM test.R
## meta rank product analysis
metaRank <- function(rank1, rank2){
  ## meta analysis between different data type,
  ## e.g. RNASeq and Array
  sort(rank1*rank2)
}

#########################################
## GOstats, DAVID API in python and shell
## brainarray human annotation 16.0.0 version
source("../code/GOHyperGAll.R")
library(hgu133plus2hsentrezgcdf)       # entrez
library(hgu133plus2hsentrezg.db)       # chrom and genesname info
library(annotate); library(GOstats)
library(GO.db)
library(reactome.db)
goann <- as.list(GOTERM)
zz <- eapply(GOTERM, function(x) x@Ontology);
table(unlist(zz))
goterms <- unlist(eapply(GOTERM, function(x) x@Term))
goterms[grep("molecular_function", goterms)]

go_df <- data.frame(GOID=unlist(eapply(GOTERM, function(x) x@GOID)), Term=unlist(eapply(GOTERM, function(x) x@Term)),
                    Ont=unlist(eapply(GOTERM, function(x) x@Ontology)))
head(go_df)
#########################################

#########################################
## methylation(priority)

#########################################
mt <- read.table("methl27_all_expression.txt") ## methlation 27, 203 patients
## use methylation 45 is better
png("meth27_dist.png")
plot.ecdf(as.matrix(mt), col= "blue", lty=3)
abline(h=0.6, col="red" , lty=2)
dev.off()

mt.patient <- substr(gsub('.','-',colnames(mt),fixed=T), 1, 15)
length(intersect(mt.patient, colnames(mutable)))  ## 51 patients, all matched
length(intersect(mt.patient, expseq.patient)) ## 150 patients

## QC among all patients
## for beta value, no need to take log2
boxplot(mt, ylim=c(0, 1), col="gray",outline=F)   ## motif predict CpG methylation level, too biased, quantile????

## QC for the single array pass
pdf("../temp/QC_single_met.pdf")
apply(mt, 2, hist)
dev.off()
hist((mt[,1]))
library(minfi)  ## analyze methylation array
png("MDS_meth27.png", height=1200, width = 1200)
mdsPlot(as.matrix(mt), sampNames = colnames(mt))
dev.off()

## use UCSC cpg islands
library(minfi, verbose=F)
cpg <- read.delim("../data/hg19_cgi", header = F, col.names=LETTERS[1:4])

#########################################
## SNP(cnv(priority) and noncnv) and CN(waiting)
#########################################
cnv.patient <- c()

###########################################
## cluster and icluster and classification
###########################################
## 1. hierarchical
## 2. iCluster
#########################################

################################################
## archilles pipeline for synthetically lethality genes
################################################

#########################################
########################################
## feedback Analysis
## using the ones overlapped with important mutated genes
## Driver list, cosmic and archilles BRAF
## cutoff 0.05 for non, upregulated and downregulated genes
## TODO, BRAF mutation
###############################################
## 1. overlap with TCGA validated freq>= 5 mutated genes
###############################################
## Focus on important driver genes for correlation features
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
braf$genes ## 40th
matchDriver(driver)
interdriverw.seq <- intersect(as.vector(driver),names(seq_match_common$APC$w))
interdriverw.array <- intersect(as.vector(driver),names(array_match_common$APC$w))
interdriverw <- intersect(as.vector(driver), test.up)
seq_match_common$APC$w[interdriverw]

###############################################
## 2. overlap with archilles
###############################################
#### 1. BRAF
interdriverw <- intersect(as.vector(driver),names(seq_match_common$APC$w))

driver.diff <- intersect(driver, indexs0.3)
braf.diff <- intersect(braf$genes, indexs0.3)
## focus on braf
## correlation
seq.braf <- as.matrix(all_norm[braf.diff, 1:51]) ## seq
seq.braf <- t(scale(t(seq.braf)))
library(RColorBrewer)
png("BRAF_cor_exp.png", height = 1200, width = 1200)
par(mar=c(4,3,4,3), mai = c(1,1,0.5,1))
nf <- layout(matrix(c(2,1),byrow=TRUE), c(6, 6), c(1, 8), TRUE)
layout.show(nf)
x <- 10*(1:nrow(seq.braf))
y <- 10*(1:ncol(seq.braf))
image(x, y,seq.braf, col=colorRampPalette(brewer.pal(9, "Blues"))(100), axes=F, xlab = "", ylab = "",
      main = "BRAF essential genes expression")
axis(1, at=seq(min(x), max(x), length=length(x)), labels = rownames(seq.braf), las=2, ## side = 2,
     outer=F, tick = F, cex.axis=1)
axis(2, at=seq(min(y), max(y), length=length(y)), labels = colnames(seq.braf), las=1, ## side = 2,
     outer=F, tick = F, cex.axis=1)
seq.brafcor <- apply(seq.braf, 1, function(x) cor(x, seq.braf["BRAF", ]))
par(mar=c(5,3,5,3),mai=c(0.2,1,0.3,1))
seq.brafcor <- as.matrix(seq.brafcor)
x <- 10*(1:nrow(seq.brafcor))
y <- 10*(1:ncol(seq.brafcor))
image(x, y, seq.brafcor, col = colorRampPalette(brewer.pal(9, "Blues"))(100), xlab = "", ylab = "",
      axes=F, main="BRAF correlation with essential gene")
axis(2, at=seq(min(y), max(y), length=length(y)), labels = "BRAF", las=1, ## side = 2,
     outer=F, tick = F, cex.axis=1)
dev.off()
#### 2. others

###############################################
## 3. overlap with cosmic
###############################################


########################################
## optionally, outlier of exon
## no need to normailized by exon lengths of genes
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
