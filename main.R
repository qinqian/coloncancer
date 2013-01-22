###################
## SL predictor ##
###################
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

## qc of outlier point
source("../code/nscore.r")
library(copa)
library(drivernet) ## very similar to what we are doing
seq_outlier <- getpatientoutliermatrix(all_norm[, 1:51]) ## drivernet
length(which(seq_outlier == true))  # outlier 18116
length(which(seq_outlier == false)) # 637642
## test for nst
test1 <- nscore(as.vector(all_norm[,5]))
bak <- backtr(test1$nscore, test1)# tails="separate")
cor(all_norm[,5], bak)

## use nst for all
nst <- function(x){
  y <- nscore(as.vector(x))
  back <- backtr(y$nscore, y) ## back transform the nscored transform
  cor=cor(x, back)
}

all.brcor <- apply(all_norm, 2, function(x) nst(x)) ## back transform has high correlation
all.nst <- apply(all_norm, 2, function(x) nscore(x)$nscore) ## transformed matrix
rownames(all.nst) <- rownames(all_norm)
write.table(all.nst, file="nstransformed_norm_filterspear0.3.txt", quote=f,sep="\t")
write.table(all_norm, file="nottransformed_norm_filterspear0.3.txt", quote=f,sep="\t")
all.transcor <- apply(all_norm, 2, function(x) cor(nscore(x)$nscore,x)) ## back transform has high correlation

################################################
## step3, differential analysis of overlap
## mutation, arrays and seq
## mutclass, data.orig$data2, mutable, filtered 156 high frequent mutated genes
seq_array_mut.index <- apply(mutable[filtered,], 1, function(x) {which(as.numeric(x) >= 1)})
## call exp.patientclass to draw plot in batch
## classify patients' expression by mutation frequent genes and draw qq-plot, heatmap
########################################
## downstream  matchseq analysis
################################
## todo mutation types classification involve
################################
## topmatch is the patient mutation data (matched with seq and array)
## adjust mutation to 51 matched patients
## all_norm, with correlation filtered 0.3 spearman correlation
## 1:51 seq, 52: 102 array all_norm
## cutoff only for p.value, no fold change cutoff yet

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
## to compare normalization and correlation filter effectiveness
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
## remove the array genes which overlaps with constant rnaseq points, to match in collective_gene
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
  write.table(tu, file=paste("t_collective_up", gene, ".xls", sep=""), quote=f, sep = "\t")
  td <- intersect(names(arraye[[gene]]$td),names(seqe[[gene]]$td))
  write.table(td, file=paste("t_collective_down", gene, ".xls",sep=""), quote=f, sep = "\t")
  write.table(data.frame(table(arrayt, seqt)), file=paste(gene,"0.05ttest", sep = "", collapse = ""), quote=f, sep="\t")
  arrayw <- c()
  arrayw[1:length(expr[,1])] <- "none"
  arrayw[arraye[[gene]]$wu] <- "up"
  arrayw[arraye[[gene]]$wd] <- "down"
  seqw <- c()
  seqw[1:length(expr[,1])] <- "none"
  seqw[seqe[[gene]]$wu] <- "up"
  seqw[seqe[[gene]]$wd] <- "down"
  wu <- rownames(expr)[intersect(arraye[[gene]]$wu,seqe[[gene]]$wu)]
  write.table(wu, file=paste("w_collective_up", gene, ".xls", sep=""), quote=f, sep = "\t")
  wd <- intersect(names(arraye[[gene]]$wd),names(seqe[[gene]]$wd))
  write.table(wd, file=paste("w_collective_down", gene, ".xls",sep=""), quote=f, sep = "\t")
  write.table(data.frame(table(arrayw, seqw)), file=paste(gene,"0.05wtest", sep = "", collapse = ""), quote=f, sep="\t")
  output <- list(t=table(arrayt, seqt), w=table(arrayw, seqw), tl=list(up=tu, down=td), wl=list(up=wu, down=wd))
  return(output)
}

# filtered
filter_overlap <- lapply(filtered, collective_gene, expr=all_norm, seqe=seq_match_common, arraye=array_match_common)
names(filter_overlap) <- filtered
## apc.non <- collective_gene("apc", expr=seqremove, seqe=seq_match_common_before, arraye = array_match_common_before)## seqe=seq_match_common, arraye=array_match_common)

######################################
## limmma
## seq_match_limma <- limma.patientclass(patient.match.mutatedseqarray,
##                                       log2(all[seqexp.cor, 1:(length(all[1,])/2)]),
##                                       head(topmatch$genes), 1)
## arraylimmamatchseq.result <- limma.patientclass(patient.match.mutatedseqarray,
##                                                 log2(all[seqexp.cor, 1:(length(all[1,])/2)]),
##                                                 head(topmatch$genes), 1)

########################################
## gsea using expression diff
## for apc testing
##
source("../code/gohypergall.r", verbose = f)
test.up <- filter_overlap$apc$wl$up
test.down <- filter_overlap$apc$wl$down
write.table(test.up, file="test_list_up", quote=f)
write.table(test.down, file="test_list_down", quote=f)
gohyper2gsea(myfile=c("", ""), type="all") ## input file names, output gmt genesets db, may be use default gmt

## using readLines and cat to creat gcl or rnk, and pheno

#################
### may cross enrichment, e.g, make a gmt file for interest gene downstream, then use other genes' downstream to
### enrich in this one.
gsea.batch <- function(gene, data.type, analyze.type){
  source("../code/gsea-p-r/gsea.1.0.r", verbose=f, max.deparse.length=9999)
  o <- paste(gene, data.type, analyze.type,"_gsea" sep = ".", collapse = "")
  system("mkdir ", o)
  gsea(                                                  # input/output files :-------------------------------------------
       input.ds =  "./datasets/p53.gct",                 # input gene expression affy dataset file in res or gct format
       input.cls = "./datasets/p53.cls",                 # input class vector (phenotype) file in cls format
       gs.db =     "../gsea.1.0.r/genesetdatabases/c2.gmt",          # gene set database in gmt format, use curated and oncogenic datasets
       output.directory      = o,                        # directory where to store output and results (default: "")
       ##  program parameters :-------------------------------------------------------------------------------------------------------------------------
       doc.string            = paste(o,"analysis", sep = ""),        # documentation string used as a prefix to name result files (default: "gsea.analysis")
       non.interactive.run   = f,               # run in interactive (i.e. r gui) or batch (r command line) mode (default: f)
       reshuffling.type      = "sample.labels", # type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels"
       nperm                 = 1000,            # number of random permutations (default: 1000)
       weighted.score.type   =  1,              # enrichment correlation-based weighting: 0=no weight (ks), 1= weigthed, 2 = over-weigthed (default: 1)
       nom.p.val.threshold   = -1,              # significance threshold for nominal p-vals for gene sets (default: -1, no thres)
       fwer.p.val.threshold  = -1,              # significance threshold for fwer p-vals for gene sets (default: -1, no thres)
       fdr.q.val.threshold   = 0.25,            # significance threshold for fdr q-vals for gene sets (default: 0.25)
       topgs                 = 20,              # besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
       adjust.fdr.q.val      = f,               # adjust the fdr q-vals (default: f)
 gs.size.threshold.min = 15,              # minimum size (in genes) for database gene sets to be considered (default: 25)
       gs.size.threshold.max = 500,             # maximum size (in genes) for database gene sets to be considered (default: 500)
       reverse.sign          = f,               # reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: f)
       preproc.type          = 0,               # preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
       random.seed           = 760435,          # random number generator seed. (default: 123456)
       perm.type             = 0,               # for experts only. permutation type: 0 = unbalanced, 1 = balanced (default: 0)
       fraction              = 1.0,             # for experts only. subsampling fraction. set to 1.0 (no resampling) (default: 1.0)
       replace               = f,               # for experts only, resampling mode (replacement or not replacement) (default: f)
       save.intermediate.results = f,           # for experts only, save intermediate results (e.g. matrix of random perm. scores) (default: f)
       old.gsea              = f,               # use original (old) version of gsea (default: f)
       use.fast.enrichment.routine = t          # use faster routine to compute enrichment for random permutations (default: t)
       )
## overlap and leading gene subset assignment analysis of the gsea results
  gsea.analyze.sets(
    directory  = o,                                # directory where to store output and results (default: "")
    topgs = 20,                                    # number of top scoring gene sets used for analysis
    height = 16,
    width = 16
    )
}

#### regulator potential
## tfrm test.r
## meta rank product analysis
metarank <- function(rank1, rank2){
  ## meta analysis between different data type,
  ## e.g. rnaseq and array
  sort(rank1*rank2)
}

#############################################
## gostats, david api in python and shell,GO pipeline
## brainarray human annotation 16.0.0 version
source("../code/gohypergall.r")
library(hgu133plus2hsentrezgcdf)       # entrez
library(hgu133plus2hsentrezg.db)       # chrom and genesname info
library(annotate); library(gostats)
library(go.db)
library(reactome.db)
goann <- as.list(goterm)
zz <- eapply(goterm, function(x) x@ontology);
table(unlist(zz))
goterms <- unlist(eapply(goterm, function(x) x@term))
goterms[grep("molecular_function", goterms)]
go_df <- data.frame(goid=unlist(eapply(goterm, function(x) x@goid)), term=unlist(eapply(goterm, function(x) x@term)),
                    ont=unlist(eapply(goterm, function(x) x@ontology)))
head(go_df)
#########################################

#########################################
## methylation(priority)
#########################################
## calculate C proportion in -25k to 25k(promoter region)
## separate patients by genes' promoter methylation level, high and low methylation
## overlap with CpG, CpG gene and nonCpG gene, and CpG high methylation, CpG low methylation
mt <- read.table("methl27_all_expression.txt", header=T, row.names=1) ## methlation 27, 203 patients
## mt45 <- read.table("./methl450_all_expression.txt", header=T, row.names=1) ## too large
mt.genes <- rownames(mt)
mt.genes <- gsub("(^.+)_(\\w+)_(\\d+)", "\\1", mt.genes) ## to remove position info
mt <- cbind(mt, genes=as.factor(mt.genes))
require(plyr)
test <- head(mt, 1000) ## test program
## using subset to focus on more group measurement, to divide and conquer

mt.bygenes <- list()
genes <- unique(mt.genes)

for (i in seq(along.with=genes)){
  mt.bygenes[[genes[i]]] <- apply(subset(mt, grepl(paste("^",genes[i],"_", sep="",collapse=""), rownames(mt)), select=1:ncol(mt))[, -ncol(mt)],
                                  2, mean)
}

mt.bygenes <- t(data.frame(mt.bygenes))
write.table(mt.bygenes, file="methy27k_by_genes.txt", sep="\t", quote=F)
mt.bygenes = read.table("methy27k_by_genes.txt", header=T, row.names=1)
length(intersect(rownames(mt.bygenes), rownames(mutable)))
length(intersect(rownames(mt.bygenes), rownames(all_norm)))


## what about hCG_1817306 and intergenic methylation regions
## using plyr, better used for less than 3 group standard
## mt.bygene <- ddply(test, .(genes), summarize, "TCGA.AA.A01P.11_66"=mean(TCGA.AA.A01P.11_66, na.rm=T))
## mt.bygene1 <- ddply(test, .(genes), summarize, "TCGA.AA.3673.01_36"=mean(TCGA.AA.3673.01_36, na.rm=T))
## mt.bygene2 <- ddply(test, .(genes), summarize, "TCGA.AA.3696.01_36"=mean(TCGA.AA.3696.01_36, na.rm=T))

## testmerge <- merge(mt.bygene, mt.bygene1)
## tesetmerge2 <- merge(testmerge, mt.bygene2)

## ok <- list()
## okk <- colnames(test[-length(colnames(test))])

## for (i in seq(along.with=okk)){
##   ok[[1]] <-  ddply(mt, .(genes), summarize, TCGA.AA.3696.01_36=mean(TCGA.AA.3696.01_36, na.rm=T))
## }

## use methylation 45 is better
png("meth27_dist.png")
plot.ecdf(as.matrix(mt), col= "blue", lty=3)
abline(h=0.6, col="red" , lty=2)
dev.off()
mt.patient <- substr(gsub('.','-',colnames(mt), fixed=T), 1, 15)
colnames(mt) <- mt.patient

length(intersect(mt.patient, colnames(mutable)))  ## 51 patients, all matched
length(intersect(mt.patient, expseq.patient)) ## 150 patients
## classify

## qc among all patients
## for beta value, no need to take log2
boxplot(mt, ylim=c(0, 1), col="gray",outline=f)   ## motif predict cpg methylation level, too biased, quantile????

## qc for the single array pass
pdf("../temp/qc_single_met.pdf")
apply(mt, 2, hist)
dev.off()
hist((mt[,1]))
library(minfi)  ## analyze methylation array
png("mds_meth27.png", height=1200, width = 1200)
mdsplot(as.matrix(mt), sampnames = colnames(mt))
dev.off()

## use ucsc cpg islands
library(minfi, verbose=f)
cpg <- read.delim("../data/hg19_cgi", header = f, col.names=letters[1:4])

meth.diff <- function(
  x="data",
  method=function(x) {print}, ## fisher.test or cancer/normal
  cutoff="")
{
  if (is.matrix(x) ||
      is.data.frame(x))
    print(1)
  else
    {cat(class(x));cat("try other format\n")}
}


#########################################
## snp(cnv(priority) and noncnv) and cn(waiting)
#########################################
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(LearnBayes)
library(rsgcc)  ## need GTK2+
snpcnv <- read.table("snp_cnv_1bp_final.txt", header = T, row.names = 1)
snpcnv_m<- as.matrix(snpcnv)
colnames(snpcnv_m) <- substr(gsub(".","-", colnames(snpcnv_m), fixed=T), 1, 15)
boxplot(snpcnv_m[, colnames(mutable)])
## y denotes seg.mean for each gene, x denotes gene
plot(sort(snpcnv_m[,444]), pch=".", col="red")
cn.exp.genes <- (intersect(rownames(snpcnv_m), rownames(all_norm)))
cn.exp.patient <- intersect(colnames(mutable), colnames(snpcnv_m))
cn.exp <- snpcnv_m[cn.exp.genes, cn.exp.patient]
### copy number filtering amplification > 1.5
seq.cn <- all_norm[cn.exp.genes, 1:51]
array.cn <- all_norm[cn.exp.genes, 52:102]
cn.seq <- cbind(seq.cn, cn.exp)
cn.array <- cbind(array.cn, cn.exp)
cn.seq.MIC <- apply(cn.seq, 1, function(x) MIC(x[1:51], x[52:102]))
cn.seq.MICv <- unlist(cn.seq.MIC)
write.table(cn.seq.MICv, file="cn_seq_MIC.txt", sep="\t", quote=F)
cn.seq.MICv = read.table("cn_seq_MIC.txt")
seq.mean <- apply(cn.seq, 1, mean)
cn.mean <- apply(cn.exp, 1, mean)
cn.miss<-grep("SCGB1C1", names(seq.mean))  ## some "NA" missing values
cn.seq.MICdf <- cbind(cn_mic=cn.seq.MICv, seq=seq.mean[-cn.miss], cn=cn.mean[-cn.miss])
ggplot2::qplot( cn.mean[-cn.miss],seq.mean[-cn.miss], main="copy number vs RNAseq expression in colon cancer")
p <- qplot(cn.mean[-cn.miss],seq.mean[-cn.miss], colour=cn.seq.MICv$x, main="copy number vs RNAseq expression in colon cancer")
p
qplot(cn.seq.MICv, binwidth=0.05)
cn.array.MIC <- apply(cn.array, 1, function(x) MIC(x[1:51], x[52:102]))

###########################################
## cluster and icluster and classification
###########################################
## 1. hierarchical
## 2. icluster
#########################################

################################################
## archilles pipeline for synthetically lethality genes
################################################

#########################################
########################################
## feedback analysis
## using the ones overlapped with important mutated genes
## driver list, cosmic and archilles braf
## cutoff 0.05 for non, upregulated and downregulated genes
## todo, braf mutation
###############################################
## 1. overlap with tcga validated freq>= 5 mutated genes
###############################################
## focus on important driver genes for correlation features
## cancer driver
## archilles braf
braf <- read.delim("../../archilles/archilles_braf_han.txt", header=f)
names(braf) <- "genes"
driver <- read.xls("../../archilles/drivers, oct 16 2012.xlsx")
driver <- driver$hugo

matchdriver <- function(genes,..){
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
matchdriver(braf$genes)
braf$genes ## 40th
matchdriver(driver)
interdriverw.seq <- intersect(as.vector(driver),names(seq_match_common$apc$w))
interdriverw.array <- intersect(as.vector(driver),names(array_match_common$apc$w))
interdriverw <- intersect(as.vector(driver), test.up)
seq_match_common$apc$w[interdriverw]

###############################################
## 2. overlap with archilles
###############################################
#### 1. braf
interdriverw <- intersect(as.vector(driver),names(seq_match_common$apc$w))

driver.diff <- intersect(driver, indexs0.3)
braf.diff <- intersect(braf$genes, indexs0.3)
## focus on braf
## correlation
seq.braf <- as.matrix(all_norm[braf.diff, 1:51]) ## seq
seq.braf <- t(scale(t(seq.braf)))
library(RColorBrewer)
brafmut <-as.vector(mutable["BRAF",])
brafmut[which(brafmut>=1)] <- "mut"
brafmut[which(brafmut==0)] <- "non"

png("BRAF_cor_exp.png", height = 1200, width = 1200)
par(mar=c(4,3,4,3), mai = c(1,1,0.5,1))
nf <- layout(matrix(c(2,1),byrow=TRUE), c(6, 6), c(1, 8), TRUE)
layout.show(nf)
x <- 10*(1:nrow(seq.braf))
y <- 10*(1:ncol(seq.braf))
image(x, y,seq.braf, col=rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100)), axes=F, xlab = "", ylab = "",
      main = "BRAF essential genes expression")
axis(1, at=seq(min(x), max(x), length=length(x)), labels = rownames(seq.braf), las=2, ## side = 2,
     outer=F, tick = F, cex.axis=1)
axis(2, at=seq(min(y), max(y), length=length(y)), labels = brafmut, las=1, ## side = 2,
     outer=F, tick = F, cex.axis=1)
par(mar=c(5,3,5,3),mai=c(0.2,1,0.3,1))
seq.brafcor <- as.matrix(seq.brafcor)
x <- 10*(1:nrow(seq.brafcor))
y <- 10*(1:ncol(seq.brafcor))
image(x, y, seq.brafcor, col = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100)), xlab = "", ylab = "",
      axes=F, main="BRAF correlation with essential gene")
axis(2, at=seq(min(y), max(y), length=length(y)), labels = "BRAF", las=1, ## side = 2,
     outer=F, tick = F, cex.axis=1)
dev.off()

seq.brafp <- apply(seq.braf, 1, function(x) cor(x, as.numeric(all_norm["BRAF", 1:51]), method="pearson"))
seq.brafcor <- apply(seq.braf, 1, function(x) cor(x, as.numeric(all_norm["BRAF", 1:51]), method="spearman"))
seq.brafMIC <- apply(seq.braf, 1, function(x) MIC(x, as.numeric(all_norm["BRAF", 1:51])))
seq.brafdist <- apply(seq.braf, 1, function(x) distance(x, as.numeric(all_norm["BRAF", 1:51])))
seq.brafkendall <- apply(seq.braf, 1, function(x) cor(x, as.numeric(all_norm["BRAF", 1:51]), method="kendall"))

qplot(seq.brafcor, seq.brafMIC)
scatterplot(as.vector(seq.brafcor), seq.brafMIC)

## par(bg="black")
## using RNA sequence
plot(seq.brafcor, col="red", lty=1, xlab="genes", ylab="correlation value", type="b", ylim=c(-0.8,0.7))
lines(seq.brafMIC, col="blue", lty=2)
lines(seq.brafdist, col="black", lty=3)
lines(seq.brafp, col="purple", lty=4)
lines(seq.brafkendall, col="green", lty=5)
title("BRAF essential genes' correlation methods comparison")
legend("bottomright", paste("cor:", c("spearman", "MIC", "distance", "pearson", "kendall")),
       inset=0.01, lty=1:5, col=c("red", "blue", "black", "purple", "green"), border="black", merge=T)


#########################################
## micRNA GA, Hiseq(wait)
#########################################
miRNA <- read.table("../data/miRNA_all_expressionGA.txt", header=T, row.names=1)
miRNArecords <- read.xls("~/Desktop/miRecords_version3.xls", sheet=1)
miRNAmatch <- subset(miRNArecords, Target.gene_species_scientific=="Homo sapiens", select=c(2,4,8))

miRNAmatch[,3] <- gsub('*','', miRNAmatch[,3],fixed=T)
miRNAmatch[,3] <- gsub('[','', miRNAmatch[,3],fixed=T)
miRNAmatch[,3] <- gsub(']','', miRNAmatch[,3],fixed=T)
miRNAmatch[,3] <- gsub('has','hsa', miRNAmatch[,3],fixed=T)
miRNAmatch[,3] <- gsub('R','r', miRNAmatch[,3],fixed=T)

mitest <- list()
for (i in seq(along.with=rownames(miRNA))){
  mitest[[i]] <- grep(rownames(miRNA)[i], miRNAmatch[,3], perl=T)
}

length(intersect(miRNAmatch[,3], rownames(miRNA)))

colnames(miRNA) <- substr(gsub('.','-',colnames(miRNA),fixed=T), 1, 15)
miRNA <- as.matrix(miRNA)
mi.match <- intersect(colnames(miRNA), colnames(mutable))
length(intersect(colnames(miRNA), colnames(mutable))) ## 44 patients

## micRNA array filter standard, one miRNA expression >0.5 account for 50% among patients
## inter quantile more than 1 when taken log2
miRNA <- miRNA[apply(miRNA > 100, 1, sum)/length(miRNA[1,])>0.6 & apply(log2(miRNA), 1, IQR) > 1.5, ]
mut.p <- colnames(mutable)[seq_array_mut.index[[ATM]]]
## mut.p %in%
## boxplot(y, x, log2(miRNA[, intersect(mut.p, colnames(miRNA))]), axes=F) ## all
miRNA <- miRNA[, mi.match]
colnames(miRNA) <- c("mut", "non")

scale.miRNA <- t(scale(t(log2(miRNA))))
miRNA_v <- as.vector(t(scale.miRNA))
miRNA <- data.frame(expand.grid(y=colnames(miRNA), x=rownames(miRNA)),v=miRNA_v)
pdf("miRNA_filter.pdf")
ggplot(miRNA, aes(y, x, fill = v, label = sprintf("%.1f", v)), xlab="", ylab="")+
  geom_tile() + geom_text() +
  scale_fill_gradient2(low = "blue", high = "red")
dev.off()

col.aes <- brewer.pal(9,"RdYlBu")
heatmap.2(scale.miRNA, trace="none", col=rev(col.aes), keysize=0.8, margins=c(10,10)) ## use margins to show fonts

par(mar=c(1, 3, 0.2, 0.2), oma=c(0.5,2,1,2), mai=c(1,2,1,2))
boxplot(log2(miRNA), outline=F, axes=F, frame=T, horizontal = T,
        at=seq(1-0.2, ncol(miRNA)+0.3, length=ncol(miRNA)), xlim=1*c(1-0.5, ncol(miRNA)+0.7)) ## ajust at=, xlim=, to ajust box widths, <1 better
p <- ggplot(miRNA[1:100,], aes(y, x, fill=y))
p+geom_boxplot()
## par(las=1) ## all horizontal
axis(2, at=seq(1-0.2, ncol(miRNA)+0.3,length=ncol(miRNA)), labels = colnames(miRNA), las=2,
     tick = T, cex.axis=0.8, outer=F)
par(mar=c(12, 0.3, 1, 1), oma=c(0.5,2,1,2), mai=c(3,1,1,0.5))
boxplot(log2(miRNA), outline=F, axes=F, frame=T, horizontal = F,
        at=seq(1-0.2, ncol(miRNA)+0.3, length=ncol(miRNA)), xlim=1*c(1-0.5, ncol(miRNA)+0.7)) ## ajust at=, xlim=, to ajust box widths, <1 better
## par(las=1) ## all horizontal
axis(1, at=1*(1:ncol(miRNA) - 0.2), labels = colnames(miRNA), las=2,
     tick = T, cex.axis=0.8, outer=F, font=4) ## font and las to define character
text(1,3, labels="test",srt=60)  ## for font angle =======
boxplot(log2(miRNA))

#############################################
## protein array, level 2, need normalization
############################################
protein <- read.table("../data/colon_protein_level2.txt", header=T,row.names=1)
colnames(protein) <- substr(gsub('.','-',colnames(protein),fixed=T), 1, 15)
protein.match <- (intersect(colnames(protein), colnames(mutable))) ## 24 patients
## BRAF classified
braf.mut <- colnames(mutable)[which(mutable["BRAF" ,] >= 1)]
## remove order
protein <- protein[-1, ]
protein <- as.matrix(protein)
ATM <- grep("^ATM.*", rownames(protein), perl=T)
protein[ATM,]

QC(protein, ggplot=FALSE)

## integrate different correlation method
method=c("MIC", "distance", "spearman", "pearman", "kendall",
  "liquid", "entropy"))

distance <- function(x, y){
  require(energy)
  dcov(x, y)
}

MIC <- function(x, y) {
      minedata <- rbind(x, y)
      source("MINE.R")
      rMINE(minedata, "matrix", 0)
      minecor <- read.csv("matrix,mv=0,cv=0.0,B=n^0.6,Results.csv")[,3]
      return(minecor)
}

classify <-
  ## classfify data by mutation data
  function(standard="", data1, data2, ...){
}

QC <- function(
  input = "", ## input data type
  class = "", ## class to denote mutation state
  ggplot = TRUE) ## ggplot or hist
{
  if (is.matrix(input)) {
    ## convert to data.frame to import ggplot
    if (ggplot){
      input_v <- as.vector(input)
      input_d <- data.frame(expand.grid(y=colnames(input), x=rownames(input)),v=input_v)
      p <- ggplot(input_d, aes(y, x, fill=y))
      p+geom_boxplot()
    }
    else {
      boxplot(input, outline=F)
    }
  }
}

param <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  return(args)
}

sl <- function(
  inputdir = "", ## input data directory
  outputdir = "", ## output
  )
{
  ## main function -- sl for synthetic gene prediction
  args <- param()
  print(args)
  library(RColorBrewer)
  library(seqinr)
  library(CorMut)
  ## limma has zscore calculations
  ## optional steps paradigm shift, kegggraph, cytoscape
  library(iCluster)
  library(som) # use som to replace iCluster knn
  library(nnet)
  require(CancerMutationAnalysis)
  library(ggplot2)
  library(gplots)
  library(car)
  library(rgl)
  ################################################
  source("Array.R", verbose=F)
  source("RNAseq.R")
  ################################################
  ## data.orig <- IO(args[1], args[2], args[3])
  data.orig <- IO("../data/",
                  "colon_cancer_TCGA_agilent_expression.xls",
                  "colon_cancer_mutation_all.maf",
                  "RNAseq_all_expression_Trim.txt")
}
