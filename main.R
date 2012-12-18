###################
## main for TCGA###
###################

## args <- commandArgs(trailingOnly = TRUE)
## print(args)
source("Array.R")
source("RNAseq.R")

## data.orig <- IO(args[1], args[2], args[3])
data.orig <- IO("../data/",
               "colon_cancer_TCGA_agilent_expression.xls",
               "colon_cancer_mutation_all.maf",
               "RNAseq_all_expression_Trim.txt")

## relieve memory
rm(data.orig$data1)

# patients mutation
mutation.patient <- substr(unique(data.orig$data2$Tumor_Sample_Barcode), 1, 15)
exp.patient <- colnames(data.orig$data1)
mutationexp.match <- match(mutation.patient, exp.patient)
mutationexp.patient <- exp.patient[mutationexp.match]
cat(exp.patient[mutationexp.match] == mutation.patient)
# write table of matched expression data
exp.match <- data.orig$data1[, mutationexp.match]
rownames(exp.match) <- data.orig$data1[,1]
write.table(exp.match, row.names=T, col.names=T, file="colon_exp_genes.txt", sep="\t", quote=F)
# Overall statistics of mutations genes in different patient's samples
Mutgenes.freq <- MutPartition(data.orig$data2)
Colon.Mutclass <- data.frame(genes=names(Mutgenes.freq),freq=as.numeric(Mutgenes.freq))
write.table(Colon.Mutclass, file="ColonMutationTop.txt", quote=F, sep="\t", col.names = T)
mutation.table <- Mut.Table(data.orig$data2, mutation.patient, Colon.Mutclass)
write.table(mutation.table, file="colon_somatic_mutations.txt", quote=F, sep="\t", col.names = T, row.names=T)

# match mutation genes symbol with expression genes symbol
length(match(unique(data.orig$data2$Hugo_Symbol),Colon.Mutclass$genes))

## exp.match.meancenter <- t(meancentering(as.matrix(exp.match)))
## For downstream genes analysis, use all expressed genes
## mutation.table denotes patient with mutation and array data
topgenelist <- topN.class(3, Colon.Mutclass, mutation.table)
## mutation.tableSeqArray denotes patient with mutation,array and RNASeq
rownames(mutation.table) <- mutation.table$genes
## Overlap RNASeq data with array and mutation
## RNAseq.matchmut <- mutation.patient %in% colnames(data.orig$data3)
RNAseq.matchmut <- match(mutation.patient[RNAseq.matchmut], colnames(data.orig$data3))
SeqArray.patient <- colnames(data.orig$data3)[RNAseq.matchmut]
RNAseq.matchexp <- match(exp.patient, colnames(data.orig$data3))
SeqArrayMut <- colnames(data.orig$data3)[RNAseq.matchmut] %in% mutationexp.patient
SeqArrayMut.patient <- colnames(data.orig$data3)[RNAseq.matchmut][SeqArrayMut]

## match for index
SeqExp <- rownames(data.orig$data3) %in% rownames(exp.match)
SeqExp.gene <- rownames(data.orig$data3)[SeqExp]
ExpSeq <- rownames(exp.match)  %in% rownames(data.orig$data3)
ExpSeq.gene <- rownames(exp.match)[ExpSeq]
SeqExpinter.gene <- intersect(rownames(data.orig$data3), rownames(exp.match))
SeqExpinter.patient <- intersect(colnames(data.orig$data3), colnames(exp.match))
## SeqExpinter.patientindex <- match(colnames(RNAseq), colnames(exp.match))
## SeqExpinter.patientindex <- SeqExpinter.patientindex[!isNA(SeqExpinter.patientindex)]
## SeqExpinter.patient <- colnames(RNAseq)[SeqExpinter.patientindex]
RNAseq.matcharray <- data.orig$data3[SeqExpinter.gene,SeqExpinter.patient]
## expression match with RNAseq
exp.matchRNAseq <- exp.match[SeqExpinter.gene, SeqExpinter.patient]
## exceptions for RNASeq annotated gene SLC35E2
exp.match["SLC35E2",]

match.expressionall <- cbind(RNAseq.matcharray, exp.matchRNAseq)
rm(RNAseq.matcharray, exp.matchRNAseq)

## get correlation of RNAseq and array overlapped genes
SeqExp.correlation <- function(expall="", cutoff=0.5){
  ## get log2 for RNAseq RPKM part,
  SeqExp.cor <- apply(expall, 1, function(x) cor(log2(x[1:(length(expall[1,])/2)]),
                                                            x[(length(expall[1,])/2+1):length(expall[1,])]))
  ## cutoff 0.5
  SeqExp.corgene <- names(sort((SeqExp.cor[SeqExp.cor>=cutoff]))) ##8423 for 0.5 correlation
  ## correlation heatmap
  SeqExp.corgene
}
SeqExp.cor <- SeqExp.correlation(match.expressionall, cutoff=0.5)

mutation.tableSeqArray <- mutation.table[, SeqExpinter.patient]
mutation.tableSeqArray <- cbind(genes=rownames(mutation.tableSeqArray), mutation.tableSeqArray)
topgenelistSeqArray <- topN.class(3, Colon.Mutclass, mutation.tableSeqArray)
colnames(match.expressionall[, 1:51]) == colnames(topgenelistSeqArray)[-c(1,2)]
## no need to normailized by exon lengths of genes
## topgenelist.norm <- topN.class()
patient.match.mutated <- apply(topgenelist, 1, function(x) {which(as.numeric(x[c(-1,-2)]) >= 1)})
names(patient.match.mutated) <- topgenelist$genes
patient.match.mutatedSeqArray <- apply(topgenelistSeqArray, 1, function(x) {which(as.numeric(x[c(-1,-2)]) >= 1)})
names(patient.match.mutatedSeqArray) <- topgenelistSeqArray$genes
## Read in RNASeq RPKM expression data
## call exp.patientclass to draw plot in batch
pdf("../results/Colon_Cancer_Scatter.pdf")
## classify patients' expression by mutation frequent genes and draw QQ-plot, heatmap
## only expression match with mutation
arraylimma.result <- exp.patientclass(patient.match.mutated, exp.match, topgenelist, topgenelist$genes, 1, "downstream", FALSE)
## RNA-seq expression match with mutation 1:51
arraylimmamatchSeq.result <- exp.patientclass(patient.match.mutatedSeqArray, match.expressionall[SeqExp.cor,1:51], topgenelistSeqArray, topgenelistSeqArray$genes, 1, "downstream", FALSE)
(arraylimmamatchSeq.result$APC[arraylimmamatchSeq.result$APC$P.Value<=0.01,])
length(arraylimma.result[arraylimma.result$APC$P.Value<=0.01])
dev.off()

## For feedback genes regulation analsysis, use overlapped genes
exp.match.withgene <- cbind(genes=rownames(exp.match), exp.match)
exp.match.overlap <- merge(Colon.Mutclass, exp.match.withgene, sort=F, all=F)
exp.match.overlaped <- exp.match.overlap[,-c(1,2)]
rownames(exp.match.overlaped) <- exp.match.overlap[,1]

pdf("../results/Colon_Cancer_Scatter_feedback.pdf")
exp.patientclass(patient.match.mutated, exp.match.overlaped, topgenelist, topgenelist$genes, 0.01, "feedback")
dev.off()


## outlier of exon
pdf("../results/exonvsmutation.pdf")
mutation.norm <- exon.explore("../temp/hg19_exon.ref", Colon.Mutclass, 1e5)
dev.off()
write.table(mutation.norm$outlier, file="../results/outlier_exon_1e5.txt", quote=F, sep="\t", col.names = F)

###############################
##### APC highlight analysis
###############################
## check APC genes mutation frequency 33 for colon cancer
sum(head(topgenelist[topgenelist[,1]=="APC",])[c(-1,-2)])

## checking data match
length(topgenelist[1,][as.numeric(topgenelist[1,]) >=1])
colnames(topgenelist)[c(-1,-2)] == colnames(exp.match)

# for APC
APC.mut <- exp.match[, patient.match.mutated$APC]
APC.non <- exp.match[, -patient.match.mutated$APC]

APC.mut <- exp.match.meancenter[, patient.match.mutated$APC]
APC.non <- exp.match.meancenter[, -patient.match.mutated$APC]

APC.reorder <- cbind(APC.non, APC.mut)
## for Heatmap.2 ColSideColors
patient.category <- rep("red", length(exp.match[1,]))
patient.category[patient.match.mutated$APC] <- "blue"

## common t.test and wilcoxon test
t.p.value <- apply(exp.match, 1,
                   function(x) {t.test(x[patient.match.mutated$APC], x[-patient.match.mutated$APC])$p.value})
w.p.value <- apply(exp.match, 1,
                   function(x) {wilcox.test(x[patient.match.mutated$APC], x[-patient.match.mutated$APC])$p.value})
length(w.p.value[w.p.value<0.01])
length(t.p.value[t.p.value<0.01])

## limma differential expressed genes extracted
control=c(rep(1,length(APC.non[1,])), rep(0,length(APC.mut[1,])))
treat=c(rep(0,length(APC.non[1,])), rep(1,length(APC.mut[1,])))
APC.result.gm <- rna.diff(APC.reorder,
                          control, treat,"gm")

## write out to results folder
write.table(sort(t.p.value), file="../results/ttest_colon.txt", quote=F, sep="\t", col.names = F)
write.table(sort(w.p.value), file="../results/wilcoxon_colon.txt", quote=F, sep="\t", col.names = F)
# Gaussian distribution QQ plot
qqnorm(t.p.value, col='red', main="APC unnormalized by exon"); qqline(t.p.value, col='blue')

# Uniform distribution QQ plot
t.p.value.unif <- runif(length(t.p.value), min=min(t.p.value), max=max(t.p.value))
w.p.value.unif <- runif(length(w.p.value), min=min(w.p.value), max=max(w.p.value))

qqplot(t.p.value.unif, t.p.value, col=2)
lines(loess(sort(t.p.value)~sort(t.p.value.unif)), col="Blue", lty=1, lwd=3)
qqplot(w.p.value.unif, w.p.value, col=2)
lines(loess(sort(w.p.value)~sort(w.p.value.unif)), col="Blue", lty=1, lwd=3)

## scatterplot for uniform t.p.value and observed t.p.value
library(car)
## library(rgl)
## scatter3d(sort(seq(length(t.p.value))), sort(t.p.value.unif), sort(t.p.value))
## smoothScatter(sort(t.p.value))
scatterplot(sort(t.p.value.unif),sort(t.p.value))

# common scatterplot
plot(sort(t.p.value.unif), sort(t.p.value),
     main="APC uniform p.value against observed p.value", pch=2, col="blue", lty=3)
# linear regression
abline(lm(sort(t.p.value)~sort(t.p.value.unif)), col="red", lty=1, lwd=4)
lines(loess(sort(w.p.value)~sort(w.p.value.unif)), col="yellow", lty=4, lwd=3)

