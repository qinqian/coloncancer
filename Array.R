## Author: Qin Qian
## Time-stamp: < modified by qinqianhappy :2012-12-17 21:16:56 >
## TCGA process
## Usage: Downstream genes and feedback loop analysis on exp and mutation data

## Focusing on the mutated genes with expression data
## filter expression data overlapped with only mutated genes now
## TODO: modify to use location info to analyze gene expression
# limma for arrays
library(ggplot2)
library(gplots)
require(limma)
library(RColorBrewer)
library(car)
## gene ontology
library(GOstats)

## optional steps paradigm shift, kegggraph, cytoscape
library(copa)
library(iCluster)
library(som) # use som to replace iCluster knn
require(CancerMutationAnalysis)

IO <- function(workdirectory, datatype1, datatype2, datatype3,...){
  # datatype2, e.g. mutations
  # datatype1, e.g. expression
  setwd(workdirectory)
  data1 <- read.table(datatype1, header=T, sep="\t", allowEscapes = FALSE)
  data2 <- read.table(datatype2, header=T, sep="\t", allowEscapes = FALSE)
  data3 <- read.table(datatype3, header=T, row.names=1, allowEscapes = F)
  colnames(data3) <- substr(gsub('.','-',colnames(data3),fixed=T), 1, 15)
  colnames(data1) <- substr(gsub('.','-',colnames(data1),fixed=T), 1, 15)
  data <- list(data1=data1, data2=data2, data3=data3)
  data
}

CaMP <- function(n, PROB, RANK) {
  ## function from Census coding sequencing, Science, 2006
  ## n stands for gene number evaluated
  ## PROB stands for the probability of having the observed mutation number
  ## RANK represents its numerical positions
  value <- -log10(n * PROB/ RANK)
  value
  }

MutPartition <- function (data, top=10, ...) {
  ## partition Mutation data according to patients and samples
  sort(table(data[,1]), decreasing = T)
}

## rna differential expression 
rna.diff <- function(expr, control, treat, method, pcutoff, fcutoff){
  if (method=="foldT"){
    # fold change and T-test
    foldchange=apply(expr, 1, function(x) mean(x[treat])-mean(x[control]))
    T.p.value=apply(expr, 1,
                    function(x) t.test(x[treat], x[control], var.equal=T)$p.value)
    #fdr=p.adjust(T.p.value, method="BH") # BH, Bonferroni, fdr
    fdr = T.p.value
    # genes
    genes.up = expr[which(fdr<0.05 & foldchange>0)]
    genes.down = expr[which(fdr<0.05 & foldchange<0)]
    genes.id = c(which(fdr<0.05 & foldchange>0), which(fdr<0.05 & foldchange<0))
    return(expr[genes.id,])
  }
  # limma group means
  if (method=="gm"){
    ## gm.design<-model.matrix(~ 0+factor(c(1,1,2,2)))
    gm.design<-cbind(control, treat)
    print(gm.design)
    colnames(gm.design) <- c("control","treat")
    #gm.design = cbind(control = control, treat = treat)
    print(gm.design)
    gm.fit = lmFit(expr, gm.design)
    print(gm.fit)
    gm.matrix = makeContrasts(CvsT = control-treat, levels=gm.design)
    gm.fit = contrasts.fit(gm.fit, gm.matrix)
    gm.fit = eBayes(gm.fit)
    # coef=2 for between-group difference for two group comparison, coef=1 for intercept
    # coef=1,2,3 for multiple group between-group differences
    diff.expr <- topTable(gm.fit,  number=length(expr[,1]), adjust.method="BH",
                          p.value=pcutoff,coef='CvsT')
    return(diff.expr)
  }
  # limma parameters
}

# separate patient with mutated genes
Mut.Table <- function(mutdata, patients, mutfreq,...){
  patient.Mutfreq.merge <- mutfreq
  for(patientid in seq(along.with=patients)){
    mutation.patient2genes <- mutdata[grep(paste("^", patients[patientid],"*",sep=""), mutdata$Tumor_Sample_Barcode, perl=T), c(1,16)]
    ## patient.mutationfreq <- table(unique(mutation.patient2genes[,1])) # genes mutation freq
    patient.mutationfreq <- table(mutation.patient2genes[,1]) # genes mutation freq

    ## ???? for each patient, overview gene mutation
    mutation.patient.freq <- data.frame(genes=names(patient.mutationfreq), patient=as.vector(patient.mutationfreq))
    mutation.patient.names <- apply(matrix(colnames(mutation.patient.freq)), 1, sub, pattern="patient", replacement=patients[patientid])
    colnames(mutation.patient.freq) <- mutation.patient.names
    patient.Mutfreq.merge <- merge(patient.Mutfreq.merge, mutation.patient.freq)
  }
  patient.Mutfreq.merge
}

## use top n mutations genes for patient classification
topN.class <- function(N, mutation.freq, mutation.table, cutoff=0.05){
  topN <- head(mutation.freq, N)
  print(topN)
  top.genes.patient <- merge(topN, mutation.table, all=F, sort=F)
  top.genes.patient
  #head(top.genes.patient[top.genes.patient$genes=='APC', ])
}

exp.patientclass <- function(patient.match.mutated="",exp.match="", topgenelist="", mut.genes="", cutoff="", type="", draw=""){
  limma.result.gm <- list()
  for(gene in mut.genes) {
      cat(gene);
      mut <- exp.match[, patient.match.mutated[[gene]]]
      non <- exp.match[, -patient.match.mutated[[gene]]]
      exp.reorder <- cbind(mut, non)
      ## design matrix
      control=c(rep(1,length(non[1,])), rep(0,length(mut[1,])))
      treat=c(rep(0,length(non[1,])), rep(1,length(mut[1,])))
      limma.result.gm[[gene]] <- rna.diff(exp.reorder,
                                  control, treat,"gm", cutoff)
      ## call rna.diff
      if (draw){
              patient.category <- rep("red", length(exp.match[1,]))
      patient.category[patient.match.mutated[[gene]]] <- "blue"

      boxplot(exp.match, col=patient.category,outline=F)

      t.p.value <- apply(exp.match, 1,
                         function(x) {t.test(x[patient.match.mutated[[gene]]], x[-patient.match.mutated[[gene]]])$p.value})
      w.p.value <- apply(exp.match, 1,
                         function(x) {wilcox.test(x[patient.match.mutated[[gene]]], x[-patient.match.mutated[[gene]]])$p.value})
      length(w.p.value[w.p.value<cutoff])
      length(t.p.value[t.p.value<cutoff])

      ## write out to results folder
      write.table(sort(t.p.value), file=paste("../results/ttest_colon", gene, type, "mutated_exp.txt", collapse="", sep=""), quote=F, sep="\t", col.names = F)
      write.table(sort(w.p.value), file=paste("../results/wilcoxon_colon", gene, type, "mutated_exp.txt", collapse="", sep=""), quote=F, sep="\t", col.names = F)
      ## Gaussian distribution QQ plot
      # qqnorm(t.p.value, col='red', main=paste(gene,"unnormalized by exon")); qqline(t.p.value, col='blue')
      ## Uniform distribution QQ plot
      t.p.value.unif <- runif(length(t.p.value), min=min(t.p.value), max=max(t.p.value))
      w.p.value.unif <- runif(length(w.p.value), min=min(w.p.value), max=max(w.p.value))

      ## qqplot(t.p.value.unif, t.p.value, col=2, main=paste(gene,"qqplot t.test p value versus uniform"))
      ## lines(loess(sort(t.p.value)~sort(t.p.value.unif)), col="Blue", lty=1, lwd=1)
      ## abline(lm(sort(t.p.value)~sort(t.p.value.unif)), col="Blue", lty=1, lwd=1)
      ## qqplot(w.p.value.unif, w.p.value, col=2, main=paste(gene, "qqplot wilcoxon test p value versus uniform"))
      ## lines(loess(sort(w.p.value)~sort(w.p.value.unif)), col="Blue", lty=1, lwd=1)
      ## abline(lm(sort(w.p.value)~sort(w.p.value.unif)), col="Blue", lty=1, lwd=1)
      ## scatterplot for uniform t.p.value and observed t.p.value
      ## library(rgl)
      ## scatter3d(sort(seq(length(t.p.value))), sort(t.p.value.unif), sort(t.p.value))
      ## smoothScatter(sort(t.p.value))
      scatterplot(log2(sort(t.p.value.unif)),log2(sort(t.p.value)), main=paste("Colon cancer t.test scatterplot",gene))
      scatterplot(log2(sort(w.p.value.unif)),log2(sort(w.p.value)), main=paste("Colon cancer wilcoxon.test scatterplot",gene))
      ## common scatterplot
      ## plot(log2(sort(t.p.value.unif)), log2(sort(t.p.value)),
      ##      main=paste(gene, "common plot of uniform p.value against observed p.value"), pch=2, col="blue", lty=3, lwd=1)
      ## ## linear regression
      ## abline(lm(sort(t.p.value)~sort(t.p.value.unif)), col="red", lty=1, lwd=1)
      ## lines(loess(sort(t.p.value)~sort(t.p.value.unif)), col="yellow", lty=4, lwd=1)
      ## plot(sort(w.p.value.unif), sort(w.p.value),
      ##      main=paste(gene, "common plot of uniform p.value against observed p.value"), pch=2, col="blue", lty=3)
      ## ## linear regression
      ## abline(lm(sort(w.p.value)~sort(w.p.value.unif)), col="red", lty=1, lwd=4)
      ## lines(loess(sort(w.p.value)~sort(w.p.value.unif)), col="yellow", lty=4, lwd=3)
      ## heatmap part
      ## wilcoxon
      w.index <- which(w.p.value<=cutoff)
      #pdf("../results/expression_diff_APC.pdf"); dev.off()
      head(exp.match[which(w.p.value<cutoff), ])
      ## bar colors are patients
      # hclust
      expw <- as.matrix(exp.match[as.numeric(w.index),])
      hc.rows <- hclust(dist(expw))
      hc.cols <- hclust(dist(t(expw)))
      heatmap.2(as.matrix(exp.match[as.numeric(w.index),]), col=brewer.pal(9,"Blues"), main= paste(gene, "wilcoxon cutoff p.value", cutoff),
                trace='none', notecex=0.2, ColSideColors = patient.category, dendrogram = "both", cexRow=0.4, scale="row")

      # t.test
      t.index <- which(t.p.value < cutoff)
      head(exp.match[which(t.p.value<cutoff), ])
      ## Col colors are patients
      heatmap.2(as.matrix(exp.match[as.numeric(t.index),]), col=brewer.pal(9,"Blues"),main= paste(gene, "t.test cutoff p.value", cutoff),
                trace='none', notecex=0.2, ColSideColors = patient.category, dendrogram = "both", cexRow=0.4, scale="row")
      }
  }
  limma.result.gm
}

meancentering <- function(expdata){
  ## x - rowmeans
  ## rowMeans for all mean value of data.frame
  data <- apply(expdata, 1, function(x) (x-mean(x)))
  data
}

## explore mutation frequency vs exon lengths
exon.explore <- function(exondata, Mutation.freq, cutoff){
  ## explore the relationship between exon lengths and mutation frequency
  ## exondata extracted from UCSC hg19_allfields.txt
  ## Mutation.freq from TCGA somatic mutation frequency
  exonlength <- read.table(exondata, sep="\t", header=F)
  colnames(exonlength) <- c("genes", "length")
  mutation.exonlength <- merge(Mutation.freq, exonlength, sort=F, all=F)
  par(mfrow = c(2,2))
  exon.norm <- mutation.exonlength$freq/as.numeric(mutation.exonlength$length)
  mutation.exonlength <- cbind(mutation.exonlength, log2(exon.norm))

  ## mutation.new.rank <- rank(mutation.exonlength$normbyexonlength)
  ## mutation.new.rank
  ## pick the outlier, cutoff 1e+05 exon length for colon cancer
  ## 5 ??? for  validated mutation frequency cutoff


  print(head(mutation.exonlength))
  outlier <- mutation.exonlength[mutation.exonlength$length >= cutoff,]
  plot(freq~length,data=mutation.exonlength, col="blue", main="original mutation frequency vs exon lengths")
  for (i in seq(along = as.vector(outlier$freq))) {
    print(as.vector(outlier$freq[i]))
    print(as.vector(outlier$length[i]))
    print(as.vector(outlier$genes[i]))
    text(as.vector(outlier$length)[i]+10, as.vector(outlier$freq)[i], paste(outlier$genes[i]), cex=.8)
    }

  plot(freq~log2(length),data=mutation.exonlength, col="blue", main="original mutation frequency vs log(exon lengths)")
  ## wrong for mut freq / exonlength and exon length dependency
##   plot(log2(exon.norm)~log2(length), data=mutation.exonlength, col="blue",
##        main="log2(normalized mutation frequency) vs log2(exon lengths)")
##   abline(lm(log2(exon.norm)~log2(length),data=mutation.exonlength))
##   lines(loess(log2(exon.norm)~log2(length),data=mutation.exonlength), col="Blue", lty=1, lwd=3)

##   ## quantile regression not fit for this analysis
##   taus <- c( .50, .75, .95)
##   rqs <- as.list(taus)
##   plot(log2(exon.norm)~log2(length), data=mutation.exonlength, col="blue",
##        main="quantile regression log2(normalized mutation frequency) vs log2(exon lengths)",
##        ylim= c(-20, 20))
##   library(quantreg)
##   for (i in seq(along=taus)){
##     rqs[[i]] = rq(log2(mutation.exonlength$length)~mutation.exonlength[,4], tau=taus[i])
##     lines(log2(mutation.exonlength$length), fitted(rqs[[i]]), col=i+1, lwd=0.3, lty=1)}
##   legend("topright", paste("tau=", taus), inset=.04, lty=1, col=2:(length(taus)+1), border="white")
    
##   print(rqs)
##   outliertest <- outlierTest(lm(mutation.exonlength[,4]~log2(mutation.exonlength$length)), cutoff = 0.01)
## ##     rstudent unadjusted p-value Bonferonni p
## ## 1 6.375859         1.9770e-10   1.0265e-06
## ## 3 5.167763         2.4579e-07   1.2761e-03
## ## 5 5.000464         5.9076e-07   3.0672e-03
## ## 2 4.884012         1.0707e-06   5.5591e-03
## ## 6 4.797797         1.6491e-06   8.5624e-03
##   output <- list(exon=mutation.exonlength, outlier = outlier, test=outliertest,rq=rqs)
##   output
}

