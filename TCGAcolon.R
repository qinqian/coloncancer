## Author: Qin Qian
## Time-stamp: < modified by qinq :2012-12-06 01:22:19 >
## TCGA process
## Usage: 

## library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
print(args)

IO <- function(workdirectory, datatype1, datatype2, ...){
  # datatype2, e.g. mutations
  # datatype1, e.g. expression
  setwd(workdirectory)
  data1 <- read.table(datatype1, header=T, sep="\t", allowEscapes = FALSE)
  data2 <- read.table(datatype2, header=T, sep="\t", allowEscapes = FALSE)
  data <- list(data1=data1, data2=data2)
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


data.orig <- IO(args[1], args[2], args[3])
data.orig <- IO("../data/",
                "colon_cancer_TCGA_agilent_expression.xls",
                "colon_cancer_mutation_all.maf")

# substitute header . to -
colnames(data.orig$data1) <- substr(gsub('.','-',colnames(data.orig$data1),fixed=T), 1, 15)

# patients mutation
mutation.patient <- substr(unique(data.orig$data2$Tumor_Sample_Barcode), 1, 15)
exp.patient <- colnames(data.orig$data1)
exp.patient <- colnames(data.orig$data1)
mutationexp.match <- match(mutation.patient,exp.patient)
#print(exp.patient[mutationexp.match] == mutation.patient)

# write table of matched expression data
exp.match <- data.orig$data1[, mutationexp.match]
rownames(exp.match) <- data.orig$data1[,1]
write.table(exp.match, row.names=T, col.names=T, file="colon_exp_genes.txt", sep="\t", quote=F)

MutPartition <- function (data, top=10, ...) {
  ## partition Mutation data according to patients and samples
  sort(table(data[,1]), decreasing = T)
}

# Overall statistics of mutations genes in different patient's samples
Mutgenes.freq <- MutPartition(data.orig$data2)
Colon.Mutclass <- data.frame(genes=names(Mutgenes.freq),freq=as.numeric(Mutgenes.freq))
write.table(Colon.Mutclass, file="ColonMutationTop.txt", quote=F, sep="\t", col.names = T)

# separate patient with mutation genes
Mut.Table <- function(mutdata, patients, mutfreq,...){
  patient.Mutfreq.merge <- mutfreq
  for(patientid in seq(along.with=patients)){
    mutation.patient2genes <- mutdata[grep(paste("^", patients[patientid],"*",sep=""), mutdata$Tumor_Sample_Barcode, perl=T), c(1,16)]
    patient.mutationfreq <- table(unique(mutation.patient2genes[,1])) # genes
    mutation.patient.freq <- data.frame(genes=names(patient.mutationfreq), patient=as.vector(patient.mutationfreq))
    mutation.patient.names <- apply(matrix(colnames(mutation.patient.freq)), 1, sub, pattern="patient", replacement=patients[patientid])
    colnames(mutation.patient.freq) <- mutation.patient.names
    patient.Mutfreq.merge <- merge(patient.Mutfreq.merge, mutation.patient.freq)
  }
  patient.Mutfreq.merge
}
mutation.table <- Mut.Table(data.orig$data2, mutation.patient, Colon.Mutclass)
colnames(mutation.table) <- substr(gsub('.','-',colnames(mutation.table),fixed=T), 1, 15)
write.table(mutation.table, file="colon_somatic_mutations.txt", quote=F, sep="\t", col.names = T, row.names=T)

# match mutation genes symbol with expression genes symbol
length(match(unique(data.orig$data2$Hugo_Symbol),Colon.Mutclass$genes))


## use top n mutations genes for patient classification
topN.class <- function(N, mutation.freq, mutation.table){
topN <- head(mutation.freq, N)$genes
top.genes.patient <- merge(top100, mutation.table, all=F, sort=F)
top.genes.patient
#head(top.genes.patient[top.genes.patient$genes=='APC', ])
}

toptable <- topN.class(100, Colon.Mutclass, mutation.table)
patient.match.mutated <- apply(toptable, 1, function(x) {which(as.numeric(x[c(-1,-2)]) >= 1)})
names(patient.match.mutated) <- toptable$genes

## checking data match
length(toptable[1,][as.numeric(toptable[1,]) >=1])
colnames(top.genes.patient)[c(-1,-2)] == colnames(exp.match)

exp.patientclass <- function(){
}

# for APC 
APC.mut <- exp.match[, patient.match.mutated$APC]
APC.non <- exp.match[, -patient.match.mutated$APC]
t.p.value <- apply(exp.match, 1,
                   function(x) {t.test(x[patient.match.mutated$APC], x[-patient.match.mutated$APC], conf=.9)$p.value})
qqnorm(t.p.value); qqline(t.p.value)
