setwd('miRNA')

file_man <- read.delim("file_manifest.txt")
file_man <- file_man[-c(1,2),]

head(file_man$Sample)

setwd("./miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/")

dirs <- dir(".")
dir_get <- grep("mirna_quantification", dirs, value=T)

patient <- c()

for (f in seq(along.with=dir_get)) {
  patient[f] <- as.vector(subset(file_man, grepl(dir_get[f], as.vector(file_man[,7])), select=5)[1,1])
}


files_result<- sapply(matrix(dir_get), read.table, header=T, sep="\t")

dir_get == colnames(files_result)
patient == substr(dir_get, 35 , 49)

colnames(files_result) <- patient

result.per<- apply(files_result, 2, function(x) x$reads_per_million_miRNA_mapped)
rownames(result.per) <- files_result[,1]$miRNA_ID

patient_all <- unique(colnames(result.per))


## (result.per[,names(which(table(patient) >= 2))])

dup <- names(which(table(patient) >= 2))

dup.index <- apply(matrix(dup), 1, grep, patient)

colnames(dup.index) <- dup

dup.mean <- list()

for (i in seq(along.with=dup.index[1,])) {
  dup.mean[[colnames(dup.index)[i]]] <- as.vector(apply(result.per[,dup.index[,i]], 1, mean))
}

dup.df <- data.frame(dup.mean)

result <- result.per[ , -as.vector(dup.index)]
result <- cbind(result, dup.df)
     
setdiff(patient, patient_all) ## the same

write.table(result, file="~/TCGA/miRNA_Hiseq.txt", quote=F, sep="\t")
