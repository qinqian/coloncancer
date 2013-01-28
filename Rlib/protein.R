
file_man <- read.delim("file_manifest.txt")

file_man <- file_man[-c(1,2),]
head(file_man$Sample)


prefix <- "Expression-Protein/MDA__MDA_RPPA_Core/Level_2"
setwd(prefix)
dirs <- dir(".")

files_result<- sapply(matrix(dirs), read.table, header=T, sep="\t", row.names=1)

(file_man$File.Name) == colnames(files_result)
colnames(files_result) <- (file_man$Sample)

write.table(files_result, file="~/TCGA/colon_protein_level2.txt", quote=F, sep="\t")
