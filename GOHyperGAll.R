######################################################################################################
### GOHyperGAll: Global Hypergeometric Test Using Custom Gene-to-GO Mappings Plus GO Slim Analysis ###
######################################################################################################
## Author: Thomas Girke
## Last update: Feb 7, 2008
## Utility: To test a sample population of genes for over-representation of GO terms, the 
## function 'GOHyperGAll' computes for all GO nodes a hypergeometric distribution test and 
## returns the corresponding raw and Bonferroni corrected p-values. A subsequent filter function 
## performs a GO Slim analysis using default or custom GO Slim categories. 
## The associated publication is available in Plant Physiol (2008) 147, 41-57.
## Note: GOHyperGAll provides similar utilities as the GOHyperG function in the GOstats package 
## from BioConductor. The main difference is that GOHyperGAll simplifies the usage of custom 
## chip-to-gene and gene-to-GO mappings.
##
## How it works:
## (A) Generate the required data objects (slow, but needs to be done only once)
## (B) Define GOhyperG_All function
## (C) Subsetting and plotting of results by assigned nodes or goSlim categories
## To demo the script and import all required functions, run the following source() command:
##         source("http://bioinfo.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/GOHyperGAll.txt")
## HTML Instructions: 
##       http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/R_BioCondManual.html#GOHyperGAll

#########################################
## (A) Generate the required data objects        
#########################################
## (A.1) Generate sample data frames with assigned gene-to-GO mappings,
## one for MF, one for BP and one for CC mappings
## custom mappings can be used here, but need to have the same format as GO_XX_DFs in the following examples
## (A.1.1) Obtain mappings from geneontology.org
readGOorg <- function(myfile = "gene_association.tair", colno = c(5,11,9), org) {
        go_org <- read.delim(myfile, na.strings = "", header=F, comment.char = "!", sep="\t")
        go_org <- go_org[ , colno]
        names(go_org) <- c("GOID", "GeneID", "GOCAT")
	if(org == "Arabidopsis") {
		go_org[,"GeneID"] <- gsub(".*(AT.G\\d\\d\\d\\d\\d).*", "\\1", as.character(go_org[,2]), perl=T)
        	go_org <- go_org[grep("^AT.G\\d\\d\\d\\d\\d", as.character(go_org$GeneID), perl=T),]
        	go_org <- go_org[!duplicated(paste(go_org[,"GOID"], gsub("\\.\\d{1,}", "", as.character(go_org[,"GeneID"]), perl=T), sep="_")),]
        }
        go_org <- na.omit(go_org)
        GO_MF_DF <<- go_org[go_org[,3]=="F",]
        write.table(GO_MF_DF, file="GO_MF_DF", quote=T, sep="\t")
        cat("\n", "Object 'GO_MF_DF' created containing assigned gene-to-MFGO mappings. To use custom mappings, generate data frame with the same structure in col 1-2.", "\n")
        GO_BP_DF <<- go_org[go_org[,3]=="P",]
        write.table(GO_BP_DF, file="GO_BP_DF", quote=T, sep="\t")
        cat("\n", "Object 'GO_BP_DF' created containing assigned gene-to-MFGO mappings. To use custom mappings, generate data frame with the same structure in col 1-2.", "\n")
        GO_CC_DF <<- go_org[go_org[,3]=="C",]
        write.table(GO_CC_DF, file="GO_CC_DF", quote=T, sep="\t")
        cat("\n", "Object 'GO_CC_DF' created containing assigned gene-to-MFGO mappings. To use custom mappings, generate data frame with the same structure in col 1-2.", "\n")
	
	## Generates "go_df" data frame containing the commonly used components for all GO nodes: GOID, GO Term and Ontology Type.
        require(GOstats); require(GO.db)
        go_df <- data.frame(GOID=unlist(eapply(GOTERM, function(x) x@GOID)), Term=unlist(eapply(GOTERM, function(x) x@Term)), Ont=unlist(eapply(GOTERM, function(x) x@Ontology))) 
        go_df <- na.omit(go_df)
	go_df <<- go_df
	write.table(go_df, file="go_df", quote=T, sep="\t")
        cat("\n", "Object 'go_df' created containing for all GO nodes the commonly used components: GOID, GO Term and Ontology Type", "\n")
}

## (A.1.1b) Convert GO objects from GOHyperGAll into GSEA format
## This step is not required for for GOHyperGAll approach. It simply converts XX_node_affy_list or GO_XX_DF files into *.gmt 
## formatted files that can be imported into GSEA from the Broad Institute.
GOhyper2GSEA <- function(myfile=c("MF_node_affy_list", "BP_node_affy_list", "CC_node_affy_list"), type="all") {
	if(type=="all"){  
		for(i in 1:length(myfile)) { 
			mynames <- gsub("(..)_.*", "GO_\\1_ALL", myfile)
			load(file=myfile[i])
			GO_List <- eval(parse(text=myfile[i]))
			myindex <- sapply(GO_List, length)
			GO_List <- GO_List[myindex!=0]
			GO_List <- lapply(GO_List, paste, collapse="\t")
			exportDF <- data.frame(GS_ID=names(GO_List), Desc=rep("NA", length(GO_List)), GeneID=as.vector(unlist(GO_List)))
			write.table(exportDF, file=paste(mynames[i], ".gmt", sep=""), row.names=F, col.names=F, quote=F, sep="\t")
			cat("Saved file:", mynames[i], "\n")
		}
	}
	if(type=="terminal") {
		for(i in 1:length(myfile)) { 
			mynames <- gsub("_DF", "_TERM", myfile)
			GO_DF <- read.delim(file=myfile[i])
			GO_List <- tapply(as.vector(GO_DF[,2]), as.factor(as.vector(GO_DF[ ,1])), as.vector)
			GO_List <- lapply(GO_List, paste, collapse="\t")
			exportDF <- data.frame(GS_ID=names(GO_List), Desc=rep("NA", length(GO_List)), GeneID=as.vector(unlist(GO_List)))
			write.table(exportDF, file=paste(mynames[i], ".gmt", sep=""), row.names=F, col.names=F, quote=F, sep="\t")
			cat("Saved file:", mynames[i], "\n")
		}
	}
}

## (A.1.2) Obtain mappings from BioC
sampleDFgene2GO <- function(lib="ath1121501.db") {
        require(GOstats); require(GO.db); require(Annotate); require(lib, character.only=T)
	mylibbase <- gsub(".db", "", lib) 
        affyGOMF <- eapply(get(paste(mylibbase, "GO", sep="")), getOntology, "MF") # generates list with GeneID components containing MFGOs
        GO_MF_DF <<- data.frame(GOID=unlist(affyGOMF), GeneID=rep(names(affyGOMF), as.vector(sapply(affyGOMF, length))), Count=rep(as.vector(sapply(affyGOMF, length)), as.vector(sapply(affyGOMF, length))))
        write.table(GO_MF_DF, file="GO_MF_DF", quote=F, sep="\t")
        cat("\n", "Object 'GO_MF_DF' created containing assigned gene-to-MFGO mappings. To use custom mappings, generate data frame with the same structure in col 1-2.", "\n")
        affyGOBP <- eapply(get(paste(mylibbase, "GO", sep="")), getOntology, "BP") # generates list with GeneID components containing BPGOs
        GO_BP_DF <<- data.frame(GOID=unlist(affyGOBP), GeneID=rep(names(affyGOBP), as.vector(sapply(affyGOBP, length))), Count=rep(as.vector(sapply(affyGOBP, length)), as.vector(sapply(affyGOBP, length))))
        write.table(GO_BP_DF, file="GO_BP_DF", quote=F, sep="\t")
        cat("\n", "Object 'GO_BP_DF' created containing assigned gene-to-BPGO mappings. To use custom mappings, generate data frame with the same structure in col 1-2.", "\n")
        affyGOCC <- eapply(get(paste(mylibbase, "GO", sep="")), getOntology, "CC") # generates list with GeneID components containing CCGOs
        GO_CC_DF <<- data.frame(GOID=unlist(affyGOCC), GeneID=rep(names(affyGOCC), as.vector(sapply(affyGOCC, length))), Count=rep(as.vector(sapply(affyGOCC, length)), as.vector(sapply(affyGOCC, length))))
        write.table(GO_CC_DF, file="GO_CC_DF", quote=F, sep="\t")
        cat("\n", "Object 'GO_CC_DF' created containing assigned gene-to-CCGO mappings. To use custom mappings, generate data frame with the same structure in col 1-2.", "\n")
	
	## Generates "go_df" data frame containing the commonly used components for all GO nodes: GOID, GO Term and Ontology Type.
        require(GOstats); require(GO.db)
        go_df <- data.frame(GOID=unlist(eapply(GOTERM, function(x) x@GOID)), Term=unlist(eapply(GOTERM, function(x) x@Term)), Ont=unlist(eapply(GOTERM, function(x) x@Ontology))) 
        go_df <- na.omit(go_df)
	go_df <<- go_df
	write.table(go_df, file="go_df", quote=T, sep="\t")
        cat("\n", "Object 'go_df' created containing for all GO nodes the commonly used components: GOID, GO Term and Ontology Type", "\n")
}
cat("\n", "(A.1) To use the GOHyperGAll() function, one needs 4 data frames containing the gene-to-GO mappings MF, BP, CC and the GO terms. \n       Demo data sets from GO.org or BioC can be created with one of these commands: \n \t (A.1.1) For annotations from geneontology.org: \n \t readGOorg(myfile = \"gene_association.tair\", colno = c(5,11,9), org = \"Arabidopsis\") \n \t \t myfile: download annotation table from geneontology.org and unzip it. Then point function to file name. \n \t \t colno: required column numbers; default 'c(5,11,9)' should work in most cases \n \t \t org: \"Arabidopsis\" or any string for other organisms \n \t (A.1.2) For annotations from BioC: \n \t sampleDFgene2GO(lib=\"ath1121501\") \n \t \t lib: defines annotation library, default is \"ath1121501\"", "\n")

## (A.2) Generate list containing gene-to-GO-OFFSPRING associations including assiged nodes
## This is very slow (3x3 minutes), but needs to be done only once! 
gene2GOlist <- function(rootUK=T) { # If the argument 'rootUK' is set to TRUE then the root nodes are treated as terminal nodes to account for the new unknown terms 
        require(GOstats); require(GO.db)
        for(i in c("MF","BP","CC")) {
                if(i=="MF") {
                        go_offspr_list <- as.list(GOMFOFFSPRING) }
                if(i=="BP") {
                        go_offspr_list <- as.list(GOBPOFFSPRING) }
                if(i=="CC") {
                        go_offspr_list <- as.list(GOCCOFFSPRING) }
                go_offspr_list <- lapply(go_offspr_list, unlist); go_offspr_list <- lapply(go_offspr_list, as.vector) # clean-up step for the list
                go_offspr_list_temp <- lapply(names(go_offspr_list), function(x) c(x, go_offspr_list[[x]]) ) # include list component (GOID) names in corresponding (GOID) vectors
                names(go_offspr_list_temp) <- names(go_offspr_list) # names list components after go_offspr_list
                go_offspr_list <- go_offspr_list_temp
                go_offspr_list <- lapply(go_offspr_list, function(x) x[!is.na(x)]) # remove NAs in vectors
                
		## Treat root nodes as terminal nodes to account for the new unknown terms. This step removes the offspring information from the root nodes.
		if(rootUK==T) {
			if(i=="MF") { go_offspr_list[["GO:0003674"]] <- c("GO:0003674") } 
			if(i=="BP") { go_offspr_list[["GO:0008150"]] <- c("GO:0008150") }
			if(i=="CC") { go_offspr_list[["GO:0005575"]] <- c("GO:0005575") }
		}

		## Retrieve gene/affy IDs for GOID vectors
		if(i=="MF") {
                        MF_node_affy_list <<- lapply(go_offspr_list, function(x) unique(as.vector(GO_MF_DF[GO_MF_DF$GOID %in% x, 2]))) 
                        save(MF_node_affy_list, file="MF_node_affy_list") }
                if(i=="BP") {
                        BP_node_affy_list <<- lapply(go_offspr_list, function(x) unique(as.vector(GO_BP_DF[GO_BP_DF$GOID %in% x, 2]))) 
                        save(BP_node_affy_list, file="BP_node_affy_list") }
                if(i=="CC") {
                        CC_node_affy_list <<- lapply(go_offspr_list, function(x) unique(as.vector(GO_CC_DF[GO_CC_DF$GOID %in% x, 2])))  
                        save(CC_node_affy_list, file="CC_node_affy_list") }
                cat("\n", paste("Object '", i, "_node_affy_list'", sep=""), "with gene-to-GO-OFFSPRING associations created and saved in your working directory.", "\n") 
        }
}
cat("\n", "(A.2) The corresponding gene-to-GO-OFFSPRING associations are created from the three data frames with the following \n       command. This creates 3 list objects with the required MF, BP and CC associations. \n \t gene2GOlist(rootUK=T)", "\n")

## (A.3) Generate AffyID-to-GeneID mappings when working with chip feature IDs
## This function creates a AffyID-to-GeneID mapping data frame using by default the TAIR mappings for the Arabidopsis ATH1 chip. 
## Once the decoding data frame 'affy2locusDF' is created, the function returns for a query set of AffyIDs the corresponding GeneIDs.
## To use the function for the mappings of other chips, one needs to create the corresponding decoding data frame 'affy2locusDF'.
AffyID2GeneID <- function(map = "ftp://ftp.arabidopsis.org/home/tair/Microarrays/Affymetrix/affy_ATH1_array_elements-2008-5-29.txt", affyIDs, probe2gene=1) {
        if(!exists("affy2locusDF")) {
                cat("\n", "Downloading AffyID-to-GeneID mappings, creating object 'affy2locusDF' and saving it in your working directory", "\n")
                affy2locus <- read.delim(map, na.strings = "", fill=TRUE, header=T, sep="\t")[,-c(2:4,7:9)]
                names(affy2locus) <- c("AffyID", "AGI", "Desc")
                row.names(affy2locus) <- as.vector(affy2locus[,1])
                my_list <- apply(affy2locus[,-c(3)], 1, list); my_list <- lapply(my_list, unlist)
                my_list <- lapply(my_list, function(x) as.vector(x[-1]))
                my_list <- lapply(my_list, strsplit, ";"); my_list <- lapply(my_list, unlist)
                affy2locusDF <- data.frame(unlist(my_list))
                affy2locusDF <- data.frame(rep(names(unlist(lapply(my_list, length))), as.vector(unlist(lapply(my_list, length)))), affy2locusDF)
                names(affy2locusDF) <- c("AffyID", "GeneID")
                affy2locusDF <<- affy2locusDF
        	write.table(affy2locusDF, file="affy2locusDF", quote=F, sep="\t")
        }
        if(!missing(affyIDs)) {
		if(probe2gene==1) { # For probe sets that match several loci, only the first locus ID will be used
			affy2locusDF <- affy2locusDF[!duplicated(affy2locusDF$AffyID),]
		}	
        	GeneIDs <- unique(as.vector(affy2locusDF[affy2locusDF[,1] %in% affyIDs, 2]))
        	return(GeneIDs)
        }
}
cat("\n", "(A.3) To work with AffyIDs, the function AffyID2GeneID() can be used to import custom AffyID-to-GeneID mappings. \n \t AffyID2GeneID(map = \"ftp://ftp.arabidopsis.org/home/tair/Microarrays/Affymetrix/affy_ATH1_array_elements-2006-07-14.txt\") \n \t \t map: location of custom AffyID-to-GeneID mappings", "\n")

## (A.4) Next time things are much faster by reading the 6 data objects from file
loadData <- function() {
                need_affy2gene <- dir()
        	if(length(need_affy2gene[need_affy2gene=="affy2locusDF"]) > 0) {
        		affy2locusDF <<- read.table(file="affy2locusDF", header=T, colClasses = "character")
        	}
                GO_MF_DF <<- read.table(file="GO_MF_DF", header=T, colClasses = "character")
                GO_BP_DF <<- read.table(file="GO_BP_DF", header=T, colClasses = "character")
                GO_CC_DF <<- read.table(file="GO_CC_DF", header=T, colClasses = "character")
                if(any(dir() %in% "go_df")) { go_df <<- read.table(file="go_df", header=T, colClasses = "character") }
}
cat("\n", "(A.4) In future R sessions one can can omit the previous 3 steps (A.1-A.3) by importing all 6 (7) data objects like this: \n \t loadData(); load(file=\"MF_node_affy_list\"); load(file=\"BP_node_affy_list\"); load(file=\"CC_node_affy_list\")", "\n")

############################
## (B) GOhyperG_All function
############################
## (B.1) Define GOhyperG_All function
GOHyperGAll <- function(gocat="MF", sample, Nannot=2) {
## Generates data frame (go_df) containing the commonly used components for all GO nodes: GOID, GO Term and Ontology Type. This step is only required if "go_df" hasn't been imported with the above load() function.
        require(GOstats); require(GO.db)
        if(!exists("go_df")) {
        	go_df <- data.frame(GOID=unlist(eapply(GOTERM, function(x) x@GOID)), Term=unlist(eapply(GOTERM, function(x) x@Term)), Ont=unlist(eapply(GOTERM, function(x) x@Ontology))) 
        	go_df <<- na.omit(go_df) 
        	cat("\n", "Object 'go_df' created containing for all GO nodes the commonly used components: GOID, GO Term and Ontology Type", "\n")
        }
## (m): Obtain for every node in GO tree their number of associated genes or chip features
        if(gocat=="MF") {node_affy_list <- MF_node_affy_list}
        if(gocat=="BP") {node_affy_list <- BP_node_affy_list}
        if(gocat=="CC") {node_affy_list <- CC_node_affy_list}
        node_stats_df <- data.frame(NodeSize=sapply(node_affy_list, length))
        node_stats_df <- data.frame(GOID=row.names(node_stats_df), node_stats_df)
        row.names(node_stats_df) <- 1:length(node_stats_df[,1])
        m <- as.vector(node_stats_df$NodeSize)        
                
## (x): Obtain for every node in GO tree the number of matching genes in sample set
        node_sample_stats <- sapply(node_affy_list, function(x) { sum(unlist(x) %in% sample) } )
        node_sample_stats <- as.vector(node_sample_stats)
        x <- node_sample_stats        

## (n): Obtain the number of unique genes at GO nodes with direct annotations
        if(gocat=="MF") { GO_DF <- GO_MF_DF }
        if(gocat=="BP") { GO_DF <- GO_BP_DF }
        if(gocat=="CC") { GO_DF <- GO_CC_DF }
        n <- length(unique(GO_DF[, 2]))

## (k): Obtain number of unique genes in test sample that have GO mappings
        k <- length(unique(GO_DF[GO_DF[,2] %in% sample, 2]))

## Obtain gene/chip keys matching at GO nodes
        match_key <- sapply(node_affy_list, function(x) { x[unlist(x) %in% sample] } )
        match_key <- sapply(match_key, function(x) { paste(x, collapse=" ") } )
	match_key <- as.vector(match_key)
	key <- match_key; key[key==""] <- "NA"

## Apply phyper function
        phyp_v <- phyper(x-1, m, n-m , k, lower.tail = FALSE)

## P-value correction according to Bioinformatics, 20, 3710-3715
	Ncorrect <- table(GO_DF[GO_DF$GeneID %in% sample, 1]) # Obtain the GO nodes with direct annotations from sample set 
	Ncorrect <- sum(Ncorrect >= Nannot) # Count only those that have 2 or more annotations from sample set
	if(Ncorrect<=1) { 
		adj_phyp_v <- phyp_v # no adjustment necessary if Ncorrect <= 1
	} else {
		adj_phyp_v <- phyp_v * Ncorrect # Calculates simple Bonferroni correction. 
		adj_phyp_v[adj_phyp_v >= 1] <- 1
		# adj_phyp_v <- sapply(phyp_v, p.adjust, method=Padj, n = Ncorrect) # Runs p.adjust(). This is disabled because most adjustment methods require that the length of the p-value vector is >= n.  
	}
	
## Generate output data format
	result_df <- data.frame(node_stats_df, SampleMatch=x, Phyper=phyp_v, Padj=adj_phyp_v, SampleKeys=key)
        result_df <- merge(result_df, go_df, x.by="GOID", y.by="GOID", all.x=T)
        result_df <- result_df[order(result_df$Phyper), ]
	result_df <- result_df[,c(1:5,7:8,6)]
        result_df
}
cat("\n", "(B.1) The function GOHyperGAll() runs the phyper test against all nodes in the GO network. \n Usage: \n \t GOHyperGAll(gocat=\"MF\", sample=test_sample, Nannot=2)[1:20,] \n \t \t gocat: \"MF\", \"BP\" or \"CC\" \n \t \t Nannot: minimum number of direct annotations for p-value adjustment \n \t \t test_sample <- unique(as.vector(GO_MF_DF[1:40,2])) # for GeneIDs\n \t \t test_sample <- AffyID2GeneID(affyIDs=affy_sample, probe2gene=1) # for AffyIDs \n \t \t affy_sample <- c(\"266592_at\", \"266703_at\", \"266199_at\", \"246949_at\", \"267370_at\", \"267115_s_at\", \"266489_at\", \"259845_at\", \"266295_at\", \"262632_at\")", "\n")

####################################################################################
## (C) Subsetting of results from GOHyperGAll by assigned nodes or goSlim categories
####################################################################################
## (C.1) Define subsetting function
        GOHyperGAll_Subset <- function(GOHyperGAll_result, sample=test_sample, type="goSlim", myslimv) { # type: "goSlim" or "assigned"; optional argument "myslimv" to privde custom goSlim vector
                if(type=="goSlim") {
                        if(missing(myslimv)) {
                                slimv <- c("GO:0003674", "GO:0008150", "GO:0005575", "GO:0030246","GO:0008289","GO:0003676","GO:0000166","GO:0019825","GO:0005515","GO:0003824","GO:0030234","GO:0003774","GO:0004871","GO:0005198","GO:0030528","GO:0045182","GO:0005215","GO:0006519","GO:0007154","GO:0016043","GO:0006412","GO:0006464","GO:0006810","GO:0007275","GO:0007049","GO:0005975","GO:0006629","GO:0006139","GO:0019748","GO:0015979","GO:0005618","GO:0005829","GO:0005783","GO:0005768","GO:0005794","GO:0005739","GO:0005777","GO:0009536","GO:0005840","GO:0005773","GO:0005764","GO:0005856","GO:0005634","GO:0005886","GO:0005576") # contains new unknown terms: "GO:0003674", "GO:0008150", "GO:0005575"
                                # slimv <- c("GO:0005554", "GO:0000004", "GO:0008372", "GO:0030246","GO:0008289","GO:0003676","GO:0000166","GO:0019825","GO:0005515","GO:0003824","GO:0030234","GO:0003774","GO:0004871","GO:0005198","GO:0030528","GO:0045182","GO:0005215","GO:0006519","GO:0007154","GO:0016043","GO:0006412","GO:0006464","GO:0006810","GO:0007275","GO:0007049","GO:0005975","GO:0006629","GO:0006139","GO:0019748","GO:0015979","GO:0005618","GO:0005829","GO:0005783","GO:0005768","GO:0005794","GO:0005739","GO:0005777","GO:0009536","GO:0005840","GO:0005773","GO:0005764","GO:0005856","GO:0005634","GO:0005886","GO:0005576") # contains old unknown terms: "GO:0005554", "GO:0000004", "GO:0008372" 
                                } else { 
                                slimv <- myslimv } 
                        GOHyperGAll_subset <- GOHyperGAll_result[GOHyperGAll_result[,1] %in% slimv, ]
                }
                if(type=="assigned") {
                        termGO <- c(as.vector(GO_MF_DF[GO_MF_DF$GeneID %in% sample, 1]), 
                                as.vector(GO_BP_DF[GO_BP_DF$GeneID %in% sample, 1]), 
                                as.vector(GO_CC_DF[GO_CC_DF$GeneID %in% sample, 1]))
                        subset_v <- unique(termGO)
                        GOHyperGAll_subset <- GOHyperGAll_result[GOHyperGAll_result[,1] %in% subset_v, ]
                }
                GOHyperGAll_subset
        }
cat("\n", "(C.1) The function GOHyperGAll_Subset() allows subsetting of the GOHyperGAll() results by assigned GO nodes or custom goSlim categories.", "\n", "Usage:", "\n", "\t GOHyperGAll_result <- GOHyperGAll(gocat=\"MF\", sample=test_sample, Nannot=2)", "\n", "\t GOHyperGAll_Subset(GOHyperGAll_result, sample=test_sample, type=\"goSlim\")", "\n", "\t \t type: \"goSlim\" or \"assigned\"", "\n", "\t \t myslimv: optional argument allows usage of a custom goSlim vector", "\n")

## Apply subsetting function
        # GOHyperGAll_Subset(GOHyperGAll_result, sample=test_sample, type="goSlim")

## (C.2) Plotting of subsetted results
        # subset <- GOHyperGAll_Subset(GOHyperGAll_result, sample=test_sample, type="goSlim")
        # pie(subset[subset$SampleMatch>0 ,3], labels=as.vector(subset[subset$SampleMatch>0 ,1]), main=unique(as.vector(subset[subset$SampleMatch>0 ,6])))
cat("\n", "(C.2) Plot pie chart of subsetted results: \n \t subset <- GOHyperGAll_Subset(GOHyperGAll_result, sample=test_sample, type=\"goSlim\") \n \t pie(subset[subset$SampleMatch>0 ,3], labels=as.vector(subset[subset$SampleMatch>0 ,1]), main=unique(as.vector(subset[subset$SampleMatch>0, 7])))", "\n")

#########################################################
## (D) Reduce GO Term Redundancy in 'GOHyperGAll_results' 
#########################################################
## (D.1) The function 'GOHyperGAll_Simplify' subsets the data frame 'GOHyperGAll_result' by a user 
## specified adjusted p-value cutoff and removes from it all GO nodes with overlapping children sets 
## (OFFSPRING). Only the best scoring nodes remain in the data frame. 
## The argument 'correct' is experimental. It aims to favor the selection of distal (information rich) 
## GO terms that have at the same time a large number of sample matches. The following calculation is used 
## for this adjustment: phyper x Number_of_children / SampleMatch
## Define GOHyperGAll_Simplify()
GOHyperGAll_Simplify <- function(GOHyperGAll_result, gocat="MF", cutoff=0.001, correct=T) { # gocat: "MF", "BP" or "CC"; cutoff: p-value cutoff; correct: TRUE or FALSE 
	if(gocat!=as.vector(GOHyperGAll_result$Ont[!is.na(GOHyperGAll_result$Ont)])[1]) { stop("The GO categories in GOHyperGAll_Simplify() and GOHyperGAll_result need to match") }
	testDF <- GOHyperGAll_result[GOHyperGAll_result$Padj<=cutoff,]
	testDF <- data.frame(testDF, test=rep(0, times=length(testDF[,1])))
	testDF <- testDF[!is.na(testDF$Ont),]
	GOIDv <- NULL
	GO_OL_Matchv <- NULL
	while(sum(testDF$test==0)>0) {
		clusterv <- NULL
		test <- as.vector(testDF[,1])
		for(j in 1:length(test)) {
			if(gocat=="MF") { mymatch <- sum(unique(na.omit(c(test[j], as.list(GOMFOFFSPRING)[[test[j]]])) %in% na.omit(c(test[1], as.list(GOMFOFFSPRING)[[test[1]]]))))
				if(mymatch==1) { mymatch <- length(as.list(GOMFOFFSPRING)[[test[j]]])  }
			}
			if(gocat=="BP") { mymatch <- sum(unique(na.omit(c(test[j], as.list(GOBPOFFSPRING)[[test[j]]])) %in% na.omit(c(test[1], as.list(GOBPOFFSPRING)[[test[1]]]))))
				if(mymatch==1) { mymatch <- length(as.list(GOBPOFFSPRING)[[test[j]]])  }
			}
			if(gocat=="CC") { mymatch <- sum(unique(na.omit(c(test[j], as.list(GOCCOFFSPRING)[[test[j]]])) %in% na.omit(c(test[1], as.list(GOCCOFFSPRING)[[test[1]]]))))
				if(mymatch==1) { mymatch <- length(as.list(GOCCOFFSPRING)[[test[j]]])  }
			}
			clusterv <- c(clusterv, mymatch)
		}
		clusterv[clusterv==0] <- NA
		testDF <- data.frame(testDF[,-9], test=clusterv)
		if(correct==T) { 
			testDF <- data.frame(testDF, decide=testDF$Padj * (testDF$test/testDF$SampleMatch)) 
			} else {
			testDF <- data.frame(testDF, decide=testDF$Padj) }
		GOIDv <- c(GOIDv, as.vector(testDF[order(testDF[,10]),][1,1]))
		GO_OL_Matchv <- c(GO_OL_Matchv, length(unique(unlist(strsplit(as.vector(testDF[!is.na(testDF$test),8]), " ")))))
		testDF <- testDF[is.na(testDF$test),]
		testDF <- testDF[order(testDF[,5]),-c(9,10)]
		testDF <- data.frame(testDF, test=rep(0, times=length(testDF[,1])))
		cat(GOIDv, "\n")
	}
	simplifyDF <- data.frame(GOID=GOIDv, GO_OL_Match=GO_OL_Matchv)
	simplifyDF
}

## Apply GOHyperGAll_Simplify
## simplifyDF <- GOHyperGAll_Simplify(GOHyperGAll_result, gocat="MF", cutoff=0.001, correct=T)
## data.frame(GOHyperGAll_result[GOHyperGAll_result[,1] %in% simplifyDF[,1], -8], GO_OL_Match=simplifyDF[,2])

########################################
## (D.2) Batch Analysis of Gene Clusters
########################################
## The function 'GOCluster_Report' performs the three GO analyses in batch mode: 'GOHyperGAll', 
## 'GOHyperGAll_Subset' or 'GOHyperGAll_Simplify'. It processes many groups of genes (e.g. 
## gene expression clusters) and organizes the results in a single data frame.
## The gene sets need to be provided in a data frame of this format:
## 	probeID/geneID	ClusterID	ClusterSize
##	id1		CL1		2
##	id2		CL1		2
##	id3		CL2		1
##	...		...		...   
##
## Define 'GOCluster_Report()'

GOCluster_Report <- function(CL_DF=CL_DF, id_type="affy", method="all", CLSZ=10, cutoff=0.001, gocats=c("MF", "BP", "CC"), myslimv="default", correct=TRUE, recordSpecGO=NULL, ...) { # CLSZ: minimum cluster size; method: "all", "slim" or "simplify"; gocat: "MF", "BP" or "CC"; cutoff: adjusted p-value cutoff; recordSpecGO: argument to include one specific GOID in each of the 3 ontologies, e.g: recordSpecGO=c("GO:0003674", "GO:0008150", "GO:0005575")
        cluster_loop <- unique(as.vector(CL_DF[CL_DF[,3]>=CLSZ,2]))
        if(length(cluster_loop[grep("CL", cluster_loop)])>0) {
                cluster_loop <- paste("CL", sort(as.numeric(gsub("CL","", as.character(cluster_loop)))), sep="") 
        }
	if(method=="all") {
        	containerDF <- data.frame(CLID=NULL, CLSZ=NULL, GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL)
		for(i in cluster_loop) {
			cat("\n", "Processing cluster no", i, "with method: \"all\" (GOHyperGAll) \n")
			if(id_type=="affy") {
				affy_sample <- CL_DF[CL_DF[,2]==i, 1]
				test_sample <- AffyID2GeneID(affyIDs=affy_sample, probe2gene=1)
			}
			if(id_type=="gene") {
				test_sample <- as.vector(CL_DF[CL_DF[,2]==i, 1])
			}	
			containerDF2 <- data.frame(GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL)
			count <- 0
			for(j in gocats) {
				count <- count+1
				GOHyperGAll_result <- GOHyperGAll(gocat=j, sample=test_sample, ...)
				tempDF <- GOHyperGAll_result[GOHyperGAll_result$Padj <= cutoff, ]
				if(length(tempDF[,1])==0) { # If filter returns empty data frame, then include at least the first two best scoring GO entries 
                                        tempDF <- GOHyperGAll_result[1:2,]
                                } 
				containerDF2 <- rbind(containerDF2, tempDF)
				if(length(recordSpecGO)>0) {
					containerDF2 <- rbind(containerDF2, GOHyperGAll_result[GOHyperGAll_result[,1]==recordSpecGO[count],])
				}
				no_annot <- test_sample[!test_sample %in% eval(parse(text=paste("GO_", j, "_DF", sep="")))[,2]]
				no_annot <- no_annot[no_annot!="no_match"]
				containerDF2 <- rbind(containerDF2, data.frame(GOID=paste("no_annot_", j, sep=""), NodeSize=NA, SampleMatch=length(no_annot), Phyper=NA, Padj=NA, Term=NA, Ont=j, SampleKeys=paste(no_annot, collapse=", ")))
			}
			tempDF2 <- data.frame(CLID=rep(i, times=length(containerDF2[,1])), CLSZ=rep(unique(as.vector(CL_DF[CL_DF[,2]==i,3])), times=length(containerDF2[,1])), containerDF2)
			containerDF <- rbind(containerDF, tempDF2)
		}
		return(containerDF)
	}	
	if(method=="slim") {
        	containerDF <- data.frame(CLID=NULL, CLSZ=NULL, GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL)
		for(i in cluster_loop) {
			cat("\n", "Processing cluster no", i, "with method: \"slim\" (GOHyperGAll_Subset) \n")
			if(id_type=="affy") {
				affy_sample <- CL_DF[CL_DF[,2]==i, 1]
				test_sample <- AffyID2GeneID(affyIDs=affy_sample, probe2gene=1)
			}
			if(id_type=="gene") {
				test_sample <- as.vector(CL_DF[CL_DF[,2]==i, 1])
			}	
			containerDF2 <- data.frame(GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL)
			count <- 0
			for(j in gocats) {
				count <- count+1
				GOHyperGAll_result <- GOHyperGAll(gocat=j, sample=test_sample, ...)
                        	if(any(myslimv == "default")) {
                               		slimv <- c("GO:0003674", "GO:0008150", "GO:0005575", "GO:0030246","GO:0008289","GO:0003676","GO:0000166","GO:0019825","GO:0005515","GO:0003824","GO:0030234","GO:0003774","GO:0004871","GO:0005198","GO:0030528","GO:0045182","GO:0005215","GO:0006519","GO:0007154","GO:0016043","GO:0006412","GO:0006464","GO:0006810","GO:0007275","GO:0007049","GO:0005975","GO:0006629","GO:0006139","GO:0019748","GO:0015979","GO:0005618","GO:0005829","GO:0005783","GO:0005768","GO:0005794","GO:0005739","GO:0005777","GO:0009536","GO:0005840","GO:0005773","GO:0005764","GO:0005856","GO:0005634","GO:0005886","GO:0005576") # contains new unknown terms: "GO:0003674", "GO:0008150", "GO:0005575"
                                } else { 
                                	slimv <- myslimv 
				} 
				tempDF <- GOHyperGAll_Subset(GOHyperGAll_result, sample=test_sample, type="goSlim", myslimv=slimv)
				containerDF2 <- rbind(containerDF2, tempDF)
				if(length(recordSpecGO)>0) {
					containerDF2 <- rbind(containerDF2, GOHyperGAll_result[GOHyperGAll_result[,1]==recordSpecGO[count],])
				}
				no_annot <- test_sample[!test_sample %in% eval(parse(text=paste("GO_", j, "_DF", sep="")))[,2]]
				no_annot <- no_annot[no_annot!="no_match"]
				containerDF2 <- rbind(containerDF2, data.frame(GOID=paste("no_annot_", j, sep=""), NodeSize=NA, SampleMatch=length(no_annot), Phyper=NA, Padj=NA, Term=NA, Ont=j, SampleKeys=paste(no_annot, collapse=", ")))
			}
			tempDF2 <- data.frame(CLID=rep(i, times=length(containerDF2[,1])), CLSZ=rep(unique(as.vector(CL_DF[CL_DF[,2]==i,3])), times=length(containerDF2[,1])), containerDF2)
			containerDF <- rbind(containerDF, tempDF2)
		}
		return(containerDF)
	}	
	if(method=="simplify") {
        	containerDF <- data.frame(CLID=NULL, CLSZ=NULL, GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL, GO_OL_Match=NULL)
		for(i in cluster_loop) {
			cat("\n", "Processing cluster no", i, "with method: \"simplify\" (GOHyperGAll_Simplify) \n")
			if(id_type=="affy") {
				affy_sample <- CL_DF[CL_DF[,2]==i, 1]
				test_sample <- AffyID2GeneID(affyIDs=affy_sample, probe2gene=1)
			}
			if(id_type=="gene") {
				test_sample <- as.vector(CL_DF[CL_DF[,2]==i, 1])
			}	
			containerDF2 <- data.frame(GOID=NULL, NodeSize=NULL, SampleMatch=NULL, Phyper=NULL, Padj=NULL, Term=NULL, Ont=NULL, SampleKeys=NULL, GO_OL_Match=NULL)
			count <- 0
			for(j in gocats) {
				count <- count+1
				GOHyperGAll_result <- GOHyperGAll(gocat=j, sample=test_sample, ...)
				simplifyDF <- GOHyperGAll_Simplify(GOHyperGAll_result, gocat=j, cutoff=cutoff, correct=correct)
				if(length(simplifyDF)==0) { # If simplifyDF() returns empty data frame, then include at least the first two best scoring GO entries 
					simplifyDF <- GOHyperGAll_Simplify(GOHyperGAll_result[1:2,], gocat=j, cutoff=1, correct=T) 
				}
				tempDF <- data.frame(GOHyperGAll_result[GOHyperGAll_result[,1] %in% simplifyDF[,1], ], GO_OL_Match=simplifyDF[,2])
				containerDF2 <- rbind(containerDF2, tempDF)
				if(length(recordSpecGO)>0) {
					containerDF2 <- rbind(containerDF2, data.frame(GOHyperGAll_result[GOHyperGAll_result[,1]==recordSpecGO[count],], GO_OL_Match=GOHyperGAll_result[GOHyperGAll_result[,1]==recordSpecGO[count],3]))
				}
				no_annot <- test_sample[!test_sample %in% eval(parse(text=paste("GO_", j, "_DF", sep="")))[,2]]
				no_annot <- no_annot[no_annot!="no_match"]
				containerDF2 <- rbind(containerDF2, data.frame(GOID=paste("no_annot_", j, sep=""), NodeSize=NA, SampleMatch=length(no_annot), Phyper=NA, Padj=NA, Term=NA, Ont=j, SampleKeys=paste(no_annot, collapse=", "), GO_OL_Match=length(no_annot)))
			}
			tempDF2 <- data.frame(CLID=rep(i, times=length(containerDF2[,1])), CLSZ=rep(unique(as.vector(CL_DF[CL_DF[,2]==i,3])), times=length(containerDF2[,1])), containerDF2)
			containerDF <- rbind(containerDF, tempDF2)
		}
		containerDF <- containerDF[, c(1:9,11,10)]
		return(containerDF)
	}	
}

## Apply GOCluster_Report
## BatchResult <- GOCluster_Report(CL_DF=CL_DF, method="all", id_type="gene", CLSZ=10, cutoff=0.001, gocats=c("MF", "BP", "CC"), recordSpecGO=c("GO:0003674", "GO:0008150", "GO:0005575"))
cat("\n", "(C.3) Batch analysis of many gene clusters: \n \t BatchResult <- GOCluster_Report(CL_DF=CL_DF, method=\"all\", id_type=\"gene\", CLSZ=10, cutoff=0.001, gocats=c(\"MF\", \"BP\", \"CC\"), recordSpecGO=c(\"GO:0003674\", \"GO:0008150\", \"GO:0005575\"))", "\n")


