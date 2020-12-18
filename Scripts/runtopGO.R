#30th September 2020
#A script to read in gene lists and run GO term analysis using TopGO

#LOAD
#BiocManager::install("topGO")
library("topGO")
#BiocManager::install("Rgraphviz")
#library("Rgraphviz")

#Read in files
setwd("input/GOI/")
filelist <- list.files(pattern ="*.txt")

#Use this list to read in all files together and store as a data variable
go.data <- lapply(filelist, function(x){
  read.table(x, header = F, sep = "\t")
})

#Little function to extract the run name from the lnLResults file and store it as a variable 
#within the dataframe
getname <- function(x){
  one <- strsplit(x, ".", fixed =T)
  name <- one[[1]][1]
  return(name)
}

for (i in 1:length(go.data)){
  names(go.data[[i]]) <- getname(filelist[i])
}


#Now to run the tests

#Map GOterms for whole transcriptome
setwd("../..")
geneID2GO<-readMappings("Genome_Misc/Amel_HAv3.1_GOverse.txt")

#Determine gene universe for analysis against genes under selection
geneUniverse<-names(geneID2GO)

#set function for determining biological process orthology terms
BP.craft <- function(x){
  one <- as.data.frame(x)
  goi <- as.character(one[,1])
  geneList <- factor(as.integer(geneUniverse %in% goi))
  names(geneList) <- geneUniverse
  myGOdata<-new("topGOdata",description="AmelGO BP", ontology="BP",
                allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)
  resultClassic<-runTest(myGOdata,algorithm = "classic",statistic = "fisher")
  resultElim <- runTest(myGOdata,algorithm = "elim", statistic = "fisher")
  results <- GenTable(myGOdata, Fisher.elim = resultElim, 
                      Fisher.classic = resultClassic,
                      orderBy = "Fisher.classic" , topNodes = 200)
  results$Ontology <- paste("Biological Process")
  results <- results[results$Fisher.elim < 0.0501,]
}

#apply

bp.list <- as.list(lapply(go.data, BP.craft))

#setwd("output/topGO/")
#for (i in 1:length(bp.list)){
#  filename <- paste(getname(filelist[i]), "_GOterms_BP.tsv", sep="")
#  write.table(bp.list[[i]], filename, col.names = T, row.names = F, sep = "\t", quote = F)
#}

#set function for determining cellular component orthology terms
CC.craft <- function(x){
  one <- as.data.frame(x)
  goi <- as.character(one[,1])
  geneList <- factor(as.integer(geneUniverse %in% goi))
  names(geneList) <- geneUniverse
  myGOdata<-new("topGOdata",description="AmelGO CC", ontology="CC",
                allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)
  resultClassic<-runTest(myGOdata,algorithm = "classic",statistic = "fisher")
  resultElim <- runTest(myGOdata,algorithm = "elim", statistic = "fisher")
  results <- GenTable(myGOdata, Fisher.elim = resultElim, 
                      Fisher.classic = resultClassic,
                      orderBy = "Fisher.classic" , topNodes = 200)
  results$Ontology <- paste("Cellular Component")
  results <- results[results$Fisher.elim < 0.0501,]
}

#apply

cc.list <- as.list(lapply(go.data, CC.craft))

#for (i in 1:length(cc.list)){
#  filename <- paste(getname(filelist[i]), "_GOterms_CC.tsv", sep="")
#  write.table(cc.list[[i]], filename, col.names = T, row.names = F, sep = "\t", quote = F)
#}

#set function for determining molecular function orthology terms
MF.craft <- function(x){
  one <- as.data.frame(x)
  goi <- as.character(one[,1])
  geneList <- factor(as.integer(geneUniverse %in% goi))
  names(geneList) <- geneUniverse
  myGOdata<-new("topGOdata",description="AmelGO MF", ontology="MF",
                allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)
  resultClassic<-runTest(myGOdata,algorithm = "classic",statistic = "fisher")
  resultElim <- runTest(myGOdata,algorithm = "elim", statistic = "fisher")
  results <- GenTable(myGOdata, Fisher.elim = resultElim, 
                      Fisher.classic = resultClassic,
                      orderBy = "Fisher.classic" , topNodes = 200)
  results$Ontology <- paste("Molecular Function")
  results <- results[results$Fisher.elim < 0.0501,]
}

#apply

mf.list <- as.list(lapply(go.data, MF.craft))


#for (i in 1:length(mf.list)){
#  filename <- paste(getname(filelist[i]), "_GOterms_MF.tsv", sep="")
#  write.table(mf.list[[i]], filename, col.names = T, row.names = F, sep = "\t", quote = F)
#}

setwd("~/Documents/Projects/ProjectZero/Proj0_Analysis/")

for (i in 1:length(go.data)){
  bp <- as.data.frame(bp.list[i])
  mf <- as.data.frame(mf.list[i])
  cc <- as.data.frame(cc.list[i])
  out <- rbind(bp, mf, cc)
  filename <- paste(getname(filelist[i]), "_GOterms.tsv", sep="")
  write.table(out, file = paste("output/topGO/", filename, sep = ""), 
              col.names = T, row.names = F, sep = "\t", quote = F)
}



