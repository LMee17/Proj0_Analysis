#3rd July 2020
#Script to read immunity codeml branch site analysis results and produce output tables
#of significance, gene descriptions, phylogenetic branch analysed, etc

#Required
library(dplyr)

#####Read input immunity files######
#This will require temporarily changing the working directory to the input folder where the 
#data is stored.
setwd("input/immunity/")

#List the file names
filelist <- list.files(pattern ="*.txt")

#Use this list to read in all files together and store as a data variable
data <- lapply(filelist, function(x){
  read.table(x, header = T, sep = "\t")
})

#Little function to extract the run name from the lnLResults file and store it as a variable 
#within the dataframe
getname <- function(x){
  one <- strsplit(x, ".", fixed =T)
  name <- one[[1]][1]
  return(name)
}

#Apply this function to the 3 dataframes using the original filelist
for (i in 1:length(data)){
  data[[i]]$SocOrigin <- getname(filelist[i])
  data[[i]]<-na.omit(data[[i]])
  data[[i]]<-data[[i]][!duplicated(data[[i]]$Gene),]
}

######Add details######
#Immune classifications
imm.verse.slim <- imm.verse[,c(1,5)]
data <- lapply(data, function(x) merge(x, imm.verse.slim, by.x="Gene", by.y="GeneID"))

#Gene details
gene.verse.slim <- gene.verse[,c(1,3)]
gene.verse.slim <- gene.verse.slim[!duplicated(gene.verse.slim$GeneID),]
data <- lapply(data, function(x) merge(x, gene.verse.slim, by.x="Gene", by.y="GeneID"))

