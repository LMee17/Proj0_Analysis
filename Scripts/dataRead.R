#12th Nov 2021
#Script to read  codeml branch site analysis results and produce output tables
#of significance, gene descriptions, phylogenetic branch analysed, etc

#Required
library(dplyr)

#####Read input  files######

#List the file names
filelist <- list.files(pattern ="*.txt", path = "input/codemlFull/")

#Use this list to read in all files together and store as a data variable
data <- lapply(filelist, function(x){
  read.table(paste("input/codemlFull/", x, sep =""), header = T, sep = "\t")
})


#Little function to extract the run name from the lnLResults file and store it as a variable 
#within the dataframe
getname <- function(x){
  one <- strsplit(x, ".", fixed =T)
  name <- one[[1]][1]
  return(name)
}



#Apply this function to the dataframes using the original filelist
for (i in 1:length(data)){
  data[[i]]$SocOrigin <- getname(filelist[i])
  data[[i]]<-na.omit(data[[i]])
  data[[i]]<-data[[i]][!duplicated(data[[i]]$Gene),]
}

head(data[[7]])

######Add details######
data <- lapply(data, function(x) merge(x, verse, by= "Gene" , all.x = T))
