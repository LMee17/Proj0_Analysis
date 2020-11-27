#Read in all results files (saved in input folder, extension lnLResults.txt)
#This will require temporarily changing the working directory to the input folder where 
#the data is stored.
setwd("~/Documents/Projects/ProjectZero/Proj0_Analysis/input/solitary/random/")

#List the file names
filelist <- list.files(pattern ="*.txt")

#Use this list to read in all files together and store as a data variable
data.sol <- lapply(filelist, function(x){
  read.table(x, header = T, sep = "\t")
})

#Little function to extrac the run name from the lnLResults file and store it as a variable 
#within the dataframe
getname <- function(x){
  one <- strsplit(x, ".", fixed =T)
  name <- one[[1]][1]
  return(name)
}

#Apply this function to the 3 dataframes using the original filelist and remove NA rows
for (i in 1:length(data.sol)){
  data.sol[[i]]$SocOrigin <- getname(filelist[i])
  data.sol[[i]]<-na.omit(data.sol[[i]])
}

##Remove "Random" from SocOrigin variable
for (i in 1:length(data.sol)){
  data.sol[[i]]$SocOrigin <- gsub("Random", "", data.sol[[i]]$SocOrigin)
}


#Classify as random
for (i in 1:length(data.sol)){
  data.sol[[i]]$Class <- paste("Random")
}

###Gene Details

data.sol <- lapply(data.sol, function(x) merge(x, gene.verse.slim, by.x="Gene", by.y="GeneID"))

##Tests for Significance

###pvalues

#Determine which genes are significant using ChiSq, one degree of freedom on the likelihood 
#ratio scores already present.


#Function to run ChiSq on LRT values
pval <- function(x){
  pchisq(x, df = 1, lower.tail = F)
}

#Apply to the list
for (i in 1:length(data.sol)){
  data.sol[[i]]$pvalue <- pval(data.sol[[i]]$LRT)  
}


###Adjust pvalue

#Use Benjamini-Hochburg procedure to counteract FDR

for (i in 1:length(data.sol)){
  data.sol[[i]]$adj_pvalue <- p.adjust(data.sol[[i]]$pvalue)  
}


##Write Up

data.sol.ran <- bind_rows(data.sol)


setwd("../../../")
write.table(data.sol.ran, "output/CodeML_RandomResults_Sol_All.tsv", sep = "\t", row.names = F, 
            quote = F)
