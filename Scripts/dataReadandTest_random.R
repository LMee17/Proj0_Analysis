#Read in all results files (saved in input folder, extension lnLResults.txt)
#This will require temporarily changing the working directory to the input folder where 
#the data is stored.
setwd("input/random/")

#List the file names
filelist <- list.files(pattern ="*.txt")

#Use this list to read in all files together and store as a data variable
data <- lapply(filelist, function(x){
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
for (i in 1:length(data)){
  data[[i]]$SocOrigin <- getname(filelist[i])
  data[[i]]<-na.omit(data[[i]])
  data[[i]]<-data[[i]][!duplicated(data[[i]]$Gene),]
}

##Remove "Random" from SocOrigin variable
for (i in 1:length(data)){
  data[[i]]$SocOrigin <- gsub("Random", "", data[[i]]$SocOrigin)
}


#Classify as random
for (i in 1:length(data)){
  data[[i]]$Class <- paste("Random")
}

###Gene Details

data <- lapply(data, function(x) merge(x, gene.verse.slim, by.x="Gene", by.y="GeneID"))

##Tests for Significance

###pvalues

#Determine which genes are significant using ChiSq, one degree of freedom on the likelihood 
#ratio scores already present.


#Function to run ChiSq on LRT values
pval <- function(x){
  pchisq(x, df = 1, lower.tail = F)
}

#Apply to the list
for (i in 1:length(data)){
  data[[i]]$pvalue <- pval(data[[i]]$LRT)  
}


###Adjust pvalue

#Use Benjamini-Hochburg procedure to counteract FDR

for (i in 1:length(data)){
  data[[i]]$adj_pvalue <- p.adjust(data[[i]]$pvalue)  
}


###Remove non-coding elements

ran.nc <- read.table("../../Genome_Misc/Random_noncoding_list.txt")
gene<-strsplit(as.character(ran.nc$V1),".",fixed = T)
ran.nc$Gene<-sapply(gene,"[",1)
ran.nc$V1 <- NULL
ran.nc$Gene <- as.factor(ran.nc$Gene)

for (i in 1:length(data)){
  data[[i]]$Gene <- as.factor(data[[i]]$Gene)
  data[[i]] <- data[[i]][!data[[i]]$Gene %in% ran.nc$Gene,]
}


##Write Up

data.ran <- bind_rows(data)


setwd("~/Documents/Projects/ProjectZero/Proj0_Analysis/")
write.table(data.ran, "output/CodeML_RandomResults_All.tsv", sep = "\t", row.names = F, 
            quote = F)
