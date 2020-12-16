#15th December 2020
#A script to pull out genes under selection per lineage per class (noncanon/random) of gene.

#Start with making an object with the name of each lineage and each class

branches <- unique(as.factor(droplevels(data.complete$SocOrigin)))
classes <- cbind("NonCanon", "Random")

#This loop will go through each lineage and extract the genes per class
#keep only those that are significant, remove the rest of the information
#and write a list in input/GOI
#it will also check for genes that get left behind, meaning have no protein
#counterpart. These are likely to be because the codon alignment was originally
#based on a noncoding transcript and need to be removed from the analysis entirely

for (i in 1:length(branches)){
  soc <- branches[i]
  for (i in 1:length(classes)){
    class <- classes[i]
    one <- data.complete[data.complete$SocOrigin == paste(soc) & 
                           data.complete$Class == paste(class),]
    two <- one[one$adj_pvalue < 0.051,]
    three <- merge(two, iso.verse, by = "Gene")
    three <- as.data.frame(three$Resolved_Isoform)
    check <- nrow(two) == nrow(three)
    if (check == FALSE) {
      miss <- as.data.frame(two$Gene[!two$Gene %in% three$Gene])
      miss.file <- paste(soc, class, "missinggene_G0analysis.txt", sep=".")
      write.table(miss, file = paste("output/GO_missingGenes/", miss.file, sep = ""),
                  col.names = F, row.names = F, quote = F)
    } 
    filename <- paste(soc, class,sep = "_")
    write.table(three, file = paste("input/GOI/", filename, ".txt", sep=""), 
                col.names = F, row.names = F, quote = F)
  }
}


