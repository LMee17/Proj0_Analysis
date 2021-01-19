#18th December
#A script to get the genes under selection from the paml analyses into a readable format 
#by UpSetR


library(tidyr)
library(dplyr)

gene.base <- as.factor(gene.verse$GeneID[!duplicated(gene.verse$GeneID)])
gene.base <- as.data.frame(gene.base)
names(gene.base) <- "GeneID"
head(gene.base)

soc <- data.complete$SocOrigin
soc <- soc[!duplicated(soc)]

gene.undersel <- gene.base
j <- 1


for (i in 1:length(soc)){
  nombre <- soc[i]
  j <- j + 1
  gene.sig.can <- data.complete$Gene[data.complete$SocOrigin == paste(nombre)
                                & data.complete$adj_pvalue < 0.05
                                & data.complete$Class == "Canon"]  
  gene.undersel[,j] <- ifelse(gene.undersel$GeneID %in% gene.sig.can, paste(1), paste(0))
  names(gene.undersel)[j] <- paste(nombre, "_Canon", sep = "")
  j <- j + 1
  gene.sig.non <- data.complete$Gene[data.complete$SocOrigin == paste(nombre)
                                & data.complete$adj_pvalue < 0.05
                                & data.complete$Class == "NonCanon"]
  gene.undersel[,j] <- ifelse(gene.undersel$GeneID %in% gene.sig.non, paste(1), paste(0))
  names(gene.undersel)[j] <- paste(nombre, "_NonCanon", sep = "")
  j <- j + 1
  gene.sig.ran <- data.complete$Gene[data.complete$SocOrigin == paste(nombre)
                                & data.complete$adj_pvalue < 0.05
                                & data.complete$Class == "Random"]
  gene.undersel[,j] <- ifelse(gene.undersel$GeneID %in% gene.sig.ran, paste(1), paste(0))
  names(gene.undersel)[j] <- paste(nombre, "_Random", sep = "")
}

head(gene.undersel)

write.table(gene.undersel, "output/upsetr_tables/GenesUnderSel_All.tsv", 
            sep = "\t", quote = F, row.names = F)

