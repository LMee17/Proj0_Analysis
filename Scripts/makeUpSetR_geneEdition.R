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

#Separate by sociality
social <- c("AllOrigin", "CorbSoc", "Halictid", "Xylo")
complex <- c("AllComplex", "Apis", "Meli")
solitary <- c("AllSol", "Habro", "Nova", 'Mega')


#obv there is no need to have the entire gene base for the canon or noncanon genes here
#so I will separate the gene base into the classes.
#this method also drops the "Epi" runs
can.base <- data.complete$Gene[data.complete$Class == "Canon"]
can.base <- can.base[!duplicated(can.base)]
can.gus <- as.data.frame(can.base)
names(can.gus) <- "GeneID"

non.base <- data.complete$Gene[data.complete$Class == "NonCanon"]
non.base <- non.base[!duplicated(non.base)]
non.gus <- as.data.frame(non.base)
names(non.gus) <- "GeneID"

ran.base <- data.complete$Gene[data.complete$Class == "Random"]
ran.base <- ran.base[!duplicated(ran.base)]
ran.gus <- as.data.frame(ran.base)
names(ran.gus) <- "GeneID"
head(ran.gus)

for (i in 1:length(soc)){
  nombre <- soc[i]
  if (nombre %in% solitary){
    can <- as.character(data.complete$Gene[data.complete$SocOrigin == paste(nombre)
                                       & data.complete$adj_pvalue < 0.05
                                       & data.complete$Class == "Canon"])
    non <- as.character(data.complete$Gene[data.complete$SocOrigin == paste(nombre)
                              & data.complete$adj_pvalue < 0.05
                              & data.complete$Class == "NonCanon"])
    ran <- as.character(data.complete$Gene[data.complete$SocOrigin == paste(nombre)
                              & data.complete$adj_pvalue < 0.05
                              & data.complete$Class == "Random"])
    if (exists("sol.can") == TRUE){
      sol.can <- c(sol.can, can)
      sol.non <- c(sol.non, non)
      sol.ran <- c(sol.ran, ran)
    } else {
      sol.can <- can
      sol.non <- non
      sol.ran <- ran
    }
  }
  if (nombre %in% social){
    can <- as.character(data.complete$Gene[data.complete$SocOrigin == paste(nombre)
                              & data.complete$adj_pvalue < 0.05
                              & data.complete$Class == "Canon"])
    non <- as.character(data.complete$Gene[data.complete$SocOrigin == paste(nombre)
                              & data.complete$adj_pvalue < 0.05
                              & data.complete$Class == "NonCanon"])
    ran <- as.character(data.complete$Gene[data.complete$SocOrigin == paste(nombre)
                              & data.complete$adj_pvalue < 0.05
                              & data.complete$Class == "Random"])
    if (exists("soc.can") == TRUE){
      soc.can <- c(soc.can, can)
      soc.non <- c(soc.non, non)
      soc.ran <- c(soc.ran, ran)
    } else {
      soc.can <- can
      soc.non <- non
      soc.ran <- ran
    }
  }
  if (nombre %in% complex){
    can <- as.character(data.complete$Gene[data.complete$SocOrigin == paste(nombre)
                              & data.complete$adj_pvalue < 0.05
                              & data.complete$Class == "Canon"])
    non <- as.character(data.complete$Gene[data.complete$SocOrigin == paste(nombre)
                              & data.complete$adj_pvalue < 0.05
                              & data.complete$Class == "NonCanon"])
    ran <- as.character(data.complete$Gene[data.complete$SocOrigin == paste(nombre)
                              & data.complete$adj_pvalue < 0.05
                              & data.complete$Class == "Random"])
    if (exists("com.can") == TRUE){
      com.can <- c(com.can, can)
      com.non <- c(com.non, non)
      com.ran <- c(com.ran, ran)
    } else {
      com.can <- can
      com.non <- non
      com.ran <- ran
    }
  }
}


com.can <- com.can[!duplicated(com.can)]
sol.can <- sol.can[!duplicated(sol.can)]
soc.can <- soc.can[!duplicated(soc.can)]
com.non <- com.non[!duplicated(com.non)]
sol.non <- sol.non[!duplicated(sol.non)]
soc.non <- soc.non[!duplicated(soc.non)]
com.ran <- com.ran[!duplicated(com.ran)]
sol.ran <- sol.ran[!duplicated(sol.ran)]
soc.ran <- soc.ran[!duplicated(soc.ran)]


can.gus$Solitary <- ifelse(can.gus$GeneID %in% sol.can, paste(1), paste(0))
can.gus$Social <- ifelse(can.gus$GeneID %in% soc.can, paste(1), paste(0))
can.gus$Advanced <- ifelse(can.gus$GeneID %in% com.can, paste(1), paste(0))

non.gus$Solitary <- ifelse(non.gus$GeneID %in% sol.non, paste(1), paste(0))
non.gus$Social <- ifelse(non.gus$GeneID %in% soc.non, paste(1), paste(0))
non.gus$Advanced <- ifelse(non.gus$GeneID %in% com.non, paste(1), paste(0))

ran.gus$Solitary <- ifelse(ran.gus$GeneID %in% sol.ran, paste(1), paste(0))
ran.gus$Social <- ifelse(ran.gus$GeneID %in% soc.ran, paste(1), paste(0))
ran.gus$Advanced <- ifelse(ran.gus$GeneID %in% com.ran, paste(1), paste(0))

for (i in 2:ncol(can.gus)){
  can.gus[,i] <- as.numeric(can.gus[,i])
}
sum(can.gus$Solitary)
sum(can.gus$Social)
sum(can.gus$Advanced)

for (i in 2:ncol(non.gus)){
  non.gus[,i] <- as.numeric(non.gus[,i])
}
sum(non.gus$Solitary)
sum(non.gus$Social)
sum(non.gus$Advanced)

for (i in 2:ncol(ran.gus)){
  ran.gus[,i] <- as.numeric(ran.gus[,i])
}
sum(ran.gus$Solitary)
sum(ran.gus$Social)
sum(ran.gus$Advanced)

##I'm going to make other versions that have omega values attached for shits and gigs

test.non <- merge(non.gus, omega, by.x="GeneID", by.y="Gene")
test.non

