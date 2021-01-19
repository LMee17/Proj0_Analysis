#18th December
#A script to get the GO data from topGO analyses into a readable format by UpSetR

#Each results table exists as a separate .tsv within output/topGO
#I need to have them all in one dataframe, with the column header being the lineage / test
#and the row.names being GO terms.
#Then there would be a 0 / 1 system if that GO term was found to be significantly
#over-represented in those set of genes under selection (1 obviously being present). 

#How the fuck am I going to do this.

library(tidyr)
library(dplyr)

go.base <- read.table("Genome_Misc/Amel_HAv3.1_GOverse.txt", sep = "\t")
head(go.base)

go.base <- go.base %>%
  mutate(V2 = strsplit(V2, ",")) %>%
  unnest(V2)

head(go.base)
go.base <- as.data.frame(go.base$V2)
head(go.base)
go.base <- as.data.frame(go.base[!duplicated(go.base),])
head(go.base)
names(go.base) <- "GOterm"
head(go.base)

#I think may be the easiest way of doing this would be to store each lineage's enriched 
#GOterms into an object, then make a logic argument based on presence / absence of the 
#GO term of that row (go.all.fin$GOterm) using the %in% function

#read in files into a list 
filelist <- list.files(path = "output/topGO/", pattern = "*tsv")

#write function to extract run name
getname2 <- function(x){
  one <- strsplit(x, "_GOterms", fixed =T)
  name <- one[[1]][1]
  return(name)
}

#read files into list
go.list <- lapply(filelist, function(x){
  read.table(file = paste("output/topGO/",x, sep=""), header = T, sep ="\t", quote ="")
})

#run through the filelist, in the same order as go.list, extract the run name
#using the getname2 function, extract the significant GOterms from that run,
go.terms <- go.base

for (i in 1:length(filelist)){
  nombre <- getname2(filelist[i])
  go.sig <- go.list[[i]]$GO.ID
  go.terms[,(i+1)] <- ifelse(go.terms$GOterm %in% go.sig, paste(1), paste(0))
  names(go.terms)[i+1] <- paste(nombre)
}

for (i in 1:length(filelist)){
  go.terms[,(i+1)] <- as.numeric(go.terms[,(i+1)])
}

#split into immunity and random
imm.terms <- go.terms[,c(1,2,4,6,8,10,12:13,15,17,19,21,22,24,26,27,29)]

ran.terms <- go.terms[,c(1,3,5,7,9,11,14,16,18,20,23,25,28,30)]


#it would also be prudent to separate these into socialities. 
social <- c("AllOrigin", "CorbSoc", "Halictid", "Xylo")
complex <- c("AllComplex", "Apis", "Meli")
solitary <- c("AllSol", "Habro", "Nova", 'Mega')

getname3 <- function(x){
  one <- strsplit(x, "_", fixed = T)
  name <- one[[1]][1]
  return(name)
}

#a function to separate immune from random
class.craft <- function(x){
  one <- strsplit(x, "_", fixed = T)
  name <- one[[1]][2]
  return(name)
}

#set bases of the sociality table (v1 = all go terms)

soc.terms <- go.base

#set iteration timers to 0
soc.imm.i <- 0
soc.ran.i <- 0
com.imm.i <- 0
com.ran.i <- 0
sol.imm.i <- 0
sol.ran.i <- 0

#this bit of code runs through the filelist and first sorts the iteration into immune or
#random. It then checks if it is in any of the three
#sociality vectors as defined above. If there is a hit, the timer for that sociality is 
#increased by 1. If its the first of that sociality, a new object is made to record
#the significantly enriched GO terms. If it is a later iteration, the next set of GOterms
#is added to that which already exists. 
for (i in 1:length(filelist)){
  class <- class.craft(filelist[i])
  run <- getname3(filelist[i])
  if (class == "Random"){
    if (run %in% social){
      soc.ran.i <- soc.ran.i +1
      if (soc.ran.i > 1){
        go.sig.soc.ran <- c(go.sig.soc.ran, go.list[[i]]$GO.ID)
        soc.check.ran <- c(soc.check.ran, paste(run))
      } else {
        go.sig.soc.ran <- go.list[[i]]$GO.ID
        soc.check.ran <- paste(run)
      }
    }
    if (run %in% complex){
      com.ran.i <- com.ran.i + 1
      if (com.ran.i > 1){
        go.sig.com.ran <- c(go.sig.com.ran, go.list[[i]]$GO.ID)
        com.check.ran <- c(com.check.ran, paste(run))
      } else {
        go.sig.com.ran <- go.list[[i]]$GO.ID
        com.check.ran <- paste(run)
      }
    }
    if (run %in% solitary){
      sol.ran.i <- sol.ran.i + 1
      if (sol.ran.i > 1){
        go.sig.sol.ran <- c(go.sig.sol.ran, go.list[[i]]$GO.ID)
        sol.check.ran <- c(sol.check.ran, paste(run))
      } else {
        go.sig.sol.ran <- go.list[[i]]$GO.ID
        sol.check.ran <- paste(run)
      }
    }
  } else {
    if (run %in% social){
      soc.imm.i <- soc.imm.i +1
      if (soc.imm.i > 1){
        go.sig.soc.imm <- c(go.sig.soc.imm, go.list[[i]]$GO.ID)
        soc.check.imm <- c(soc.check.imm, paste(run))
      } else {
        go.sig.soc.imm <- go.list[[i]]$GO.ID
        soc.check.imm <- paste(run)
      }
    }
    if (run %in% complex){
      com.imm.i <- com.imm.i + 1
      if (com.imm.i > 1){
        go.sig.com.imm <- c(go.sig.com.imm, go.list[[i]]$GO.ID)
        com.check.imm <- c(com.check.imm, paste(run))
      } else {
        go.sig.com.imm <- go.list[[i]]$GO.ID
        com.check.imm <- paste(run)
      }
    }
    if (run %in% solitary){
      sol.imm.i <- sol.imm.i + 1
      if (sol.imm.i > 1){
        go.sig.sol.imm <- c(go.sig.sol.imm, go.list[[i]]$GO.ID)
        sol.check.imm <- c(sol.check.imm, paste(run))
      } else {
        go.sig.sol.imm <- go.list[[i]]$GO.ID
        sol.check.imm <- paste(run)
      }
    }
  }
}


#remove duplications. Probably not necessary but good to have a way to check numbers
#later
go.sig.sol.imm <- go.sig.sol.imm[!duplicated(go.sig.sol.imm)]
go.sig.soc.imm <- go.sig.soc.imm[!duplicated(go.sig.soc.imm)]
go.sig.com.imm <- go.sig.com.imm[!duplicated(go.sig.com.imm)]
go.sig.sol.ran <- go.sig.sol.ran[!duplicated(go.sig.sol.ran)]
go.sig.soc.ran <- go.sig.soc.ran[!duplicated(go.sig.soc.ran)]
go.sig.com.ran <- go.sig.com.ran[!duplicated(go.sig.com.ran)]

#Make a table with everything, by column (unfortunately, this is so against my coding ethics)
soc.terms$Canon_input <- go.terms$Canon_input
soc.terms$NonCan_input <- go.terms$NonCan_input
soc.terms$Hym_Noncan <- go.terms$Hym_NonCan
soc.terms$Solitary_Imm <- ifelse(soc.terms$GOterm %in% go.sig.sol.imm, paste(1), paste(0))
soc.terms$Social_Imm <- ifelse(soc.terms$GOterm %in% go.sig.soc.imm, paste(1), paste(0))
soc.terms$Complex_Imm <- ifelse(soc.terms$GOterm %in% go.sig.com.imm, paste(1), paste(0))
soc.terms$Solitary_Random <- ifelse(soc.terms$GOterm %in% go.sig.sol.ran, paste(1), paste(0))
soc.terms$Social_Random <- ifelse(soc.terms$GOterm %in% go.sig.soc.ran, paste(1), paste(0))
soc.terms$Complex_Random <- ifelse(soc.terms$GOterm %in% go.sig.com.ran, paste(1), paste(0))

for (i in 2:ncol(soc.terms)){
  soc.terms[,i] <- as.numeric(soc.terms[,i])
}

#Immune table - include noncan and can input + all noncan runs, separated by sociality
soc.imm.terms <- soc.terms[,c(1:7)]

#Random Table
soc.ran.terms <- soc.terms[,c(1,8:10)]

#Write Up
#All lineages, all classes
write.table(go.terms, "output/upsetr_tables/AllLineage_AllClass_goUpSetR.tsv", sep = "\t", quote = F,
            row.names = F)
#All lineages, immune
write.table(imm.terms, "output/upsetr_tables/AllLineage_Immune_goUpSetR.tsv", 
            sep = "\t", quote = F, row.names = F)
#All lineages, random
write.table(ran.terms, "output/upsetr_tables/AllLineage_Random_goUpSetR.tsv", 
            sep = "\t", quote = F, row.names = F)
#All Socialities, all classes
write.table(soc.terms, "output/upsetr_tables/AllSocialities_AllClass_goUpSetR.tsv",
            sep = "\t", quote = F, row.names = F)
#All Socialities, immune
write.table(soc.imm.terms, "output/upsetr_tables/AllSocialities_Immune_goUpSetR.tsv",
            sep = "\t", quote = F, row.names = F)
#All Socialities, random
write.table(soc.ran.terms, "output/upsetr_tables/AllSocialities_Random_goUpSetR.tsv",
            sep = "\t", quote = F, row.names = F)
