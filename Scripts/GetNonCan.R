#9th November 2021
#Getting a potentially new noncanon list of genes

#I've taken all genes that are significantly (Benjamini-Hochburg corrected pvalues < 0.05)
#up, down or differentially regulated in the metatranscriptomic Doublet paper

prelim <- read.table("ImmResources/Prelim_NonCanList.txt", header = F)
head(prelim)

#remove duplicates and see how many we're left with

prelim <- as.data.frame(prelim[!duplicated(prelim),])
names(prelim) <- "Gene"

#344! Shit.
#Let's see how many have the same LOC ids in the current build

amel.id <- read.table("Genome_Misc/Amel_HAv3.1_MasterIDtable.tsv", header = T, sep = "\t")
head(amel.id)

prelim.2 <- merge(prelim, amel.id, by.x = "Gene", by.y="GeneID", all.x = T)
prelim.2 <- prelim.2[rowSums(is.na(prelim.2)) != ncol(prelim.2), ]

prelim.2 <- prelim.2[order(prelim.2$TranscriptID),]
head(prelim.2)
tail(prelim.2, n = 100)

#A fuck ton of them

#Try again

prelim <- read.table("ImmResources/Prelim_NonCanList2.txt", header = F)
prelim <- as.data.frame(prelim[!duplicated(prelim),])
names(prelim) <- "Gene"

prelim.2 <- merge(prelim, amel.id, by.x = "Gene", by.y="GeneID", all.x = T)
prelim.2 <- prelim.2[rowSums(is.na(prelim.2)) != ncol(prelim.2), ]

prelim.2 <- prelim.2[order(prelim.2$TranscriptID),]
tail(prelim.2, n = 100)

length(prelim.2$Gene[!duplicated(prelim.2$Gene)])

#now I need to remove those that are actually canonically immune
pre.lim3 <- prelim.2[!prelim.2$Gene %in% imm.nov$Gene,]
pre.lim3[pre.lim3$Gene == "LOC413749",]
length(pre.lim3$Gene[!duplicated(pre.lim3$Gene)])
noncangenes <- pre.lim3$Gene[!duplicated(pre.lim3$Gene)]
grep("LOC413749", noncangenes)

nrow(dvg2[dvg2$Gene %in% noncangenes,])
dvg3 <- dvg2
dvg3$Function[dvg3$Gene %in% noncangenes] <- "Non-Canon"
head(dvg3)

write.table(dvg3, "output/AlignmentInfo2_Nov2021.tsv", 
             col.names = T, row.names = F, quote = F, sep = "\t")
