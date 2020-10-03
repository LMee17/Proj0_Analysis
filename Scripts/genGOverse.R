#17th September 2020

#The purpose of this script is simply to get an eggNOG txt file into an annotation
#file that is readable by topGO and can thus be used to set the "geneverse" in a
#GO term analysis.

#read in the file

raw <- read.table("Genome_Misc/AmelHAv3.1_eggNOG.txt", sep="\t", quote="", header=F)
head(raw)

#all I need is an id and GO terms

raw <- raw[,c(1,7)]
head(raw)

#this file has GO terms per protein ID, and so I will need to link gene and protein IDs
#I reduced all isoforms before beginning the analysis and so I need to use the 
#reduced table
#there shouldn't be any issues as I used the reduced proteome in the eggNOG annotation

reduced <- read.table("Genome_Misc/AmelHAv3.1_filtered.checked.faa.table", sep = "\t", header = T, 
                      quote = "")

head(reduced)

#combine

ann <- merge(reduced, raw, by.x = "Resolved_Isoform", by.y = "V1", all = T)
head(ann)

#check we haven't got any problems with duplicated genes
ann[duplicated(ann$Gene),]

#so these are an artefact from an issue early on making the isoform reduced table
#the () being present in the gene name has led to all proteins being brought through
#I will just remove for now

remove <- ann[duplicated(ann$Gene),]
remove <- remove$Resolved_Isoform
remove

nrow(ann[!ann$Resolved_Isoform %in% remove,])
ann2 <- ann[!ann$Resolved_Isoform %in% remove,]

#double check
ann2[duplicated(ann2$Gene),]

#remove genes
ann2 <- ann2[,c(1,3)]
head(ann2)

#good to go
write.table(ann2, "Genome_Misc/Amel_HAv3.1_GOverse.txt", sep = "\t", quote = F, 
            col.names = F,row.names = F)

