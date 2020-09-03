#26th August 
#A script to read in omega scores and annotate with gene class.

#Read in omega data
omega.raw <- read.table("input/omega/dNdSratio.fullcomplement.tsv", header = T, sep = "\t")
#Remove duplicates
omega.edit <- omega.raw[!duplicated(omega.raw$Gene),]
#Remove any NA values
omega.edit <- na.omit(omega.edit)


#Annotate with gene class from data
data.all <- bind_rows(data)
data.class <- data.all[,c(1,6)]
data.class <- data.class[!duplicated(data.class$Gene),]

#remove any erroneous figures
omega <- merge(omega.edit, data.class, by = "Gene")
omega <- omega[!omega$dN.dS > 1,]
head(omega)

#Store any missing genes
omega.error <- data.class[omega$dN.dS > 1, ]
omega.miss <- data.class[!data.class$Gene %in% omega$Gene,]

