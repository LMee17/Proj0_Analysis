#8th November 2021
#Preparing Gene lists for multihit model checks

#I need to pull out genes that are PSGs in at least one lineage
#Remove genes that are too small (less than 30 codons)

library("dplyr")

#Read in a file with all the branch-site test results (pvalues)
alldata <- read.table("output/All.BranchSite.CodeML.Data.tsv",
                      header = T, sep = "\t")


#Read in file with alignment physical info
dvg2 <- read.table("output/AlignmentInfo.tsv",
                  header = T, sep = "\t")

#combine
alldata2 <- merge(alldata, dvg2, by = "Gene", all = T)
head(alldata2)

#Make prettier
alldata3 <- alldata2[,c(1, 14, 17, 16, 13, 18, 2:12)]
head(alldata3)

#Write up (useful file)
write.table(alldata3, "output/All.Data.PSGbyBranch.AlnInfo.Nov2021.tsv",
            col.names = T, row.names = F, quote = F, sep = "\t")

#Now start filtering out
#Remove all genes under 30 codons
#Go through each branch test and keep genes that are undersel
#make empty list for results (size = number of branch site tests)

alldata3 <- alldata3[rowSums(is.na(alldata3)) != ncol(alldata3), ]
alldata3[with(alldata3, order (Codon)), ]
#It says there isn't any but just to make sure

alldata4 <- alldata3[alldata3$Codon > 30,]
alldata4 <- alldata4[rowSums(is.na(alldata4)) != ncol(alldata4), ]

#Pull out PSGs
out <- vector(mode = "list", length = 11)

for (i in 7:17){
  j <- i - 6
  out[[j]] <- alldata4$Gene[alldata4[,i] < 0.05]
  out[[j]] <- out[[j]][!is.na(out[[j]])]
}

out[[1]]


psglist <- unlist(out)
psglist <- unique(psglist)
psglist <- as.data.frame(psglist)
head(psglist)

#write up
write.table(psglist, "output/PSGlist.txt", col.names = F, row.names = F, quote = F)
