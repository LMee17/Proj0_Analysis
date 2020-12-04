#3rd December 2020
#Script to pull out the most significant genes, aggregate it by the point in which the
#phylogenetic tree codeml was assessing selection (social origin/elaboration), and 
#pull out those that are appear more than once

library("dplyr")

#unlist data
one <- bind_rows(data.sol)
#keep only the significant hits and some other details
two <- one[one$adj_pvalue < 0.5,]
two <- two[,c(1,5)]
#aggregate by social origin
two <- aggregate(two[2], two[-2], 
                     FUN = function(X) paste(unique(X), collapse=", "))
head(two)

#get gene descriptions back in
three <- one[,c(1,6:7)]
three <- three[!duplicated(three$Gene),]

#combine
undersel.sol <- merge(two, three, by = "Gene") 
head(undersel.sol)

#add a way of counting those genes that occur more than once
library("stringr")
undersel.sol$count <- str_count(undersel.sol$SocOrigin, pattern = ",")+1

#create a dataframe of said genes
goi.sol <- undersel.sol[undersel.sol$count > 1,]
