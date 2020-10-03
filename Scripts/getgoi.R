#4th September 2020
#Script to pull out the most significant genes, aggregate it by the point in which the
#phylogenetic tree codeml was assessing selection (social origin/elaboration), and 
#pull out those that are appear more than once

library("dplyr")

#unlist data
one <- bind_rows(data)
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
undersel <- merge(two, three, by = "Gene") 
head(undersel)

#add a way of counting those genes that occur more than once
library("stringr")
undersel$count <- str_count(undersel$SocOrigin, pattern = ",")+1

#create a dataframe of said genes
goi <- undersel[undersel$count > 1,]
