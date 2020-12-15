#15th December 2020
#A script to pull out genes under selection per lineage per class (noncanon/random) of gene.

#Start with making an object with the name of each lineage 

branches <- unique(as.factor(droplevels(data.all$SocOrigin)))

for (i in 1:length(branches)){
  soc <- branches[i]
  one <- data.all[data.all$SocOrigin == paste(soc) & data.all$Class == "NonCanon",]
  two <- merge(one, iso.verse, by = "Gene")
  
}

one <- data.all[data.all$SocOrigin == "Meli" & data.all$Class == "Random" & 
                  data.all$adj_pvalue < 0.051,]

two <- merge(one, iso.verse, by = "Gene")
head(two)
 

miss <- one[!one$Gene %in% two$Gene,]
miss

grep("LOC100577196", undersel.all)


#needs to be isoform reduced