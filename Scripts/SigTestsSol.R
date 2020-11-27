#20th Nov 2020
#Script take the LRT scores from codeml output and performs a one-tailed, Chisq (df=1) 
#test of significance - solitary lineage version

#Function to run ChiSq on LRT values
pval <- function(x){
  pchisq(x, df = 1, lower.tail = F)
}


#Apply to the list
for (i in 1:length(data.sol)){
  data.sol[[i]]$pvalue <- pval(data.sol[[i]]$LRT)  
}

#Use the Benjamini-Hochburg procedure to correct for multiple testing

for (i in 1:length(data.sol)){
  data.sol[[i]]$adj_pvalue <- p.adjust(data.sol[[i]]$pvalue)  
}