
#25th Aug 2020
#Run bootstrapping code on each origin and assess where the number of canon/noncanon genes
#found as under selection fall agains the resulting distribution of random selections

library("dplyr")


#Canonical first

#This is going to be a function within a function job
#First, a function to sample the population of the origin test (ie all genes that underwent
#branch site analysis in the run) with the number of genes that were found under selection
#for that run selected at random. This will be compared with the list of canon genes as
#classified in the data and the number of "hits" - number of times canoncial genes are
#selected at random - is recorded and returned.
can.sel <- function(x){
  x$Gene <- as.factor(x$Gene)
  x$Class <- as.factor(x$Class)
  pop <- x$Gene
  no <- nrow(x[x$adj_pvalue < 0.05,])
  can <- droplevels(as.factor(x$Gene[x$Class == "Canon"]))
  bite <- sample(pop, size=no,replace=T)
  hit.can <- sum(as.integer(can %in% bite))
  return(hit.can)
}



#Secondly, I need a logical function that will allow for me to check to see if the
#actual number of canon hits falls in or outside of the confidence intervals of the
#distribution produced by the sampling iterations.
call.make <- function(x,y){
  if (x == "TRUE" || y == "TRUE"){
    verdict <- "Reject"
  } else {
    verdict <- "Accept"
  }
  return(verdict)
}

#Thirdly, this function uses the above function to begin to actually run the bootstrap
#analysis. "x" is the target dataframe (or list, using lapply) and "y" is the number of
#iterations required. 
#Once the bootstrap has been iterated y amount of times, a matrix is finished with all
#the canonical gene hits that were found with each iteration. This matrix can then be
#assessed with descriptive statistics - comparing the actual number of canonical genes
#found in the analysis with the mean, upper and lower confidence intervals of the 
#distribution created by the bootstrap should allow me to flag any instances were
#the actual data does not reflect what would be expected by chance.
#The output should be a table or row of a table with SocOrigin (soc), Gene Class (class), number of
#that class under selection in the real data (can.data), mean of how many times that class
#was selected at random after y iterations (mu), the lower and upper confidence intervals (lower/upper),
#a logical response as to whether or not the actual data fall outside of these intervals (verdict),
#and the number of times the actual data frequency occurred as a division of the number
#of iterations (ie 1 occurence out of 100 iterations would be 0.01), (echo).
can.assess <- function(x,y){
  x$SocOrigin <- as.factor(x$SocOrigin)
  out <- rep(0,y)
  for (i in 1:y){out[i] <- can.sel(x)}
  mu <- mean(out)
  soc <- droplevels(as.factor(x$SocOrigin[!duplicated(x$SocOrigin)]))
  can.data <- nrow(subset(x, adj_pvalue < 0.05 & Class == "Canon"))
  class <- "Canon"
  low <- as.integer(quantile(out, probs=0.025))
  high <- as.integer(quantile(out, probs=0.975))
  echo <- sum(out==can.data)/y
  check1 <- mu < low
  check2 <- mu > high
  verdict <- call.make(check1, check2)
  final <- as.data.frame(cbind(as.character(soc), class, can.data, 
                           mu, low, high, 
                           echo, verdict))
  names(final) <- c("SocOrigin", "GeneClass", "RealDataFreqUnderSel", "BootStrapMean", 
                    "LowerCI", "UpperCI", "PropRealDataResult","NullHypothesis")
  return(final)
}

out <- rep(1,5)




can.sol.results <- lapply(X=data.sol, y=100000, FUN = can.assess)
can.sol.results <- bind_rows(can.sol.results)

#Same again, for Noncan.
non.sel <- function(x){
  x$Gene <- as.factor(x$Gene)
  x$Class <- as.factor(x$Class)
  pop <- x$Gene
  no <- nrow(x[x$adj_pvalue < 0.05,])
  can <- droplevels(as.factor(x$Gene[x$Class == "NonCanon"]))
  bite <- sample(pop, size=no,replace=T)
  hit.non <- sum(as.integer(can %in% bite))
  return(hit.non)
}


non.assess <- function(x,y){
  x$SocOrigin <- as.factor(x$SocOrigin)
  out <- rep(0,y)
  for (i in 1:y){out[i] <- non.sel(x)}
  mu <- mean(out)
  soc <- droplevels(as.factor(x$SocOrigin[!duplicated(x$SocOrigin)]))
  non.data <- nrow(subset(x, adj_pvalue < 0.05 & Class == "NonCanon"))
  class <- "NonCanon"
  low <- as.integer(quantile(out, probs=0.025))
  high <- as.integer(quantile(out, probs=0.975))
  echo <- sum(out==non.data)/y
  check1 <- mu < low
  check2 <- mu > high
  verdict <- call.make(check1, check2)
  final <- as.data.frame(cbind(as.character(soc), class, non.data, 
                               mu, low, high, 
                               echo, verdict))
  names(final) <- c("SocOrigin", "GeneClass", "RealDataFreqUnderSel", "BootStrapMean", 
                    "LowerCI", "UpperCI", "PropRealDataResult","NullHypothesis")
  return(final)
}

non.sol.results <- lapply(X=data.sol, y=100000, FUN = non.assess)
non.sol.results <- bind_rows(non.sol.results)

#Create one table
boot.sol.results <- rbind(can.sol.results, non.sol.results)
boot.sol.results

