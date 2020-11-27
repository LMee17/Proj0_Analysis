chi.craft.adj <- function(x){
  d.sig <- x[x$adj_pvalue < 0.051,]
  #ONE: Imm vs Random
  a.sig <- as.integer(summary(d.sig$Class)[1] + summary(d.sig$Class)[2])
  b.sig <- as.integer(summary(d.sig$Class)[3])
  summary(d.sig$Class)
  a.tot <- as.integer(summary(x$Class)[1]+ summary(x$Class)[2])
  b.tot <- as.integer(summary(x$Class)[3])
  a.non <- a.tot - a.sig
  b.non <- b.tot - b.sig
  Sig <- c(a.sig, b.sig)
  NonSig <- c(a.non, b.non)
  one <- cbind(Sig,NonSig)
  rownames(one) <- c("Imm", "Random")
  ImmVRandom <- chisq.test(one, simulate.p.value = TRUE, B = 10000)$p.value
  #TWO: Can versus random
  a.sig <- as.integer(summary(d.sig$Class))[1]
  a.tot <- as.integer(summary(x$Class))[1]
  a.non <- a.tot - a.sig
  Sig <- c(a.sig, b.sig)
  NonSig <- c(a.non, b.non)
  two <- cbind(Sig,NonSig)
  rownames(two) <- c("Can", "Random")
  CanVRandom <- chisq.test(two, simulate.p.value = TRUE, B = 10000)$p.value
  #THREE: NonCan versus random
  a.sig <- as.integer(summary(d.sig$Class))[2]
  a.tot <- as.integer(summary(x$Class))[2]
  a.non <- a.tot - a.sig
  Sig <- c(a.sig, b.sig)
  NonSig <- c(a.non, b.non)
  three <- cbind(Sig,NonSig)
  rownames(three) <- c("NonCan", "Random")
  NonCanVRandom <- chisq.test(three, simulate.p.value = TRUE, B = 10000)$p.value
  #FOUR: Can vs NonCan
  b.sig <- a.sig
  b.tot <- a.tot
  b.non <- a.non
  a.sig <- as.integer(summary(d.sig$Class))[1]
  a.tot <- as.integer(summary(x$Class))[1]
  a.non <- a.tot - a.sig
  Sig <- c(a.sig, b.sig)
  NonSig <- c(a.non, b.non)
  four <- cbind(Sig,NonSig)
  rownames(four) <- c("Can", "NonCan")
  CanVNonCan <- chisq.test(four, simulate.p.value = TRUE, B = 10000)$p.value
  out <- rbind(ImmVRandom, CanVRandom, NonCanVRandom, CanVNonCan)
  return(out)
}

chi.sol.adj.results <- lapply(data.sol, chi.craft.adj)

#Combine into a table
#Record tests
Test <- rownames(chi.sol.adj.results[[1]])

#Unlist and add test names as first column
chi.sol.adj.results <- as.data.frame(bind_rows(chi.sol.adj.results))
chi.sol.adj.results <- cbind(Test, chi.sol.adj.results)


#Write up

write.table(chi.sol.adj.results, "output/ChiSqTestbyOrigin_adjustedpvalues.tsv", 
            sep="\t", quote=F, row.names=F)
