#4th July 2020
#A script to run Fisher's exact test on unadjusted pvalues

fisher.craft <- function(x){
  x$Class <- as.factor(x$Class)
  d.sig <- x[x$pvalue < 0.051,]
  #ONE: Imm vs Random
  a.sig <- as.integer(summary(d.sig$Class)[1]) + as.integer(summary(d.sig$Class)[2])
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
  ImmVRandom <- fisher.test(one)$p.value
  #TWO: Can versus random
  a.sig <- as.integer(summary(d.sig$Class))[1]
  a.tot <- as.integer(summary(x$Class))[1]
  a.non <- a.tot - a.sig
  Sig <- c(a.sig, b.sig)
  NonSig <- c(a.non, b.non)
  two <- cbind(Sig,NonSig)
  rownames(two) <- c("Can", "Random")
  CanVRandom <- fisher.test(two)$p.value
  #THREE: NonCan versus random
  a.sig <- as.integer(summary(d.sig$Class))[2]
  a.tot <- as.integer(summary(x$Class))[2]
  a.non <- a.tot - a.sig
  Sig <- c(a.sig, b.sig)
  NonSig <- c(a.non, b.non)
  three <- cbind(Sig,NonSig)
  rownames(three) <- c("NonCan", "Random")
  NonCanVRandom <- fisher.test(three)$p.value
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
  CanVNonCan <- fisher.test(four)$p.value
  out <- rbind(ImmVRandom, CanVRandom, NonCanVRandom, CanVNonCan)
  return(out)
}

p.results <- lapply(data, fisher.craft)

#Combine into a table
#Record tests
Test <- rownames(p.results[[1]])

#Unlist and add test names as first column
p.results <- as.data.frame(bind_rows(p.results))
p.results <- cbind(Test, p.results)


#Write up

write.table(p.results, "output/FishersExactTestbyOrigin_unadjustedpvalues.tsv", 
            sep="\t", quote=F, row.names=F)
