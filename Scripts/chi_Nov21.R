#13th November 2021

#A script to run chisq tests looking at each gene class/ function against the background gene set

cont.list <- vector(mode = "list", length = length(data))

one <- data[[5]]
one$PSG <- ifelse(one$adj_pvalue < 0.05, "Under Selection", "Not Under Selection")

table(one$Class, one$PSG)

work
table(work$PSG, work$Function)

x<- data[[7]]
head(trial)
trial$PSG <- ifelse(trial$adj_pvalue < 0.05, "Under Selection", "Not Under Selection")
a <- table(trial$Function, trial$PSG)
chisq.test(a)

chi.craft <- function(x){
  work2 <- x
  work2$Class <- as.factor(work2$Class)
  work2$PSG <- ifelse(work2$adj_pvalue < 0.05, "Under Selection", "Not Under Selection")
  work2 <- work2[!is.na(work2$Class),]
  #ONE: all can versus background
  one <- work2[work2$Class == "Canon Immune" | work2$Class == "Background",]
  one$Class <- droplevels(one$Class)
  cont.1 <- table(one$PSG, one$Class)
  one.res <- chisq.test(cont.1, simulate.p.value = TRUE, B = 10000)
  run <- unique(work2$SocOrigin)
  prop.can <- sprintf(cont.1[[2,2]] / sum(cont.1[,2]), fmt = '%#.3f')
  prop.back1 <- sprintf(cont.1[[1,2]] / sum(cont.1[,1]), fmt = '%#.3f')
  chi1 <- sprintf(one.res$statistic[[1]], fmt = '%#.3f')
  p1 <- sprintf(one.res$p.value, fmt = '%#.3f')
  one.out <- cbind(run, prop.can, prop.back1, chi1, p1)
  #TWO: non-canon versus background
  two <- work2[work2$Class == "Non-Canon" | work2$Class == "Background", ]
  two$Class <- droplevels(two$Class)
  cont.2 <- table(two$PSG, two$Class)
  two.res <- chisq.test(cont.2, simulate.p.value = T, B = 10000)
  prop.non <- sprintf(cont.2[[2,2]] / sum(cont.2[,2]), fmt = '%#.3f')
  prop.back2 <- sprintf(cont.2[[1,2]] / sum(cont.2[,1]), fmt = '%#.3f')
  chi2 <- sprintf(two.res$statistic[[1]], fmt = '%#.3f')
  p2 <- sprintf(two.res$p.value, fmt = '%#.3f')
  two.out <- cbind(prop.non, prop.back2, chi2, p2)
  out <- as.data.frame(cbind(one.out, two.out))
  names(out) <- c("Branch(es) Tested", "CanonProp UnderSel", "BackgroundProp UnderSel",
                  "X-Squared (df = NA)", "p.value", "NonCanonProp UnderSel", "BackgroundProp UnderSel",
                  "X-Squared (df = NA)", "p.value")
  return(out)
}


chi.list <- lapply(data, chi.craft)
chi.list[[5]]

chi.res <- bind_rows(chi.list)
