##Can Vs Non Vs Random within Lineages

#I've needed a table like this for a while so let's do it and add some statistics in
#Gonna have a table with number of each gene tested, how many are under selection
#and how many aren't per class, per lineage. I will then also run prop.tests if I can
#within the dataframe ... if not make a function where I just have the pvalue pasted in
#Then there will be a column with adj pvalues cos there be some multiple testing. 
#Not sure if its ok to be doing this with the first method, where I had separate
#tests per can versus random etc. I will do this way then drop Seth a message over whether
#its appropriate.

soc <- data.complete$SocOrigin[!duplicated(data.complete$SocOrigin)]
soc <- soc[c(-3,-6)]
soc

out <- data.frame(matrix(0, ncol = 6, nrow = length(soc)))
colnames(out) <- c("CanonUnderSel", "NonCanonUnderSel", "RandomUnderSel", 
                   "CanTotal", "NonTotal", "RandomTotal")

for (i in 1:length(soc)){
  can.tmp <- data.complete$Gene[data.complete$Class == "Canon"
                                                  & data.complete$SocOrigin == paste(soc[i])
                                                  & data.complete$adj_pvalue < 0.05]
  out$CanonUnderSel[i] <- length(can.tmp[!duplicated(can.tmp)])
  non.tmp <- data.complete$Gene[data.complete$Class == "NonCanon"
                                & data.complete$SocOrigin == paste(soc[i])
                                & data.complete$adj_pvalue < 0.05]
  out$NonCanonUnderSel[i] <- length(non.tmp[!duplicated(non.tmp)])
  ran.tmp <- data.complete$Gene[data.complete$Class == "Random"
                                & data.complete$SocOrigin == paste(soc[i])
                                & data.complete$adj_pvalue < 0.05]
  out$RandomUnderSel[i] <- length(can.tmp[!duplicated(ran.tmp)])
  can2.tmp <- data.complete$Gene[data.complete$SocOrigin == paste(soc)[i]
                                 & data.complete$Class == "Canon"]
  out$CanTotal[i] <- length(can2.tmp[!duplicated(can2.tmp)])
  non2.tmp <- data.complete$Gene[data.complete$SocOrigin == paste(soc)[i]
                                & data.complete$Class == "NonCanon"]
  out$NonTotal[i] <- length(non2.tmp[!duplicated(non2.tmp)])
  ran2.tmp <- data.complete$Gene[data.complete$Class == "Random"
                                 & data.complete$SocOrigin == paste(soc[i])]
  out$RandomTotal[i] <- length(ran2.tmp[!duplicated(ran2.tmp)])
}

out <- cbind(soc, out)
out

canprop <- vector(length = length(soc))
nonprop <- vector(length = length(soc))
ranprop <- vector(length = length(soc))
proptres <- vector(length = length(soc))

for (i in 1:length(soc)){
  one <- c(test[i,2], test[i,3], test[i,4])
  two <- c((test[i,5]-test[i,2]),(test[i,6]-test[i,3]),(test[i,7]-test[i,4]))
  three <- prop.test(one,two)
  canprop[i] <- as.numeric(formatC(three$estimate[[1]], digits = 3, format = "f"))
  nonprop[i] <- as.numeric(formatC(three$estimate[[2]], digits = 3, format = "f"))
  ranprop[i] <- as.numeric(formatC(three$estimate[[3]], digits = 3, format = "f"))
  proptres[i] <- three$p.value
}


out2 <- as.data.frame(cbind(paste(soc), canprop, nonprop, ranprop, proptres))
out2

out2$adjpvalue <- p.adjust(out2$proptres, method = "fdr")
out2

titles <- c("All Advanced", "All Social", "Apis", "Social Corbiculates",
            "Lasioglossum", "Melipona", "Ceratina", "All Solitary", "Habropoda", 
            "Megachile", "Novaeangliae")
out2 <- cbind(out2, titles)
out2
out2$V1 <- NULL
out2 <- out2[,c(6,1:5)]
out2
colnames(out2) <- c("Branch(es) Tested", "Canon", "Non-Canon", "Background", "pvalue", "adjpvalue")
out2

##by sociality


soc.craft <- function(x,y){
  one <- data.complete$Gene[data.complete$SocOrigin %in% paste(x)
                            & data.complete$Class == paste(y)
                            & data.complete$adj_pvalue < 0.05]
  one <- one[!duplicated(one)]
  two <- data.complete$Gene[data.complete$SocOrigin %in% paste(x)
                            & data.complete$Class == paste(y)]
  two <- two[!duplicated(two)]
  two <- two[!two %in% one]
  three <- c(length(one), length(two))
  return(three)
}

sol.can.prop <- soc.craft(solitary, "Canon")
soc.can.prop <- soc.craft(social, "Canon")
adv.can.prop <- soc.craft(complex, "Canon")
canon.prop <- rbind(sol.can.prop, soc.can.prop, adv.can.prop)
canon.prop

sol.non.prop <- soc.craft(solitary, "NonCanon")
soc.non.prop <- soc.craft(social, "NonCanon")
adv.non.prop <- soc.craft(complex, "NonCanon")
noncanon.prop <- rbind(sol.non.prop, soc.non.prop, adv.non.prop)

out3 <- cbind(canon.prop, noncanon.prop)
out3

sol.ran.prop <- soc.craft(solitary, "Random")
soc.ran.prop <- soc.craft(social, "Random")
adv.ran.prop <- soc.craft(complex, "Random")
random.prop <- rbind(sol.ran.prop, soc.ran.prop, adv.ran.prop)

out3 <- cbind(out3, random.prop)
out3

soc.canprop <- vector(length = nrow(out3))
soc.nonprop <- vector(length = nrow(out3))
soc.ranprop <- vector(length = nrow(out3))
soc.proptres <- vector(length = nrow(out3))

for (i in 1:nrow(out3)){
  one <- c(out3[i,1], out3[i,3], out3[i,5])
  two <- c(out3[i,2], out3[i,4], out3[i,6])
  three <- prop.test(one,two)
  soc.canprop[i] <- as.numeric(formatC(three$estimate[[1]], digits = 3, format = "f"))
  soc.nonprop[i] <- as.numeric(formatC(three$estimate[[2]], digits = 3, format = "f"))
  soc.ranprop[i] <- as.numeric(formatC(three$estimate[[3]], digits = 3, format = "f"))
  soc.proptres[i] <- three$p.value
}

Sociality <- c("Solitary", "Social", "Advanced Eusocial")

out4 <- as.data.frame(cbind(Sociality, soc.canprop, soc.nonprop, soc.ranprop, soc.proptres))
out4$adjpvalue <- p.adjust(out4$soc.proptres, method = "bonferroni")                      
colnames(out4) <- c("Sociality", "Canon", "Non-Canon", "Background", "pvalue", "adjpvalue")
out4

out5 <- out[,c(1:4)]
out5
out5$CanNoSelection <- out$CanTotal - out$CanonUnderSel
out5$NonNoSelection <- out$NonTotal - out$NonCanonUnderSel
out5$RanNoSelection <- out$RandomTotal - out$RandomUnderSel

#run chisq on can versus random and non versus random per lineage, 
#then correct for multiple testing

#canon
can.chi <- vector(length = length(soc))
can.chi.pvalue <- vector(length = length(soc))
for (i in 1:length(soc)){
  one <- c(out5[i,2], out5[i,5])
  two <- c(out5[i,4], out5[i,7])
  three <- rbind(one,two)
  four <- chisq.test(three, simulate.p.value = TRUE, B = 10000)
  can.chi[i] <- four$statistic
  can.chi.pvalue[i] <- four$p.value
}

#noncanon
non.chi <- vector(length = length(soc))
non.chi.pvalue <- vector(length = length(soc))
for (i in 1:length(soc)){
  one <- c(out5[i,3], out5[i,6])
  two <- c(out5[i,4], out5[i,7])
  three <- rbind(one,two)
  four <- chisq.test(three, simulate.p.value = TRUE, B = 10000)
  non.chi[i] <- four$statistic
  non.chi.pvalue[i] <- four$p.value
}


out6 <- as.data.frame(cbind(paste(soc), can.chi, can.chi.pvalue))
out6$can.chi.padj <- p.adjust(out6$can.chi.pvalue, method = "BH")
out6 <- cbind(out6, non.chi, non.chi.pvalue)
out6$non.chi.padj <- p.adjust(out6$non.chi.pvalue, method = "BH")
out6

four
out6.can <- out6[,c(1:4)]
out6.non <- out6[,c(1,5:7)]
out6.non

names(out6.can) <- c("Branch(es) Tested", "X-Squared", "p.value", "adj.p.value")
names(out6.non) <- c("Branch(es) Tested", "X-Squared", "p.value", "adj.p.value")

p.adjust(out6$can.chi, method = "fdr")

four$statistic


four
one <- as.data.frame(can.chi)
two <- as.data.frame(non.chi)
test.chi <- c(can.chi, non.chi)
test <- p.adjust(test.chi, method = "fdr")
test2
cbind(test.chi, test)

out
out2
out3
out4
out5
out6

out6.can$Canon.Prop.PSG <- out2$Canon
out6.can$Background.Prop.PSG <- out2$Background
out6.can$BY <- p.adjust(out6.can$p.value, method = "BY")


test.p <- as.data.frame(c(out6.can$p.value, out6.non$p.value))
test.p$BY <- p.adjust(test.p$`c(out6.can$p.value, out6.non$p.value)`, method = "BY")
test.p

out6.non$NonCanon.Prop.PSG <- out2$`Non-Canon`
out6.non$Background.Prop.PSG <- out2$Background

out6.non

out6.can$`X-Squared` <- as.numeric(formatC(out6.can$`X-Squared`, digits = 4, format = "f"))
out6.can$p.value <- as.numeric(formatC(out6.can$p.value, digits = 4, format = "f"))

for (i in 4:6){
  out6.can[,i] <- as.numeric(formatC(out6.can[,i], digits = 4, format = "f"))
}

out6.can
data.complete[data.complete$Gene == "Dl" & data.complete$SocOrigin == "Apis",]

goi.all[goi.all$count > 5,]
goi.all[goi.all$Gene == "Dl",]

imm.verse[imm.verse$GeneID == "LOC726167",]
imm.fun[imm.fun$Gene == "LOC726167",]

one <- c(0, 872)
two <- c(11, 4842)
three <- rbind(one,two)
four <- chisq.test(three, simulate.p.value = F)
four

out

data.complete[data.complete$SocOrigin == "CorbSoc"
              & data.complete$Class == "Canon"
              & data.complete$adj_pvalue < 0.05,]
data.complete[data.complete$Gene == "Tspan6",]


undersel.all[undersel.all$Gene == "Tspan6",]
