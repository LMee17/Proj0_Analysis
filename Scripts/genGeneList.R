#29th September 2020
#Generate gene peptide ID lists for GOterm analysis

#Non-Canonical immune genes (perhaps with those that do not have orthologs across Insecta ? Depends on whether I can find the raw materials in order to make a list)
#non conserved insecta genes under Hym_NonCan.txt
#All NC gene products
head(data.all)
nci <- data.all[,c(1,6)]
nci <- nci[!duplicated(nci),]
nci <- nci[nci$Class == "NonCanon",]
nci.2 <- merge(nci, reduced, by = "Gene")
nci.2<-nci.2[,3]

write.table(nci.2, "input/GOI/NonCan_input.txt", col.names = F, row.names = F, quote = F)

#All non-canonical immune genes under selection
head(undersel)
ncus <- merge(undersel, reduced, by = "Gene")
head(ncus)
ncus <- ncus[ncus$Class == "NonCanon",]
head(ncus)
ncus.2 <- as.data.frame(ncus[,6])
head(ncus.2)

write.table(ncus.2, "input/GOI/NonCan_Undersel_all.txt", col.names = F, row.names = F, quote = F)

#do this per origin / elaboration of sociality
origin <- c("AllOrigin", "Xylo", "Halictid", "CorbSoc")
ori <- data.all[data.all$adj_pvalue < 0.05 & data.all$Class == "NonCanon",]
ori <- ori[ori$SocOrigin %in% origin,]
ori.2 <- merge(ori, reduced, by = "Gene")
ori.2 <- as.data.frame(ori.2[,10])

head(ori.2)

write.table(ori.2, "input/GOI/NonCan_Undersel_origin.txt", col.names = F, row.names = F, quote = F)

complex <- c("Apis", "Meli", "AllComplex")
com <- data.all[data.all$adj_pvalue < 0.05 & data.all$Class == "NonCanon",]
com <- com[com$SocOrigin %in% complex,]
com.2 <- merge(com, reduced, by = "Gene")
com.2 <- as.data.frame(com.2[,10])
head(com.2)

write.table(com.2, "input/GOI/NonCan_Undersel_complex.txt", col.names = F, row.names = F, quote = F)

episodic <- c("CorbSocEpi", "AllOriginEpi")
epi <- data.all[data.all$adj_pvalue < 0.05 & data.all$Class == "NonCanon",]
epi <- epi[epi$SocOrigin %in% episodic,]
head(epi)
epi.2 <- merge(epi, reduced, by = "Gene")
epi.2 <- as.data.frame(epi.2[,10])
head(epi.2)

write.table(epi.2, "input/GOI/NonCan_UnderSel_epi.txt", col.names = F, row.names = F, quote = F)

#Canonical immune genes under selection (though I suspect this won't 
#be all that interesting)

can <- undersel[undersel$Class == "Canon",]
can <- merge(can, reduced, by = "Gene")
head(can)
can.2 <- as.data.frame(can[,6])

write.table(can.2, "input/GOI/Can_UnderSel.txt", col.names = F, row.names = F, quote = F)                                   

#All random genes under selection (sanity check to compare to other papers that 
#assessed selection on this genome set)
origin <- c("AllOrigin", "Xylo", "Halictid", "CorbSoc")
ori <- data.all[data.all$adj_pvalue < 0.05 & data.all$Class == "Random",]
ori <- ori[ori$SocOrigin %in% origin,]
ori.2 <- merge(ori, reduced, by = "Gene")
ori.2 <- as.data.frame(ori.2[,10])

head(ori.2)

write.table(ori.2, "input/GOI/Random_Undersel_origin.txt", col.names = F, row.names = F, quote = F)

complex <- c("Apis", "Meli", "AllComplex")
com <- data.all[data.all$adj_pvalue < 0.05 & data.all$Class == "Random",]
com <- com[com$SocOrigin %in% complex,]
com.2 <- merge(com, reduced, by = "Gene")
com.2 <- as.data.frame(com.2[,10])
head(com.2)

write.table(com.2, "input/GOI/Random_Undersel_complex.txt", col.names = F, row.names = F, quote = F)

episodic <- c("CorbSocEpi", "AllOriginEpi")
epi <- data.all[data.all$adj_pvalue < 0.05 & data.all$Class == "Random",]
epi <- epi[epi$SocOrigin %in% episodic,]
head(epi)
epi.2 <- merge(epi, reduced, by = "Gene")
epi.2 <- as.data.frame(epi.2[,10])
head(epi.2)

write.table(epi.2, "input/GOI/Random_UnderSel_epi.txt", col.names = F, row.names = F, quote = F)

