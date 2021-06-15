#######Making plots#####
######UpSetR#####

library("UpSetR")
library("ggplot2")

####genes under selection####
upset(can.gus, mb.ratio = c(0.7, 0.3),
      sets.x.label = "Canon Gene Set Size",
      order.by = c("freq", "degree"), decreasing = c(T,F))
#Canon_GUS.png

upset(test.can, mb.ratio = c(0.7, 0.3),
      group.by = "sets",
      order.by = c("freq", "degree"), decreasing = c(T,F),
      boxplot.summary = "dN.dS")
#Can_GUS_omega_byset.png

upset(test.can, mb.ratio = c(0.7, 0.3),
      order.by = c("freq", "degree"), decreasing = c(T,F),
      boxplot.summary = "dN.dS")
#Can_GUS_omega.png

upset(non.gus, mb.ratio = c(0.7, 0.3),
      sets.x.label = "Canon Gene Set Size",
      order.by = c("freq", "degree"), decreasing = c(T,F))
#NonCanon_GUS.png

upset(test.non, mb.ratio = c(0.7, 0.3),
      group.by = "sets",
      order.by = c("freq", "degree"), decreasing = c(T,F),
      boxplot.summary = "dN.dS")
#Non_GUS_omega_byset.png

upset(test.non, mb.ratio = c(0.7, 0.3),
      order.by = c("freq", "degree"), decreasing = c(T,F),
      boxplot.summary = "dN.dS")
#Non_GUS_omega.png

upset(ran.gus, mb.ratio = c(0.7, 0.3),
      sets.x.label = "Canon Gene Set Size",
      order.by = c("freq", "degree"), decreasing = c(T,F))
#random_GUS.png

upset(test.ran, mb.ratio = c(0.7, 0.3),
      group.by = "sets",
      order.by = c("freq", "degree"), decreasing = c(T,F),
      boxplot.summary = "dN.dS")
#Ran_GUS_omega_byset.png

upset(test.ran, mb.ratio = c(0.7, 0.3),
      order.by = c("freq", "degree"), decreasing = c(T,F),
      boxplot.summary = "dN.dS")
#Ran_GUS_omega.png

####Making plots broken down by lineage (GUS)####
#start with canon
head(gene.undersel)
can.undersel <- gene.undersel[,c(1:2,5,11,14,20,23,26,29,32,35,38)]
can.undersel <- can.undersel[can.undersel$GeneID %in% can.base,]
head(can.undersel)
names(can.undersel) <- c("GeneID", "All Advanced", "All Social", "Apis", "Social Corbiculates", 
                         "Lasioglossum", "Melipona", "Ceratina", "All Solitary",
                         "Habropoda", "Megachile", "Novaeangliae")
for (i in 2:ncol(can.undersel)){
  can.undersel[,i] <- as.integer(can.undersel[,i])
}
upset(can.undersel, mb.ratio = c(0.7,0.3),
      sets = c("Apis", "Melipona", "All Advanced",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Novaeangliae", "All Solitary"),
      sets.x.label = "Genes Under Selection",
      order.by = c("freq", "degree"), keep.order = T)
#CanGUS_lineage_whole.png

upset(can.undersel, mb.ratio = c(0.7,0.3),
      sets = c("Apis", "Melipona",
               "Lasioglossum", "Ceratina", "Social Corbiculates",
               "Megachile", "Habropoda", "Novaeangliae"),
      sets.x.label = "Genes Under Selection",
      order.by = c("freq", "degree"), keep.order = T)
#CanGUS_lineage_minusALL.png

#with the problems removed
upset(can.undersel, mb.ratio = c(0.7,0.3),
      sets = c("Apis",
               "Ceratina", "Social Corbiculates",
               "Megachile", "Habropoda", "Novaeangliae"),
      sets.x.label = "Genes Under Selection",
      order.by = c("freq", "degree"), keep.order = T)
#CanGUS_lineage_minusOutlier.png

#NonCanon
head(gene.undersel)
non.undersel <- gene.undersel[,c(1,3,6,12,15,21,24,27,30,33,36,39)]
non.undersel <- non.undersel[non.undersel$GeneID %in% non.base,]
head(non.undersel)
names(non.undersel) <- c("GeneID", "All Advanced", "All Social", "Apis", "Social Corbiculates", 
                         "Lasioglossum", "Melipona", "Ceratina", "All Solitary",
                         "Habropoda", "Megachile", "Novaeangliae")
for (i in 2:ncol(non.undersel)){
  non.undersel[,i] <- as.integer(non.undersel[,i])
}
upset(non.undersel, mb.ratio = c(0.7,0.3),
      sets = c("Apis", "Melipona", "All Advanced",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Novaeangliae", "All Solitary"),
      sets.x.label = "Genes Under Selection",
      order.by = c("freq", "degree"), keep.order = T)
#NonGUS_lineage_whole.png

upset(non.undersel, mb.ratio = c(0.7,0.3),
      sets = c("Apis", "Melipona",
               "Lasioglossum", "Ceratina", "Social Corbiculates",
               "Megachile", "Habropoda", "Novaeangliae"),
      sets.x.label = "Genes Under Selection",
      order.by = c("freq", "degree"), keep.order = T)
#NonGUS_lineage_minusALL.png

#with the problems removed
upset(non.undersel, mb.ratio = c(0.7,0.3),
      sets = c("Apis",
               "Ceratina", "Social Corbiculates",
               "Megachile", "Habropoda", "Novaeangliae"),
      sets.x.label = "Genes Under Selection",
      order.by = c("freq", "degree"), keep.order = T)
#NonGUS_lineage_minusOutlier.png

#Random
head(gene.undersel)
ran.undersel <- gene.undersel[,c(1,4,7,13,16,22,25,28,31,34,37,40)]
ran.undersel <- ran.undersel[ran.undersel$GeneID %in% ran.base,]
head(ran.undersel)
names(ran.undersel) <- c("GeneID", "All Advanced", "All Social", "Apis", "Social Corbiculates", 
                         "Lasioglossum", "Melipona", "Ceratina", "All Solitary",
                         "Habropoda", "Megachile", "Novaeangliae")
for (i in 2:ncol(ran.undersel)){
  ran.undersel[,i] <- as.integer(ran.undersel[,i])
}
upset(ran.undersel, mb.ratio = c(0.7,0.3),
      sets = c("Apis", "Melipona", "All Advanced",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Novaeangliae", "All Solitary"),
      sets.x.label = "Genes Under Selection",
      order.by = c("freq", "degree"), keep.order = T)
#RanGUS_lineage_whole.png

upset(ran.undersel, mb.ratio = c(0.7,0.3),
      sets = c("Apis", "Melipona",
               "Lasioglossum", "Ceratina", "Social Corbiculates",
               "Megachile", "Habropoda", "Novaeangliae"),
      sets.x.label = "Genes Under Selection",
      order.by = c("freq", "degree"), keep.order = T)
#RanGUS_lineage_minusALL.png

#with the problems removed
upset(ran.undersel, mb.ratio = c(0.7,0.3),
      sets = c("Apis",
               "Ceratina", "Social Corbiculates",
               "Megachile", "Habropoda", "Novaeangliae"),
      sets.x.label = "Genes Under Selection",
      order.by = c("freq", "degree"), keep.order = T)
#RanGUS_lineage_minusOutlier.png

#Apidae versus Mega versus versus Halictidae
#Canon
upset(can.undersel, mb.ratio = c(0.7,0.3),
      sets = c("Apis", "Melipona", "Social Corbiculates", "Ceratina", "Habropoda",
               "Lasioglossum", "Novaeangliae",
               "Megachile"),
      sets.x.label = "Genes Under Selection",
      order.by = c("freq", "degree"), keep.order = T)
#CanGUS_orderbyPhylo.png

#NonCanon
upset(non.undersel, mb.ratio = c(0.7,0.3),
      sets = c("Apis", "Melipona", "Social Corbiculates", "Ceratina", "Habropoda",
               "Lasioglossum", "Novaeangliae",
               "Megachile"),
      sets.x.label = "Genes Under Selection",
      order.by = c("freq", "degree"), keep.order = T)
#NonGUS_orderbyPhylo.png

#Random
upset(ran.undersel, mb.ratio = c(0.7,0.3),
      sets = c("Apis", "Melipona", "Social Corbiculates", "Ceratina", "Habropoda",
               "Lasioglossum", "Novaeangliae",
               "Megachile"),
      sets.x.label = "Genes Under Selection",
      order.by = c("freq", "degree"), keep.order = T)
#RanGUS_orderbyPhylo.png

####GO Terms: UpSetR#####
##Noncanon (immune)
head(imm.terms)
imm.terms.work <- imm.terms[,c(1:3,5:8,10:11,13:14,16:17)]
head(imm.terms.work)
names(imm.terms.work) <- c("GOterm", "All Advanced", "All Social", "All Solitary", 
                           "Apis", "Canon Immune", "Social Corbiculates", "Habropoda",
                           "Lasioglossum", "Megachile", "Melipona",
                           "Novaeangliae", "Ceratina")

upset(imm.terms.work, mb.ratio = c(0.55, 0.45),
      sets = c("Canon Immune", "Apis", "Melipona",
               "Social Corbiculates", "Lasioglossum", "Ceratina",
               "Habropoda", "Megachile", "Novaeangliae"),
      order.by = c("freq", "degree"), keep.order = T)
#GOterms_imm_bylineage.png

upset(imm.terms.work, mb.ratio = c(0.55, 0.45),
      sets = c("Canon Immune", "Apis", "Melipona", "All Advanced",
               "Social Corbiculates", "Lasioglossum", "Ceratina", "All Social",
               "Habropoda", "Megachile", "Novaeangliae", "All Solitary"),
      order.by = c("freq", "degree"), keep.order = T)
#GOterms_imm_bylineage_all.png

head(soc.imm.terms)
soc.imm.terms.work <- soc.imm.terms[,c(1:2,5,6,7)]
head(soc.imm.terms.work)
names(soc.imm.terms.work) <- c("GOterm", "Canon Immune", "Solitary", "Social", 
                               "Advanced Social")

upset(soc.imm.terms.work, mb.ratio = c(0.7, 0.3),
      sets = c("Canon Immune", "Solitary", "Social",
               "Advanced Social"),
      order.by = c("freq", "degree"), keep.order = T)
#GOterms_imm_bysociality.png

head(ran.terms)
ran.terms.work <- ran.terms[,c(1:3,5:7,9:14)]
head(ran.terms.work)
names(ran.terms.work) <- c("GOterms", "All Advanced", "All Social", "All Solitary",
                           "Apis", "Social Corbiculates", "Habropoda", "Lasioglossum",
                           "Megachile", "Melipona", "Novaeangliae", "Ceratina")

upset(ran.terms.work, mb.ratio = c(0.7, 0.3),
      sets = c("Apis", "Melipona",
                 "Social Corbiculates", "Lasioglossum", "Ceratina",
                 "Habropoda", "Megachile", "Novaeangliae"),
      order.by = c("freq", "degree"), keep.order = T)
#GOterms_random_bylineage.png

upset(ran.terms.work, mb.ratio = c(0.7, 0.3),
      sets = c("Apis", "Melipona", "All Advanced",
               "Social Corbiculates", "Lasioglossum", "Ceratina", "All Social",
               "Habropoda", "Megachile", "Novaeangliae", "All Solitary"),
      order.by = c("freq", "degree"), keep.order = T)
#GOterms_random_bylineage_all.png

head(soc.ran.terms)
soc.ran.terms.work <- soc.ran.terms
head(soc.ran.terms.work)
names(soc.ran.terms.work) <- c("GOterms", "Solitary", "Social", "Advanced Social")

upset(soc.ran.terms.work, mb.ratio = c(0.7, 0.3),
      sets = c("Solitary", "Social", "Advanced Social"),
      order.by = c("freq", "degree"), keep.order = T)
#GOterms_random_bysociality.png

####Omega Scatterplot####
library("ggpubr")

ggline(omega, x = "Class", y = "dN.dS", 
       add = c("mean_se", "jitter"), 
       order = c("Canon", "NonCanon", "Random"),
       ylab = "Omega", xlab = "Gene Class")
#omega_plainScatterplot.png

ggboxplot(omega, x = "Class", y = "dN.dS", 
          color = "Class", palette = c("blue", "green", "orange"),
          order = c("Canon", "NonCanon", "Random"),
          panel.labs = c("Canon", "Non-Canon", "Background"),
          ylab = "Omega", xlab = "Gene Class")
#omega_colourScatterplot.png

omega2 <- omega
omega2 <- na.omit(omega2)
library("dplyr")
library("stringr")
levels(omega2$Class) <- gsub("Random", "Background", levels(omega2$Class))
levels(omega2$Class) <- gsub("NonCanon", "Non-Canon", levels(omega2$Class))
summary(omega2$Class)

ggboxplot(omega2, x = "Class", y = "dN.dS", 
          color = "Class", palette = c("blue", "green", "orange"),
          order = c("Canon", "Non-Canon", "Background"),
          ylab = "dNdS", xlab = "Gene Class")
#omega_colourScatterplot2.png



######Venn Diagram####

#Let's make a Venn diagram of overlapping genes under positive selection via phylogenetic
#family

#CANON

Apidae <- c("Apis", "Habro", "Meli", "CorbSoc", "Xylo")
Megachilidae <- "Mega"
Halictidae <- c("Halictid", "Nova")

api.can.gene <- data.complete$Gene[data.complete$SocOrigin %in% Apidae
                               & data.complete$Class == "Canon"
                               & data.complete$adj_pvalue < 0.05]
api.can.gene <- api.can.gene[!duplicated(api.can.gene)]

meg.can.gene <- data.complete$Gene[data.complete$SocOrigin %in% Megachilidae
                                   & data.complete$Class == "Canon"
                                   & data.complete$adj_pvalue < 0.05]
meg.can.gene <- meg.can.gene[!duplicated(meg.can.gene)]

hal.can.gene <- data.complete$Gene[data.complete$SocOrigin %in% Halictidae
                                   & data.complete$Class == "Canon"
                                   & data.complete$adj_pvalue < 0.05]
hal.can.gene <- hal.can.gene[!duplicated(hal.can.gene)]

install.packages("VennDiagram")
library("VennDiagram")

require("RColorBrewer")
scheme <- brewer.pal(3, "Dark2")

venn.diagram(x <- list(api.can.gene, meg.can.gene, hal.can.gene),
             category.names = c("Apidae", "Megachilidae", "Halictidae"),
             filename = "Plots/can_by_family_Venn.png", 
             output = T,
             imagetype = "png",
             lwd = 2, lty = "blank", fill = scheme,
             cex = 1, fontface = "bold", fontfamily = "sans",
             cat.cex = 1, cat.fontface = "bold", cat.default.pos = "outer")

#NONCANON

api.non.gene <- data.complete$Gene[data.complete$SocOrigin %in% Apidae
                                   & data.complete$Class == "NonCanon"
                                   & data.complete$adj_pvalue < 0.05]
api.non.gene <- api.non.gene[!duplicated(api.non.gene)]

meg.non.gene <- data.complete$Gene[data.complete$SocOrigin %in% Megachilidae
                                   & data.complete$Class == "NonCanon"
                                   & data.complete$adj_pvalue < 0.05]
meg.non.gene <- meg.non.gene[!duplicated(meg.non.gene)]

hal.non.gene <- data.complete$Gene[data.complete$SocOrigin %in% Halictidae
                                   & data.complete$Class == "NonCanon"
                                   & data.complete$adj_pvalue < 0.05]
hal.non.gene <- hal.non.gene[!duplicated(hal.non.gene)]

venn.diagram(x <- list(api.non.gene, meg.non.gene, hal.non.gene),
             category.names = c("Apidae", "Megachilidae", "Halictidae"),
             filename = "Plots/noncan_by_family_Venn.png", 
             output = T,
             imagetype = "png",
             cat.pos = c(-27, 27, 135),
             lwd = 2, lty = "blank", fill = scheme,
             cex = 1, fontface = "bold", fontfamily = "sans",
             cat.cex = 1, cat.fontface = "bold", cat.default.pos = "outer")

#RANDOM

api.ran.gene <- data.complete$Gene[data.complete$SocOrigin %in% Apidae
                                   & data.complete$Class == "Random"
                                   & data.complete$adj_pvalue < 0.05]
api.ran.gene <- api.ran.gene[!duplicated(api.ran.gene)]

meg.ran.gene <- data.complete$Gene[data.complete$SocOrigin %in% Megachilidae
                                   & data.complete$Class == "Random"
                                   & data.complete$adj_pvalue < 0.05]
meg.ran.gene <- meg.ran.gene[!duplicated(meg.ran.gene)]

hal.ran.gene <- data.complete$Gene[data.complete$SocOrigin %in% Halictidae
                                   & data.complete$Class == "Random"
                                   & data.complete$adj_pvalue < 0.05]
hal.ran.gene <- hal.ran.gene[!duplicated(hal.ran.gene)]


venn.diagram(x <- list(api.ran.gene, meg.ran.gene, hal.ran.gene),
             category.names = c("Apidae", "Megachilidae", "Halictidae"),
             filename = "Plots/random_by_family_Venn.png", 
             output = T,
             cat.pos = c(-27, 27, 135),
             imagetype = "png",
             lwd = 2, lty = "blank", fill = scheme,
             cex = 1, fontface = "bold", fontfamily = "sans",
             cat.cex = 1, cat.fontface = "bold", cat.default.pos = "outer")

#For shits and gigs, all of them together

api.all.gene <- data.complete$Gene[data.complete$SocOrigin %in% Apidae
                                   & data.complete$adj_pvalue < 0.05]
api.all.gene <- api.all.gene[!duplicated(api.all.gene)]

meg.all.gene <- data.complete$Gene[data.complete$SocOrigin %in% Megachilidae
                                   & data.complete$adj_pvalue < 0.05]
meg.all.gene <- meg.all.gene[!duplicated(meg.all.gene)]

hal.all.gene <- data.complete$Gene[data.complete$SocOrigin %in% Halictidae
                                   & data.complete$adj_pvalue < 0.05]
hal.all.gene <- hal.all.gene[!duplicated(hal.all.gene)]


venn.diagram(x <- list(api.all.gene, meg.all.gene, hal.all.gene),
             category.names = c("Apidae", "Megachilidae", "Halictidae"),
             filename = "Plots/all_by_family_Venn.png", 
             output = T,
             cat.pos = c(-27, 27, 135),
             imagetype = "png",
             lwd = 2, lty = "blank", fill = scheme,
             cex = 1, fontface = "bold", fontfamily = "sans",
             cat.cex = 1, cat.fontface = "bold", cat.default.pos = "outer")


#####Stacked BarCharts GUS versus lineage / sociality#####
soc <- data.complete$SocOrigin[!duplicated(data.complete$SocOrigin)]
soc <- soc[c(-3,-6)] #remove epi runs
soc

canon <- vector(length = length(soc))
noncanon <- vector(length = length(soc))
random <- vector(length = length(soc))

for (i in 1:length(soc)){
  one <- nrow(data.complete[data.complete$SocOrigin == paste(soc[i])
                            & data.complete$adj_pvalue < 0.05
                            & data.complete$Class == "Canon",])
  two <- nrow(data.complete[data.complete$SocOrigin == paste(soc[i])
                            & data.complete$adj_pvalue < 0.05
                            & data.complete$Class == "NonCanon",])
  three <- nrow(data.complete[data.complete$SocOrigin == paste(soc[i])
                              & data.complete$adj_pvalue < 0.05
                              & data.complete$Class == "Random",])
  canon[i] <- one
  noncanon[i] <- two
  random[i] <- three
}

chart.table <- rbind(canon, noncanon, random)
colnames(chart.table) <- paste(soc)
chart.table

barplot(chart.table, col = colors()[c(23, 89, 12)],
        border = "white", space = 0.04, xlab = "Lineage")

#ok so that's gross let's not use it.
#I'll do it again without the random genes.

chart.table <- rbind(canon, noncanon)
colnames(chart.table) <- paste(soc)
chart.table

barplot(chart.table, col = colors()[c(23, 89, 12)],
        border = "white", space = 0.04, xlab = "Lineage")

#Yeah ... this is just not worth doing 

#####Sociality versus number of genes under selection####

solitary <- c("AllSol", "Habro", "Mega", "Nova")
social <- c("AllOrigin", "Halictid", "Xylo", "CorbSoc")
complex <- c("AllComplex", "Meli", "Apis")

plot.data <- cbind(as.character(soc),canon,noncanon,random)
plot.data <- as.data.frame(plot.data)
colnames(plot.data)[1] <- "Lineage"
for (i in 2:ncol(plot.data)){
  plot.data[,i] <- as.integer(plot.data[,i])
}
plot.data$Sociality[plot.data$Lineage %in% solitary] <- "Solitary"
plot.data$Sociality[plot.data$Lineage %in% social] <- "Social"
plot.data$Sociality[plot.data$Lineage %in% complex] <- "Advanced Eusocial"
plot.data$Sociality2[plot.data$Lineage %in% solitary] <- 1
plot.data$Sociality2[plot.data$Lineage %in% social] <- 2
plot.data$Sociality2[plot.data$Lineage %in% complex] <- 3
plot.data$Family[plot.data$Lineage %in% Apidae] <- "Apidae"
plot.data$Family[plot.data$Lineage %in% Halictidae] <- "Halictidae"
plot.data$Family[plot.data$Lineage == "Mega"] <- "Megachilidae"

library(ggplot2)

ggplot(plot.data, aes(x = Sociality2, y = canon)) + 
  geom_point(aes(col=Sociality), size = 3) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_smooth(method="lm") +
  labs(title = "Canon Immune Genes vs Sociality", 
       y = "No Genes Under Selection", x = "Sociality")

#Hmm. Remove All runs I think.

plot.data2 <- na.omit(plot.data)

ggplot(plot.data2, aes(x = Sociality2, y = canon)) + 
  geom_point(aes(col=Sociality), size = 3) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_smooth(method="lm") +
  labs(title = "Canon Immune Genes vs Sociality", 
       y = "No Genes Under Selection", x = "Sociality")

prettylabels <- c("Apis", "Social Corbiculates", "Lasioglossum", "Melipona", "Ceratina",
                  "Habropoda", "Megachile", "Novaeangliae")
plot.data3 <- cbind(plot.data2, prettylabels)
colnames(plot.data3)[8] <- "DP"
plot.data3

require("ggrepel")

ggplot(plot.data3, aes(x = Sociality2, y = canon)) + 
  geom_point(aes(col=Sociality), size = 3) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_smooth(method="lm") +
  geom_text_repel(aes(label = DP), box.padding = .4) +
  labs(title = "Canon Immune Genes vs Sociality", 
       y = "No Genes Under Selection", x = "Sociality")
#CanonByLineage_lmtrendline_ggplot.png

plot.data3

ggplot(plot.data3, aes(x = Sociality2, y = canon)) + 
  geom_point(aes(col=Sociality), size = 3) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
#  geom_smooth(method="lm") +
  geom_text_repel(aes(label = DP), box.padding = .4) +
  labs(title = "Canon Immune Genes vs Sociality", 
       y = "No Genes Under Selection", x = "Sociality")
#CanonByLineage_ggplot.png

ggplot(plot.data3, aes(x = Sociality2, y = noncanon)) + 
  geom_point(aes(col=Sociality), size = 3) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_smooth(method="lm") +
  geom_text_repel(aes(label = DP), box.padding = .4) +
  labs(title = "Non-Canon Immune Genes vs Sociality", 
       y = "No Genes Under Selection", x = "Sociality")
#NonCanonByLineage_lmtrendline_ggplot.png

ggplot(plot.data3, aes(x = Sociality2, y = noncanon)) + 
  geom_point(aes(col=Sociality), size = 3) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  #  geom_smooth(method="lm") +
  geom_text_repel(aes(label = DP), box.padding = .4) +
  labs(title = "Non-Canon Immune Genes vs Sociality", 
       y = "No Genes Under Selection", x = "Sociality")
#NonCanonByLineage_ggplot.png

ggplot(plot.data3, aes(x = Sociality2, y = random)) + 
  geom_point(aes(col=Sociality), size = 3) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_smooth(method="lm") +
  geom_text_repel(aes(label = DP), box.padding = .4) +
  labs(title = "Background Genes vs Sociality", 
       y = "No Genes Under Selection", x = "Sociality")
#RandomByLineage_lmtrendline_ggplot.png

ggplot(plot.data3, aes(x = Sociality2, y = random)) + 
  geom_point(aes(col=Sociality), size = 3) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  #  geom_smooth(method="lm") +
  geom_text_repel(aes(label = DP), box.padding = .4) +
  labs(title = "Background Genes vs Sociality", 
       y = "No Genes Under Selection", x = "Sociality")
#RandomByLineage_ggplot.png

plot.data3$Total <- plot.data3$canon + plot.data3$noncanon + plot.data3$random
plot.data3

ggplot(plot.data3, aes(x = Sociality2, y = Total)) + 
  geom_point(aes(col=Sociality), size = 3) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_smooth(method="lm") +
  geom_text_repel(aes(label = DP), box.padding = .4) +
  labs(title = "All Genes vs Sociality", 
       y = "No Genes Under Selection", x = "Sociality")
#AllGenesByLineage_lmtrendline_ggplot.png

ggplot(plot.data3, aes(x = Sociality2, y = Total)) + 
  geom_point(aes(col=Sociality), size = 3) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  #  geom_smooth(method="lm") +
  geom_text_repel(aes(label = DP), box.padding = .4) +
  labs(title = "All Genes vs Sociality", 
       y = "No Genes Under Selection", x = "Sociality")
#AllGenesByLineage_ggplot.png

#Let's go back to this and add in the all lineages

plot.data

ggplot(plot.data, aes(x = Sociality2, y = canon)) + 
  geom_point(aes(col=Sociality), size = 3) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_smooth(method="lm") +
  labs(title = "Canon Immune Genes vs Sociality", 
       y = "No Genes Under Selection", x = "Sociality")

prettylabels2 <- c("All Advanced Eusocial", "All Social", "Apis", "Social Corbiculates", 
                   "Lasioglossum", "Melipona", "Ceratina", "All Solitary",
                  "Habropoda", "Megachile", "Novaeangliae")
plot.data <- cbind(plot.data, prettylabels2)
colnames(plot.data)[8] <- "DP"
plot.data

require("ggrepel")
require("ggpubr")

ggplot(plot.data, aes(x = Sociality2, y = canon)) + 
  geom_point(aes(col=Sociality), size = 3, alpha = 0.6 ) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_smooth(method="lm") +
  geom_text_repel(aes(label = DP), box.padding = .4) +
  labs(title = "Canon Immune Genes vs Sociality", 
       y = "No Genes Under Selection", x = "Sociality") +
  scale_color_manual(values = c("#e69f00", "#009e73", "#0072B2")) 
#CanonByLineage_lmtrendline_ggplot.png

#Hmmm. I see why I stopped this last time. May be just supplement with a bar chart ? 

plot.data7 <- plot.data[c(1:2,8),]
plot.data7

plot.data7$DP <- factor(plot.data7$DP, levels = c("All Solitary", "All Social", 
                                                  "All Advanced Eusocial"))

d <- ggplot(plot.data7, aes(x = DP, y = canon, fill = Sociality)) +
  geom_bar(stat = "identity", show.legend = F) +
  xlab("Sociality") +
  ylab("No Genes under Positive Selection") +
  #theme(axis.text.x = element_blank()) +
  #ggtitle("Canon Immune") +
  #theme(legend.position = c(.8, .8)) +
  scale_y_continuous(limits=c(0, 5)) +
  scale_fill_manual(values = c("#e69f00", "#009e73", "#0072B2")) 
#Allruns_canon_bylineage_ggplot.png


e <- ggplot(plot.data7, aes(x = DP, y = noncanon, fill = Sociality)) +
  geom_bar(stat = "identity", show.legend = F) +
  xlab("Sociality") +
  ylab("No Genes under Positive Selection") +
  #theme(axis.text.x = element_blank()) +
  #ggtitle("Non-Canon Immune") +
  scale_y_continuous(limits=c(0, 5)) +
  scale_fill_manual(values = c("#e69f00", "#009e73", "#0072B2"))
#Allruns_noncanon_bylineage_ggplot.png

f <- ggplot(plot.data7, aes(x = DP, y = random, fill = Sociality)) +
  geom_bar(stat = "identity", show.legend = F) +
  xlab("Sociality") +
  ylab("No Genes under Positive Selection") +
  #theme(axis.text.x = element_blank()) +
  #ggtitle("Background") +
  scale_y_continuous(limits=c(0, 50)) +
  scale_fill_manual(values = c("#e69f00", "#009e73", "#0072B2"))
#Allruns_background_bylineage_ggplot.png
  
a <- ggplot(plot.data3, aes(x = Sociality2, y = canon)) + 
  geom_point(aes(col=Sociality), size = 3) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(.15, .85)) +
  geom_smooth(method="lm", stat = 'smooth', colour = "dark gray", se = T) +
  geom_text_repel(aes(label = DP), box.padding = .4) +
  labs(y = "No Genes Under Selection", x = " ") +
  scale_color_manual(values = c("#e69f00", "#009e73", "#0072B2")) 

b <- ggplot(plot.data3, aes(x = Sociality2, y = noncanon)) + 
  geom_point(aes(col=Sociality), size = 3, show.legend = F) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_smooth(method="lm", stat = 'smooth', colour = "dark gray") +
  geom_text_repel(aes(label = DP), box.padding = .3) +
  labs(y = "No Genes Under Selection", x = " ") +
  scale_color_manual(values = c("#e69f00", "#009e73", "#0072B2")) 

c <- ggplot(plot.data3, aes(x = Sociality2, y = random)) + 
  geom_point(aes(col=Sociality), size = 3, show.legend = F) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_smooth(method="lm", stat = 'smooth', colour = "dark gray") +
  geom_text_repel(aes(label = DP), box.padding = .3) +
  labs(y = "No Genes Under Selection", x = " ") +
  scale_color_manual(values = c("#e69f00", "#009e73", "#0072B2")) 



fig <- ggarrange(a, b, c, d, e, f, 
                 labels = c("Canon Immune", "Non-Canon Immune", "Background", " ", " ", " "),
                 common.legend = T)
fig

#####GGPlot stacked bar chat for gus per class#######

#I need genes under selection per lineage

data.complete.sig <- data.complete[data.complete$adj_pvalue < 0.05,]

#aggregate by social lineage to get combinations of genes under selection
data.complete.sig <- data.complete.sig[,c(1,5,6)]
#remove epi runs
data.complete.sig <- data.complete.sig[!data.complete.sig$SocOrigin == "CorbSocEpi",]
data.complete.sig <- data.complete.sig[!data.complete.sig$SocOrigin == "AllOriginEpi",]

head(data.complete.sig)
unique(data.complete.sig$SocOrigin)

adv <- c("AllComplex", "Apis", "Meli")
soc <- c("AllOrigin", "CorbSoc", "Halictid", "Xylo")
sol <- c("AllSol", "Habro", "Mega", "Nova")

data.complete.sig$SocOrigin <- as.character(data.complete.sig$SocOrigin)

data.complete.sig$Sociality[as.character(data.complete.sig$SocOrigin) %in% adv] <- "Advanced Eusocial"
data.complete.sig$Sociality[as.character(data.complete.sig$SocOrigin) %in% soc] <- "Social"
data.complete.sig$Sociality[as.character(data.complete.sig$SocOrigin) %in% sol] <- "Solitary"

head(data.complete.sig)

#group by sociality
data.complete.sig$Sociality <- as.factor(data.complete.sig$Sociality)
tmp <- data.complete.sig[,c(1,3,4)]

tmp <- aggregate(tmp[3], tmp[-3], 
                               FUN = function(X) paste(unique(X), collapse=", "))
head(tmp)
combo <- unique(tmp$Sociality)
combo


#I'm gonna have to go through and hand code these combinations in terms of sociality (sob)
soc.code <- c("Social", "Advanced Eusocial and Solitary", 
              "Solitary", "Advanced Eusocial",
              "Advanced Eusocial and Solitary", "Advanced Eusocial and Social",
              "All Socialities", "Advanced Eusocial and Social",
              "Social and Solitary", "Social and Solitary",
              "All Socialities", "All Socialities",
              "All Socialities", "All Socialities")

#check
soc.key <- as.data.frame(cbind(combo, soc.code))
soc.key

for (i in 1:nrow(soc.key)){
  tmp$SocKey[tmp$Sociality == paste(soc.key[i,1])] <- paste(soc.key[i,2])
}

head(tmp)

#Canon
can.tmp <- tmp[tmp$Class == "Canon",]

plot.data4 <- data.frame()

for (i in 1:length(unique(soc.key$combo))){
  print(soc.key[i,2])
  plot.data4[i,1] <- paste(soc.key[i,2])
  one <- can.tmp[can.tmp$SocKey == paste(soc.key[i,2]),]
  plot.data4[i,2] <- nrow(one)
}

names(plot.data4) <- c("Sociality", "Number_PSGs")
plot.data4$Class <- paste("Canon")
plot.data4 <- unique(plot.data4)

library(ggplot2)

ggplot(plot.data4, aes(fill = Sociality, y = Class, x = Number_PSGs))+
  geom_bar(position = "stack", stat = "identity") +
  theme(legend.position = "bottom")


#I want these ideally as percentages so I can compare amongst socialities.
plot.data4$PercentagePSG <- (plot.data4$Number_PSGs / sum(plot.data4$Number_PSGs))*100
plot.data4

ggplot(plot.data4, aes(fill = Sociality, y = Class, x = PercentagePSG))+
  geom_bar(position = "stack", stat = "identity") +
  theme(legend.position = "bottom")

#add noncanon
non.tmp <- tmp[tmp$Class == "NonCanon",]
two <- data.frame()

for (i in 1:length(unique(soc.key$combo))){
  print(soc.key[i,2])
  two[i,1] <- paste(soc.key[i,2])
  one <- non.tmp[non.tmp$SocKey == paste(soc.key[i,2]),]
  two[i,2] <- nrow(one)
}
names(two) <- c("Sociality", "Number_PSGs")
two$Class <- paste("NonCanon")

two <- unique(two)
two$PercentagePSG <- (two$Number_PSGs / sum(two$Number_PSGs)*100)
two

plot.data5 <- rbind(plot.data4, two)
plot.data5

ggplot(plot.data5, aes(fill = Sociality, y = Class, x = PercentagePSG))+
  geom_bar(position = "stack", stat = "identity") +
  theme(legend.position = "bottom")

#add background

back.tmp <- tmp[tmp$Class == "Random",]
two <- data.frame()

for (i in 1:length(unique(soc.key$combo))){
  print(soc.key[i,2])
  two[i,1] <- paste(soc.key[i,2])
  one <- back.tmp[back.tmp$SocKey == paste(soc.key[i,2]),]
  two[i,2] <- nrow(one)
}
names(two) <- c("Sociality", "Number_PSGs")
two$Class <- paste("Background")
two <- unique(two)
two$PercentagePSG <- (two$Number_PSGs / sum(two$Number_PSGs)*100)
two

plot.data5 <- rbind(plot.data5, two)
plot.data5

ggplot(plot.data5, aes(fill = Sociality, y = Class, x = PercentagePSG))+
  geom_bar(position = "fill", stat = "identity") +
  theme(legend.position = "bottom")


#try and order the x and y axes
library("dplyr")
library("ggplot2")

plot.data5$Sociality <- factor(plot.data5$Sociality, levels=c("Advanced Eusocial", 
                                                "Advanced Eusocial and Social", "Social",
                                                "Advanced Eusocial and Solitary", "Social and Solitary", 
                                                "Solitary", "All Socialities"))
plot.data5$Class <- factor(plot.data5$Class, levels = c("Canon", "NonCanon", "Background"))
  
ggplot(plot.data5, aes(fill = Sociality, y = Number_PSGs, x = Class)) +
  geom_bar(position = "fill", stat = "identity") +
  theme(legend.position = "bottom") +
  xlab("Proportion of Genes Under Selection")
  
#I'm gonna try and add another bar per class showing social versus solitary versus all 
#I also don't need the percentage thing cos ggplot can do it itself with the position ="fill" option
plot.data5$Key <- paste(1)
plot.data6 <- plot.data5[,c(1:3,5)]
head(plot.data6)

three <- plot.data6
three$Key <- paste(2)
#recode the socialities (sigh)
social <- c("Social", "Advanced Eusocial", "Advanced Eusocial and Social")
three$Sociality[three$Sociality %in% social] <- "Social"
mix <- c("Social and Solitary", "Advanced Eusocial and Solitary")
three$Sociality[three$Sociality %in% mix] <- "Social and Solitary"
three

plot.data6 <- rbind(plot.data6, three)
plot.data6$Class<- as.character(plot.data6$Class)
plot.data6$Class[plot.data6$Class == "NonCanon"] <- "Non-Canon"
plot.data6$Class <- as.factor(plot.data6$Class)

#No idea if this is gonna work but let's try
gg1 <- ggplot(plot.data6, aes(fill = Sociality, y = Number_PSGs, x = Key)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~Class) +
  theme(legend.position = "right") +
  theme(axis.text.x = element_blank()) +
  ylab("Proportion of Genes Under Selection") +
  xlab("Gene Class")

#set colours  
unique(plot.data6$Sociality)
paint <- c("green", "purple" , "blue", "orange", "bluish green", "black", "yellow")
names(paint) <- levels(factor(c(levels(plot.data6$Sociality))))

paint2 <- c("#E69F00", "#F0E442", "#56B4E9", "#000000", "#009E73", "#D55E00", "#0072B2")
names(paint2) <- levels(factor(c(levels(plot.data6$Sociality))))

plot.data6$Class <- factor(plot.data6$Class, levels =
                                 c("Canon", "Non-Canon", "Background"))

ggplot(plot.data6, aes(fill = Sociality, y = Number_PSGs, x = Key)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~Class) +
  theme(legend.position = "right") +
  theme(axis.text.x = element_blank()) +
  ylab("Proportion of Genes Under Selection") +
  xlab("Gene Class") + 
  scale_fill_manual(name = "Sociality", values = paint2)


levels(plot.data6$Sociality)

plot.data6$Sociality <- factor(plot.data6$Sociality, levels=c("All Socialities", "Advanced Eusocial",
                                                              "Advanced Eusocial and Social", "Social",
                                                              "Advanced Eusocial and Solitary", "Social and Solitary",
                                                              "Solitary"))
plot.data5$Sociality <- factor(plot.data5$Sociality, levels=c("All Socialities", "Advanced Eusocial",
                                                              "Advanced Eusocial and Social", "Social",
                                                              "Advanced Eusocial and Solitary", "Social and Solitary",
                                                              "Solitary"))


ggplot(plot.data6, aes(fill = Sociality, y = Number_PSGs, x = Key)) +
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(~Class) +
  theme(legend.position = "bottom") +
  xlab("Proportion of Genes Under Selection")

#####GGPBoxplot for omega values#######
omega2 <- omega
omega2$Class2[omega2$Class == "Random"] <- "Background"
omega2$Class2[omega2$Class == "Canon"] <- "Canon"
omega2$Class2[omega2$Class == "NonCanon"] <- "Non-Canon"


head(omega)
ggplot(omega2, aes(x = Class2, y = dN.dS, fill = Class2)) +
  geom_violin() +
  geom_boxplot(width = .1)

ggplot(omega, aes(x = Class, y = dN.dS)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .01)


ggplot(omega2, aes(x = Class2, y = dN.dS, fill = Class2)) +
  geom_violin() +
  geom_boxplot(width = .1) +
  coord_flip() +
  xlab("Gene Class") +
  guides(fill = guide_legend(title = "Gene Class")) +
  scale_x_discrete(limits = c("Background", "Non-Canon", "Canon")) +
  ylab("dN/dS ratio") +
#  scale_fill_brewer(palette = "PuRd")
 scale_fill_manual(values = c("#56B4E9", "#F0E442", "#CC79A7")) 
  
3500*12
  