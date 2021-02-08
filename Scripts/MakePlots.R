#######Making plots#####
######UpSetR#####

library("UpSetR")

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

####GO Terms#####
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
          ylab = "Omega", xlab = "Gene Class")
#omega_colourScatterplot.png
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
             cex = .6, fontface = "bold", fontfamily = "sans",
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer")

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
             lwd = 2, lty = "blank", fill = scheme,
             cex = .6, fontface = "bold", fontfamily = "sans",
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer")

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
             imagetype = "png",
             lwd = 2, lty = "blank", fill = scheme,
             cex = .6, fontface = "bold", fontfamily = "sans",
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer")

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
             imagetype = "png",
             lwd = 2, lty = "blank", fill = scheme,
             cex = .6, fontface = "bold", fontfamily = "sans",
             cat.cex = 0.6, cat.fontface = "bold", cat.default.pos = "outer")


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