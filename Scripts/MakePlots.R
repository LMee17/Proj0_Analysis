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