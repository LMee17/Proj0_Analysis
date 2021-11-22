#######Making plots#####
######UpSetR#####

##V2: nov 2021##

#install.packages("UpSetR")
#install.packages("VennDiagram")

library("UpSetR")
library("VennDiagram")
library("ggplot2")
require("RColorBrewer")
library("ggrepel")

######Venn Diagram####

#Let's make a Venn diagram of overlapping genes under positive selection via phylogenetic
#family

#CANON

Apidae <- c("Apis", "Habro", "Meli", "CorbSoc", "Ceratina")
Megachilidae <- "Mega"
Halictidae <- c("Lasio", "Dufourea")


api.can.gene <- data.df$Gene[data.df$SocOrigin %in% Apidae
                                   & data.df$Class == "Canon Immune"
                                   & data.df$adj_pvalue < 0.05]
api.can.gene <- api.can.gene[!duplicated(api.can.gene)]

meg.can.gene <- data.df$Gene[data.df$SocOrigin %in% Megachilidae
                                   & data.df$Class == "Canon Immune"
                                   & data.df$adj_pvalue < 0.05]
meg.can.gene <- meg.can.gene[!duplicated(meg.can.gene)]

hal.can.gene <- data.df$Gene[data.df$SocOrigin %in% Halictidae
                                   & data.df$Class == "Canon Immune"
                                   & data.df$adj_pvalue < 0.05]
hal.can.gene <- hal.can.gene[!duplicated(hal.can.gene)]

scheme <- brewer.pal(3, "Dark2")

venn.diagram(x <- list(api.can.gene, meg.can.gene, hal.can.gene),
             category.names = c("Apidae", "Megachilidae", "Halictidae"),
             filename = "Plots/Revisions/can_by_family_Venn.png", 
             output = T,
             imagetype = "png",
             lwd = 2, lty = "blank", fill = scheme,
             cex = 3, fontface = "bold", fontfamily = "sans",
             cat.cex = 2, cat.fontface = "bold", cat.default.pos = "outer", cat.dist = 0.15)

#NONCANON


api.non.gene <- data.df$Gene[data.df$SocOrigin %in% Apidae
                                   & data.df$Class == "Non-Canon"
                                   & data.df$adj_pvalue < 0.05]
api.non.gene <- api.non.gene[!duplicated(api.non.gene)]

meg.non.gene <- data.df$Gene[data.df$SocOrigin %in% Megachilidae
                                   & data.df$Class == "Non-Canon"
                                   & data.df$adj_pvalue < 0.05]
meg.non.gene <- meg.non.gene[!duplicated(meg.non.gene)]

hal.non.gene <- data.df$Gene[data.df$SocOrigin %in% Halictidae
                                   & data.df$Class == "Non-Canon"
                                   & data.df$adj_pvalue < 0.05]
hal.non.gene <- hal.non.gene[!duplicated(hal.non.gene)]

venn.diagram(x <- list(api.non.gene, meg.non.gene, hal.non.gene),
             category.names = c("Apidae", "Megachilidae", "Halictidae"),
             filename = "Plots/Revisions/noncan_by_family_Venn.png", 
             output = T,
             imagetype = "png",
             cat.pos = c(-27, 27, 135),
             lwd = 2, lty = "blank", fill = scheme,
             cex = 3, fontface = "bold", fontfamily = "sans",
             cat.cex = 2, cat.fontface = "bold", cat.default.pos = "outer")

#RANDOM

api.ran.gene <- data.df$Gene[data.df$SocOrigin %in% Apidae
                                   & data.df$Class == "Background"
                                   & data.df$adj_pvalue < 0.05]
api.ran.gene <- api.ran.gene[!duplicated(api.ran.gene)]

meg.ran.gene <- data.df$Gene[data.df$SocOrigin %in% Megachilidae
                                   & data.df$Class == "Background"
                                   & data.df$adj_pvalue < 0.05]
meg.ran.gene <- meg.ran.gene[!duplicated(meg.ran.gene)]

hal.ran.gene <- data.df$Gene[data.df$SocOrigin %in% Halictidae
                                   & data.df$Class == "Background"
                                   & data.df$adj_pvalue < 0.05]
hal.ran.gene <- hal.ran.gene[!duplicated(hal.ran.gene)]


venn.diagram(x <- list(api.ran.gene, meg.ran.gene, hal.ran.gene),
             category.names = c("Apidae", "Megachilidae", "Halictidae"),
             filename = "Plots/Revisions/background_by_family_Venn.png", 
             output = T,
             cat.pos = c(-27, 27, 135),
             imagetype = "png",
             lwd = 2, lty = "blank", fill = scheme,
             cex = 3, fontface = "bold", fontfamily = "sans",
             cat.cex = 2, cat.fontface = "bold", cat.default.pos = "outer")

#For shits and gigs, all of them together

api.all.gene <- data.df$Gene[data.df$SocOrigin %in% Apidae
                                   & data.df$adj_pvalue < 0.05]
api.all.gene <- api.all.gene[!duplicated(api.all.gene)]

meg.all.gene <- data.df$Gene[data.df$SocOrigin %in% Megachilidae
                                   & data.df$adj_pvalue < 0.05]
meg.all.gene <- meg.all.gene[!duplicated(meg.all.gene)]

hal.all.gene <- data.df$Gene[data.df$SocOrigin %in% Halictidae
                                   & data.df$adj_pvalue < 0.05]
hal.all.gene <- hal.all.gene[!duplicated(hal.all.gene)]


venn.diagram(x <- list(api.all.gene, meg.all.gene, hal.all.gene),
             category.names = c("Apidae", "Megachilidae", "Halictidae"),
             filename = "Plots/Revisions/all_by_family_Venn.png", 
             output = T,
             cat.pos = c(-27, 27, 135),
             imagetype = "png",
             lwd = 2, lty = "blank", fill = scheme,
             cex = 3, fontface = "bold", fontfamily = "sans",
             cat.cex = 2, cat.fontface = "bold", cat.default.pos = "outer")





######dNdS######

#violin plot of evolutionary rate
plot.data.1 <- data.df
plot.data.1$Function <- factor(plot.data.1$Function, levels=c("Receptor", 
                               "Signalling", "Effector",
                               "Non-Canon Immune", "Background"))

plot.data.1$Class <- factor(plot.data.1$Class, levels=c("Canon Immune", 
                                                     "Non-Canon", "Background"))


p <- ggplot(data = plot.data.1, aes(x = Function, y = dN.dS, fill = Class)) +
  geom_boxplot(outlier.shape = 16, outlier.size = .8, alpha = .75) +
  scale_fill_manual(values = c("#cccc00",
                               "#800080",
                               "#1e90ff"),
                    name = "Gene Class") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  xlab("Gene Function") +
  ylab("dN/dS Ratio")

png("Plots/Revisions/Omega_BoxPlot1.png")
p
dev.off()

#######Number of Genes Under Selection######
#Previously was simply bar charts .... with scattergraphs ? 

#scattergraph

head(plot.data.1)
plot.data.1$SocKey[plot.data.1$BranchTested == "Elaboration of Sociality"] <- 3
plot.data.1$SocKey[plot.data.1$BranchTested == "Origin of Sociality"] <- 2
plot.data.1$SocKey[plot.data.1$BranchTested == "Solitary"] <- 1

#I only want to count those that are under selection atm, and I don't want duplicate genes
pd2 <- plot.data.1[plot.data.1$PSG == 1,]

pd2$BranchTested <- factor(pd2$BranchTested, levels = c("Solitary", 
                                                        "Origin of Sociality", 
                                                        "Elaboration of Sociality"))

ggplot(data = pd2) +
  geom_bar(aes(x = BranchTested, fill = BranchTested)) + 
  facet_wrap(~ Class, scales = "free") +
  scale_fill_manual(values = c("#00578a",
                              "#00cc00",
                              "#ffa500")) +
  theme(legend.title = element_blank()) + 
  scale_x_discrete(labels = NULL, breaks = NULL) + 
  labs(x = "", y = "Number of Genes Under Positive Selection")
  #scale_x_discrete(labels = function(x) str_wrap(x, width = 7))
  
#Plots/Revisions/PSGCount_BarChart1.png

ggplot(data = pd2) +
  geom_bar(aes(x = SocOrigin, fill = Class)) 

one <- pd2[,c(5,17:18)]
one <- one[!duplicated(one),]
rownames(one) <- c(1,2,3,4,5,6,7,8,9,10,11)

countPSG <- function(x,y){
  psg <- nrow(pd2[pd2$SocOrigin == paste(x) & pd2$Class == paste(y),])
  return(psg)
}

for (i in 1:nrow(one)){
  one[i,4] <- countPSG(one[i,1], "Background")
  one[i,5] <- countPSG(one[i,1], "Canon Immune")
  one[i,6] <- countPSG(one[i,1], "Non-Canon")
}

names(one)[4:6] <- c("Background", "Canon", "NonCanon")
pd3 <- one
pd3$DN <- c("All Elaboration", "All Origin", "All Solitary", "Apis", "Ceratina", 
            "Corbiculates", "Dufourea", "Habropoda", 
            "Lasioglossum", "Megachile", "Melipona")

#I don't want this plot any more

#But maybe a piechart of the immune functions under selection ? Per branch test type
#No, there's not enough canon PSG to even bother

pd3 <- data.df[data.df$Class == "Canon Immune" & data.df$PSG == 1,]
pd3 <- pd3[,c(1,6)]
pd3 <- pd3[!duplicated(pd3),]
pd3

tot <- nrow(pd3)
tot.sig <- nrow(pd3[pd3$Function == "Signalling",])
tot.rec <- nrow(pd3[pd3$Function =="Receptor",])
tot.eff <- nrow(pd3[pd3$Function == "Effector",])
tots <- rbind(tot, tot.sig, tot.rec, tot.eff)
prop <- as.numeric(tots[2:4])/42
fnc <- unique(pd3$Function)
pd4 <- as.data.frame(cbind(prop, fnc))
pd4$prop <- as.numeric(pd4$prop)
pd4$count <- tots[2:4]

pd4$hsize <- 4

ggplot(pd4, aes(x = hsize, y = prop, fill = fnc)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(0.2, hsize + 0.5)) +
  geom_text(aes(label = count),
          position = position_stack(vjust = 0.5)) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_fill_manual(values = c("#af6eeb",
                               "#00cc00",
                               "#51e5ff"),
                    name = "Immune Function") 

#Plots/Revisions/PSGImmFunct_DonutChart1.png

######GO Upset Plot #####
#This is gonnna be a pain

#Start wiht readin gin go analysis results

filelist <- list.files(pattern = "*tsv", path = "output/TopGo/")

go.res.list <- lapply(filelist, function(x){
  read.table(file = paste("output/TopGo/", x, sep = ""), header = T, sep = "\t", quote ="")
})
for (i in 1:length(filelist)){
  go.res.list[[i]]$SocOrigin <- paste(strsplit(filelist[[i]], "_", fixed =T)[[1]][1])
  go.res.list[[i]]$Class <- paste(strsplit(filelist[[i]], "_", fixed = T)[[1]][2])
  go.res.list[[i]]$Ontology <- paste(strsplit(filelist[[i]], "_", fixed =T)[[1]][3])
}

socs <- vector() 
for (i in 1:length(filelist)){
  socs[i] <- strsplit(filelist[i], "_", fixed = T)[[1]][1]
}
socs <- unique(socs)
socs <- socs[-c(1,3)]
DN <- c("All Elaboration", "All Origin", "All Solitary", "Apis", "Ceratina", 
        "Corbiculates", "Dufourea", "Habropoda", 
        "Lasioglossum", "Megachile", "Melipona")

#remove immune onlly stuff
go.res.list <- go.res.list[-c(1,3)]

go.res.df$BranchTest[go.res.df$SocOrigin %in% elaboration] <- "Elaboration of Sociality"
go.res.df$BranchTest[go.res.df$SocOrigin %in% origin] <- "Origin of Sociality"
go.res.df$BranchTest[go.res.df$SocOrigin %in% solitary] <- "Solitary"

#make a basic framework of my future frames with all go terms
go.base <- read.table("Genome_Misc/AllGoTerms_Nov21.txt", header = F)
names(go.base) <- "GoTerm"

sol.can <- unique(go.res.df$GO.ID[go.res.df$BranchTest == "Solitary" &
                       go.res.df$Class == "Canon Immune"])
sol.non <- unique(go.res.df$GO.ID[go.res.df$BranchTest == "Solitary" &
                                    go.res.df$Class == "Non-Canon"])
sol.back <- unique(go.res.df$GO.ID[go.res.df$BranchTest == "Solitary" &
                                    go.res.df$Class == "Background"])

ori.can <- unique(go.res.df$GO.ID[go.res.df$BranchTest == "Origin of Sociality" &
                                    go.res.df$Class == "Canon Immune"])
ori.non <- unique(go.res.df$GO.ID[go.res.df$BranchTest == "Origin of Sociality" &
                                    go.res.df$Class == "Non-Canon"])
ori.back <- unique(go.res.df$GO.ID[go.res.df$BranchTest == "Origin of Sociality" &
                                    go.res.df$Class == "Background"])

ela.can <- unique(go.res.df$GO.ID[go.res.df$BranchTest == "Elaboration of Sociality" &
                                    go.res.df$Class == "Canon Immune"])
ela.non <- unique(go.res.df$GO.ID[go.res.df$BranchTest == "Elaboration of Sociality" &
                                    go.res.df$Class == "Non-Canon"])
ela.back <- unique(go.res.df$GO.ID[go.res.df$BranchTest == "Elaboration of Sociality" &
                                    go.res.df$Class == "Background"])

go.base$Can_Sol <- ifelse(go.base$GoTerm %in% sol.can, 1, 0)
go.base$Can_Ori <- ifelse(go.base$GoTerm %in% ori.can, 1, 0)
go.base$Can_Ela <- ifelse(go.base$GoTerm %in% ela.can, 1, 0)
go.base$Non_Sol <- ifelse(go.base$GoTerm %in% sol.non, 1, 0)
go.base$Non_Ori <- ifelse(go.base$GoTerm %in% ori.non, 1, 0)
go.base$Non_Ela <- ifelse(go.base$GoTerm %in% ela.non, 1, 0)
go.base$Back_Sol <- ifelse(go.base$GoTerm %in% sol.back, 1, 0)
go.base$Back_Ori <- ifelse(go.base$GoTerm %in% ori.back, 1, 0)
go.base$Back_Ela <- ifelse(go.base$GoTerm %in% ela.back, 1, 0)

names(go.base) <- c("GoTerms", "Solitary Canon Immune", "Origin Canon Immune", "Elaboration Canon Immune",
          "Solitary Non-Canon Immune", "Origin Non-Canon Immune", "Elaboration Non-Canon Immune",
                        "Solitary Background", "Origin Background", "Elaboration Background")

#add immune go terms
imm.go <- read.table("ImmResources/ImmGoTerms_Nov21.tsv", header = F, sep = "\t")
pd5 <- go.base
pd5$Immune <- ifelse(pd5$GoTerm %in% imm.go$V1, paste(TRUE), paste(FALSE))

upset(pd5,
      query.legend = "top",
      mb.ratio = c(0.55, 0.45),
      nsets = 6,
      sets = c("Elaboration Background", "Elaboration Non-Canon Immune", "Elaboration Canon Immune",
               "Origin Background", "Origin Non-Canon Immune", "Origin Canon Immune",
               "Solitary Background", "Solitary Non-Canon Immune", "Solitary Canon Immune"),
      order.by = c("freq"), keep.order = T,
      queries = list(list(query = elements,
                          params = list("Immune", "TRUE"), color = "red", active = T,
                          query.name = "Immune GO terms")))

meta.pd5 <- as.data.frame(names(pd5[2:10]))
meta.pd5$key <- c("Sol", "Ori", "Ela", "Sol", "Ori", "Ela", "Sol", "Ori", "Ela")
names(meta.pd5$`names(pd5[2:10])`) <- "setName"

upset(pd5,
      query.legend = "top",
      mb.ratio = c(0.55, 0.45),
      number.angles = 350,
      line.size = .8,
      point.size = 3,
      nsets = 9,
      sets = c("Elaboration Background", "Elaboration Non-Canon Immune", "Elaboration Canon Immune",
               "Origin Background", "Origin Non-Canon Immune", "Origin Canon Immune",
               "Solitary Background", "Solitary Non-Canon Immune", "Solitary Canon Immune"),
      sets.bar.color = c("orange", "orange", "orange", "green", "green",
                         "green", "blue", "blue", "blue"),
      order.by = c("freq"), keep.order = T,
      sets.x.label = "Number GO terms from PSG",
      set.metadata = list(data = meta.pd5, 
                          plots = list(list(type = "matrix_rows", 
                                              column = "key", assign = 10, 
                                             colors = c(Ela = "orange", 
                                                        Ori = "green",
                                                        Sol = "blue")))),
      queries = list(list(query = elements,
                          params = list("Immune", "TRUE"), color = "red", active = F,
                          query.name = "Immune GO terms")))

#Plots/Revisions/GOTerms_Upset1.png

######GC Content vs dNdS######

pd6 <- data.df
pd6$Imm <- ifelse(pd6$Class == "Background", paste(FALSE), paste(TRUE))
pd6.back <- pd6[pd6$Class == "Background",]
pd6.imm <- pd6[!pd6$Class == "Background",]
pd6.imm[3562,] <- pd6.back[1,]
pd6.non <- pd6[pd6$Class == "Non-Canon",]

pd6.imm$Function <- factor(pd6.imm$Function, levels = c("Receptor",
                                                        "Signalling", 
                                                        "Effector", 
                                                        "Non-Canon Immune", 
                                                        "Background"))

ggplot(pd6.back, aes(x = log(dN.dS), y = log(GeneAverage))) +
  geom_point(alpha = 0.01, size = 0.5) +
  geom_point(data = pd6.imm, aes(x = log(dN.dS), y = log(GeneAverage),
                                 color = Function, shape = Class),
                                  size = 1.5) + 
  stat_ellipse(data = pd6.non, colour = "red") +
  scale_color_manual(values = c("turquoise",
                                "#9acd32",
                                "#1e90ff",
                                "red",
                                "gray"), 
                     name = "Gene Function") +
  geom_line(stat="smooth",method = "lm",
            size = .5,
            linetype ="dashed",
            alpha = 0.5) + 
  xlab("log dN/dS ratio") +
  ylab("log GC content")

  
#Plots/Revisions/GCvsOmega_Scatter1.png
work4$Function <- factor(work4$Function, levels = c("Receptor",
                                                        "Signalling", 
                                                        "Effector", 
                                                        "Non-Canon Immune", 
                                                        "Background"))


ggplot(work4, aes(x = logOmega, y = logGC)) +
  geom_point(aes(color = Function, alpha = Immune, shape = Class)) +
  scale_alpha(guide = "none") +
  stat_ellipse(data = work4[work4$Function == "Non-Canon Immune",], colour = "violet") +
  scale_color_manual(values = c("#009E73",
                                "#F0E442",
                                "#56B4E9",
                                "violet",
                                "gray"), 
                     name = "Gene Function") +
  scale_shape(name = "Gene Class") +
  geom_line(stat="smooth",method = "lm",
            size = .5,
            linetype ="dashed",
            alpha = 0.5) + 
  xlab("log dN/dS ratio") +
  ylab("log GC content")

#Plots/Revisions/GCvsOmega_Scatter2.png


######GCvsSelection#####

pd7 <- work5
head(pd7)

ggplot(data = pd7, aes(x = Selection, y = AvgBranch, fill = Selection))+
  geom_boxplot(alpha = .75) +
  facet_wrap(~ BranchTested) +
  scale_fill_manual(values = c("#9370db", "#7fffd4"),
    guide = "none") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  xlab("") +
  ylab("GC Content (%)")

#Plots/Revisions/GCvsSel_Boxplot1.png

#######GC per Species#####
gc3$Sociality[gc3$Sociality == "Simple (facultative)"] <- "Simple (Facultative) Social"
gc3$Sociality <- factor(gc3$Sociality, levels = c("All",
                                                  "Solitary",
                                                  "Simple (Facultative) Social",
                                                  "Simple Eusocial",
                                                  "Complex Eusocial"))
pd8 <- gc3
pd8$Species <- str_replace(pd8$Species, "Apis mellifera", "A.mellifera")
pd8$Species <- str_replace(pd8$Species, "Aflo", "A.florea")
pd8$Species <- str_replace(pd8$Species, "M.quadrifasciata", "Mel.quadrifasciata")
pd8$Species <- str_replace(pd8$Species, "Bter", "B.terrestris")
pd8$Species <- str_replace(pd8$Species, "Bimp", "B.impatiens")
pd8$Species <- str_replace(pd8$Species, "Lalb", "L.albipes")
pd8$Species <- str_replace(pd8$Species, "Emex", "E.mexicana")
pd8$Species <- str_replace(pd8$Species, "Ccal", "C.calcarata")
pd8$Species <- str_replace(pd8$Species, "Mrot", "Meg.rotundata")
pd8$Species <- str_replace(pd8$Species, "Hlab", "H.laboriosa")
pd8$Species <- str_replace(pd8$Species, "Dnov", "D.novaeangliae")

pd8$Species <- factor(pd8$Species, levels = c("All Species",
                                              "D.novaeangliae",
                                              "H.laboriosa", 
                                              "Meg.rotundata",
                                              "C.calcarata",
                                              "E.mexicana",
                                              "L.albipes",
                                              "B.impatiens",
                                              "B.terrestris",
                                              "Mel.quadrifasciata",
                                              "A.florea",
                                              "A.mellifera"))

ggplot(data = pd8, aes(x = Species, y = GC, fill = Sociality)) +
  geom_violin(alpha = .8, scale = "width", trim = F) +
  geom_boxplot(width = 0.2) +
  coord_flip() +
  ylab("GC Content (%)")

#Plots/Revisions/GCbySpecies_ViolinPlot1.png
