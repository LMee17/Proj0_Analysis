#9th April 2021
#I have new reduced GO terms I need to read in and massage into a format that can
#allow me to make UpSetR plots or Chord Diagrams
#So here we go

library("dplyr")
library("stringr")
library("UpSetR")
library("ggplot2")

#Read in
filelist <- list.files(path = "output/ReViGo/", pattern = "*csv")

#read files into list
rvg.list <- lapply(filelist, function(x){
  read.table(file = paste("output/ReViGo/", x, sep = ""), header = T, sep = "\t", quote ="")
})


#goterms already exist in a previous object, go.base
#preparing a list space and counter for the following loop
#also adding a key (posh) to put "proper" labels in
socs <- vector()

for (i in 1:length(filelist)){
  socs[i] <- strsplit(filelist[i], "_", fixed = T)[[1]][1]
}
socs <- unique(socs)
posh <- c("All Advanced Eusocial", "All Social", "All Solitary", "Apis", "Ceratina",
          "Dufourea", "Habropoda", "Lasioglossum", "Megachile", "Melipona", 
          "Social Corbiculates")

#add notations from filenames to dataframes
for (i in 1:length(filelist)){
  rvg.list[[i]]$SocOrigin <- paste(strsplit(filelist[[i]], "_", fixed =T)[[1]][1])
  rvg.list[[i]]$Class <- paste(strsplit(filelist[[i]], "_", fixed = T)[[1]][2])
  rvg.list[[i]]$Ontology <- paste(strsplit(filelist[[i]], "_", fixed =T)[[1]][3])
}


#I can't combine until all the columns are the same types amongst dataframes
#as some have "Null" in plotX / plotY they're character in some places and numeric in others
for (i in 1:length(filelist)){
  rvg.list[[i]]$PlotX <- as.numeric(rvg.list[[i]]$PlotX)
  rvg.list[[i]]$PlotY <- as.numeric(rvg.list[[i]]$PlotY)
}
#Nulls will be replace with "NA" but who cares cos we're removing them anyway

#Add sociality, in case I need it

adv <- c("AllAdvanced", "Melipona", "Apis")
soc <- c("AllSoc", "Ceratina", "Lasioglossum", "SocCorb")
sol <- c("AllSol", "Dufourea", "Habropoda", "Megachile")

rvg.complete <- bind_rows(rvg.list)
tail(rvg.complete)

#remove whitespace from Eliminated column
rvg.complete$Eliminated <- trimws(rvg.complete$Eliminated)

rvg.complete$Sociality[rvg.complete$SocOrigin %in% adv] <- "Advanced_Eusocial"
rvg.complete$Sociality[rvg.complete$SocOrigin %in% soc] <- "Social"
rvg.complete$Sociality[rvg.complete$SocOrigin %in% sol] <- "Solitary"


for (i in 1:nrow(rvg.complete)){
  rvg.complete$Lineage[i] <- posh[(grep(paste(rvg.complete$SocOrigin[i]), socs))]
}


##Canon####
#Pull out all Canon terms order by lineage tested
rvg.can <- go.base

for (i in 1:length(posh)){
  tmp <- rvg.complete[rvg.complete$Eliminated == "False" & 
                        rvg.complete$Lineage == paste(posh[i]) &
                        rvg.complete$Class == "Canon",]
  tmp <- tmp$TermID
  rvg.can[,i+1] <- ifelse(rvg.can$GOterm %in% tmp, paste(1), paste(0))
  names(rvg.can)[i+1] <- paste(posh[i])
}
head(rvg.can)


for (i in 2:ncol(rvg.can)){
  rvg.can[,i] <- as.numeric(rvg.can[,i])
}

upset(rvg.can,
      sets = c("Apis", "Melipona", "All Advanced Eusocial",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Dufourea", "All Solitary"),
      keep.order = T)

#What if I added immune go terms into the mix ? 
go.immune <- read.table("Genome_Misc/ImmGOterms_GUS.txt", header = F)
head(go.immune)

#And also colour the set size by sociality ....
rvg.can$Imm <- ifelse(rvg.can$GOterm %in% go.immune$V1, paste(TRUE), paste(FALSE))
head(rvg.can[rvg.can$Imm == TRUE,])
rvg.can$Imm <- as.character(rvg.can$Imm)


upset(rvg.can, 
      query.legend = "top",
      mb.ratio = c(0.7, 0.3),
      nsets = length(socs),
      sets = c("Apis", "Melipona", "All Advanced Eusocial",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Dufourea", "All Solitary"),
      sets.bar.color = c("orange", "orange", "orange", "green", "green",
                         "green", "green", "blue", "blue", "blue", "blue"),
      order.by = c("freq"), keep.order = T,
      queries = list(list(query = elements,
                          params = list("Imm", "TRUE"), color = "red", active = T,
                          query.name = "Immune GOterms")))


#and add sociality as metdata item .. 

lineage <- rvg.complete$Lineage
sociality <- rvg.complete$Sociality

rvg.meta <- cbind(lineage, sociality)
rvg.meta <- as.data.frame(unique(rvg.meta))

#try and plot

upset(rvg.can, 
      query.legend = "top",
      mb.ratio = c(0.7, 0.3),
      nsets = length(socs),
      sets = c("Apis", "Melipona", "All Advanced Eusocial",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Dufourea", "All Solitary"),
      sets.bar.color = c("orange", "orange", "orange", "green", "green",
                         "green", "green", "blue", "blue", "blue", "blue"),
      sets.x.label = "Number GO terms from PSG",
      set.metadata = list(data = rvg.meta, plots = list(list(type = "matrix_rows", 
                                                             column = "sociality", assign = 10, 
                                                             colors = c(Advanced_Eusocial = "orange", 
                                                                        Social = "green",
                                                                        Solitary = "blue"))),
                          list(type = "text", column = "sociality", assign = 10, 
                               colors = c(Advanced_Eusocial = "orange", 
                                          Social = "green",
                                          Solitary = "blue"),
                               query.name = "Sociality")),
      order.by = c("freq"), keep.order = T,
      queries = list(list(query = elements,
                          params = list("Imm", "TRUE"), color = "violet", active = T,
                          query.name = "Immune GO terms")))


#make table so that I can have a stacked barchart with immune / not immune per sociality
rvg.can.stack <- rvg.complete[rvg.complete$Class == "Canon" &
                                rvg.complete$Eliminated == "False",]
rvg.can.stack <- rvg.can.stack[,c(1,15)]
head(rvg.can.stack)
rvg.can.stack$Imm <- ifelse(rvg.can.stack$TermID %in% go.immune$V1, paste(TRUE), paste(FALSE))

head(rvg.can.stack)
socs <- unique(rvg.can.stack$Sociality)

can.stack <- data.frame(matrix(nrow = 6, ncol = 3))
for (i in 1:length(socs)){
  can.stack[i,1] <- nrow(rvg.can.stack[rvg.can.stack$Sociality == paste(socs[i]) & rvg.can.stack$Imm == FALSE,])
  can.stack[i,2] <- paste(socs[i])
  can.stack[i,3] <- paste("Not immune")
  can.stack[i+3,1] <- nrow(rvg.can.stack[rvg.can.stack$Sociality == paste(socs[i]) & rvg.can.stack$Imm == TRUE,])
  can.stack[i+3,2] <- paste(socs[i])
  can.stack[i+3,3] <- paste("Immune")
}

names(can.stack) <- c("Number_GOterms", "Sociality", "Immune")
can.stack$Sociality[can.stack$Sociality=="Advanced_Eusocial"]<-"Advanced Eusocial"

can.stack <- as.data.frame(can.stack)
can.stack$Immune <- factor(can.stack$Immune, levels = c("Not immune", "Immune"))
can.stack[7,1] <- 71+51+74
can.stack[c(7,8),2] <- "All"
can.stack[7,3] <- "Not immune"
can.stack[8,1] <- 25+10+17
can.stack[8,3] <- "Immune"
can.stack

can.stack$Sociality <- factor(can.stack$Sociality, levels =
                                c("Solitary", "Social", "Advanced Eusocial", "All"))

ggplot(can.stack, aes(fill = Immune, y = Number_GOterms, x = Sociality)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_flip()

#add colour
paint3 <- c("red", "black")
names(paint3) <- levels(factor(c(levels(can.stack$Immune))))
paint3

ggplot(can.stack, aes(fill = Immune, y = Number_GOterms, x = Sociality)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_flip() +
  scale_fill_manual(name = "Immune", values = paint3)



##NonCanon#####
#Pull out all Noncanon terms order by lineage tested
rvg.non <- go.base

for (i in 1:length(posh)){
  tmp <- rvg.complete[rvg.complete$Eliminated == "False" & 
                        rvg.complete$Lineage == paste(posh[i]) &
                      rvg.complete$Class == "Non",]
  tmp <- tmp$TermID
  rvg.non[,i+1] <- ifelse(rvg.non$GOterm %in% tmp, paste(1), paste(0))
  names(rvg.non)[i+1] <- paste(posh[i])
}
head(rvg.non)

for (i in 2:ncol(rvg.non)){
  rvg.non[,i] <- as.numeric(rvg.non[,i])
}

upset(rvg.non,
      sets = c("Apis", "Melipona", "All Advanced Eusocial",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Dufourea", "All Solitary"),
      keep.order = T)

#What if I added immune go terms into the mix ? 
#And also colour the set size by sociality ....
rvg.non$Imm <- ifelse(rvg.non$GOterm %in% go.immune$V1, paste(TRUE), paste(FALSE))
head(rvg.non[rvg.non$Imm == TRUE,])
rvg.non$Imm <- as.character(rvg.non$Imm)


upset(rvg.non, 
      query.legend = "top",
      mb.ratio = c(0.7, 0.3),
      nsets = length(socs),
      sets = c("Apis", "Melipona", "All Advanced Eusocial",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Dufourea", "All Solitary"),
      sets.bar.color = c("orange", "orange", "orange", "green", "green",
          "green", "green", "blue", "blue", "blue", "blue"),
      order.by = c("freq"), keep.order = T,
      queries = list(list(query = elements,
                          params = list("Imm", "TRUE"), color = "red", active = T,
                          query.name = "Immune GOterms")))


#and add sociality as metdata item .. 

lineage <- rvg.complete$Lineage
sociality <- rvg.complete$Sociality

rvg.meta <- cbind(lineage, sociality)
rvg.meta <- as.data.frame(unique(rvg.meta))

#try and plot

upset(rvg.non, 
      query.legend = "top",
      mb.ratio = c(0.7, 0.3),
      nsets = length(socs),
      sets = c("Apis", "Melipona", "All Advanced Eusocial",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Dufourea", "All Solitary"),
      sets.bar.color = c("orange", "orange", "orange", "green", "green",
                         "green", "green", "blue", "blue", "blue", "blue"),
      sets.x.label = "Number GOterms from PSG",
      set.metadata = list(data = rvg.meta, plots = list(list(type = "matrix_rows", 
                              column = "sociality", assign = 10, 
                              colors = c(Advanced_Eusocial = "orange", 
                                         Social = "green",
                                         Solitary = "blue"))),
                          list(type = "text", column = "sociality", assign = 10, 
                                    colors = c(Advanced_Eusocial = "orange", 
                                               Social = "green",
                                               Solitary = "blue"),
                               query.name = "Sociality")),
      order.by = c("freq"), keep.order = T,
      queries = list(list(query = elements,
                          params = list("Imm", "TRUE"), color = "red", active = T,
                          query.name = "Immune GOterms")))

#make table so that I can have a stacked barchart with immune / not immune per sociality
rvg.non.stack <- rvg.complete[rvg.complete$Class == "Non" &
                                rvg.complete$Eliminated == "False",]
rvg.non.stack <- rvg.non.stack[,c(1,15)]
head(rvg.non.stack)
rvg.non.stack$Imm <- ifelse(rvg.non.stack$TermID %in% go.immune$V1, paste(TRUE), paste(FALSE))

head(rvg.non.stack)
socs <- unique(rvg.non.stack$Sociality)

non.stack <- data.frame(matrix(nrow = 6, ncol = 3))
for (i in 1:length(socs)){
  non.stack[i,1] <- nrow(rvg.non.stack[rvg.non.stack$Sociality == paste(socs[i]) & rvg.non.stack$Imm == FALSE,])
  non.stack[i,2] <- paste(socs[i])
  non.stack[i,3] <- paste("Not immune")
  non.stack[i+3,1] <- nrow(rvg.non.stack[rvg.non.stack$Sociality == paste(socs[i]) & rvg.non.stack$Imm == TRUE,])
  non.stack[i+3,2] <- paste(socs[i])
  non.stack[i+3,3] <- paste("Immune")
}

names(non.stack) <- c("Number_GOterms", "Sociality", "Immune")
non.stack$Sociality[non.stack$Sociality=="Advanced_Eusocial"]<-"Advanced Eusocial"

non.stack <- as.data.frame(non.stack)
non.stack$Immune <- factor(non.stack$Immune, levels = c("Not immune", "Immune"))
non.stack[7,1] <- sum(non.stack[c(1:3),1])
non.stack[c(7,8),2] <- "All"
non.stack[7,3] <- "Not immune"
non.stack[8,1] <- sum(non.stack[c(4:6),1])
non.stack[8,3] <- "Immune"
non.stack

non.stack$Sociality <- factor(non.stack$Sociality, levels =
                                c("Solitary", "Social", "Advanced Eusocial", "All"))

ggplot(non.stack, aes(fill = Immune, y = Number_GOterms, x = Sociality)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_flip()

#add colour
paint3 <- c("red", "black")
names(paint3) <- levels(factor(c(levels(can.stack$Immune))))
paint3

ggplot(non.stack, aes(fill = Immune, y = Number_GOterms, x = Sociality)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_flip() +
  scale_fill_manual(name = "Immune", values = paint3)


##Background####
rvg.back <- go.base

for (i in 1:length(posh)){
  tmp <- rvg.complete[rvg.complete$Eliminated == "False" & 
                        rvg.complete$Lineage == paste(posh[i]) &
                        rvg.complete$Class == "Background",]
  tmp <- tmp$TermID
  rvg.back[,i+1] <- ifelse(rvg.back$GOterm %in% tmp, paste(1), paste(0))
  names(rvg.back)[i+1] <- paste(posh[i])
}
head(rvg.back)

for (i in 2:ncol(rvg.back)){
  rvg.back[,i] <- as.numeric(rvg.back[,i])
}

upset(rvg.back,
      sets = c("Apis", "Melipona", "All Advanced Eusocial",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Dufourea", "All Solitary"),
      keep.order = T)


#add immune go terms
rvg.back$Imm <- ifelse(rvg.back$GOterm %in% go.immune$V1, paste(TRUE), paste(FALSE))
head(rvg.back[rvg.back$Imm == TRUE,])
rvg.back$Imm <- as.character(rvg.back$Imm)


upset(rvg.back, 
      query.legend = "top",
      mb.ratio = c(0.7, 0.3),
      nsets = length(socs),
      sets = c("Apis", "Melipona", "All Advanced Eusocial",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Dufourea", "All Solitary"),
      sets.bar.color = c("orange", "orange", "orange", "green", "green",
                         "green", "green", "blue", "blue", "blue", "blue"),
      order.by = c("freq"), keep.order = T,
      queries = list(list(query = elements,
                          params = list("Imm", "TRUE"), color = "red", active = T,
                          query.name = "Immune GOterms")))

upset(rvg.back, 
      query.legend = "top",
      mb.ratio = c(0.7, 0.3),
      nsets = length(socs),
      sets = c("Apis", "Melipona", "All Advanced Eusocial",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Dufourea", "All Solitary"),
      sets.bar.color = c("orange", "orange", "orange", "green", "green",
                         "green", "green", "blue", "blue", "blue", "blue"),
      sets.x.label = "Number GO terms from PSG",
      set.metadata = list(data = rvg.meta, plots = list(list(type = "matrix_rows", 
                                                             column = "sociality", assign = 10, 
                                                             colors = c(Advanced_Eusocial = "orange", 
                                                                        Social = "green",
                                                                        Solitary = "blue"))),
                          list(type = "text", column = "sociality", assign = 10, 
                               colors = c(Advanced_Eusocial = "orange", 
                                          Social = "green",
                                          Solitary = "blue"),
                               query.name = "Sociality")),
      order.by = c("freq"), keep.order = T,
      queries = list(list(query = elements,
                          params = list("Imm", "TRUE"), color = "red", active = T,
                          query.name = "Immune GO terms")))

#make table so that I can have a stacked barchart with immune / not immune per sociality
rvg.back.stack <- rvg.complete[rvg.complete$Class == "Background" &
                                rvg.complete$Eliminated == "False",]
rvg.back.stack <- rvg.back.stack[,c(1,15)]
head(rvg.back.stack)
rvg.back.stack$Imm <- ifelse(rvg.back.stack$TermID %in% go.immune$V1, paste(TRUE), paste(FALSE))

head(rvg.back.stack)
socs <- unique(rvg.back.stack$Sociality)

back.stack <- data.frame(matrix(nrow = 6, ncol = 3))
for (i in 1:length(socs)){
  back.stack[i,1] <- nrow(rvg.back.stack[rvg.back.stack$Sociality == paste(socs[i]) & rvg.back.stack$Imm == FALSE,])
  back.stack[i,2] <- paste(socs[i])
  back.stack[i,3] <- paste("Not immune")
  back.stack[i+3,1] <- nrow(rvg.back.stack[rvg.back.stack$Sociality == paste(socs[i]) & rvg.back.stack$Imm == TRUE,])
  back.stack[i+3,2] <- paste(socs[i])
  back.stack[i+3,3] <- paste("Immune")
}

names(back.stack) <- c("Number_GOterms", "Sociality", "Immune")
back.stack$Sociality[back.stack$Sociality== "Advanced_Eusocial"]<-"Advanced Eusocial"

back.stack <- as.data.frame(back.stack)
back.stack$Immune <- factor(back.stack$Immune, levels = c("Not immune", "Immune"))
back.stack[7,1] <- sum(back.stack[c(1:3),1])
back.stack[c(7,8),2] <- "All"
back.stack[7,3] <- "Not immune"
back.stack[8,1] <- sum(back.stack[c(4:6),1])
back.stack[8,3] <- "Immune"
back.stack

back.stack$Sociality <- factor(back.stack$Sociality, levels =
                                c("Solitary", "Social", "Advanced Eusocial", "All"))

ggplot(back.stack, aes(fill = Immune, y = Number_GOterms, x = Sociality)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_flip()

#add colour
paint3 <- c("red", "black")
names(paint3) <- levels(factor(c(levels(can.stack$Immune))))
paint3

ggplot(back.stack, aes(fill = Immune, y = Number_GOterms, x = Sociality)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_flip() +
  scale_fill_manual(name = "Immune", values = paint3)


##Combine#####
can.stack$Class <- paste("Canon")
non.stack$Class <- paste("Non-Canon")
back.stack$Class <- paste("Background")

all.stack <- rbind(can.stack, non.stack, back.stack)
all.stack

all.stack$Class <- factor(all.stack$Class, levels =
                                 c("Canon", "Non-Canon", "Background"))

ggplot(all.stack, aes(fill = Immune, y = Number_GOterms, x = Sociality)) +
  geom_bar(position = "fill", stat = "identity") +
  coord_flip() +
  facet_grid(~Class) +
  scale_fill_manual(name = "Immune", values = c("light grey", "red")) +
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(angle = 30)) +
  ylab("Proportion of GO terms considered functionally immune") +
  xlab(NULL)

# Get Proportions
all.stack

prop.imm <- data.frame(matrix(nrow = 12, ncol = 3))
socs <- unique(all.stack$Sociality)
class <- unique(all.stack$Class)

summary(all.stack$Class)

r <- 1
for (s in 1:length(socs)){
  df <- all.stack[all.stack$Sociality == paste(socs[s]),]
  for (c in 1:length(class)){
    t <- sum(df$Number_GOterms[df$Class == paste(class[c])])
    i <- df$Number_GOterms[df$Class == paste(class[c]) & df$Immune == "Immune"]
    p <- (i/t)*100
    prop.imm[r,] <- c(paste(socs[s]), paste(class[c]), paste(p))
    r <- r +1
    print(r)
  }
}

all.stack

prop.imm$X3 <- round(as.numeric(prop.imm$X3), digits = 2)
prop.imm

sig.test <- function(x,y){
  one <- all.stack[all.stack$Class == paste(x),]
  two <- all.stack[all.stack$Class == paste(y),]
}


25/96
all.stack
a <- c(25,17) 
b <- c(71,180)
prop.test(a,b)

##Upset ? 
head(rvg.can)
rvg.can$Class <- paste("Canon")
rvg.non$Class <- paste("Non-Canon")
rvg.back$Class <- paste("Background")
rvg.all <- rbind(rvg.can, rvg.non, rvg.back)
head(rvg.all)

upset(rvg.back,
      sets = c("Apis", "Melipona", "All Advanced Eusocial",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Dufourea", "All Solitary"),
      keep.order = T)

upset(rvg.all, 
      query.legend = "top",
      mb.ratio = c(0.7, 0.3),
      number.angles = 10, 
      #order = c(intersection size title, intersection size tick labels, 
      #set size title, set size tick labels, set names, numbers above bars)
      text.scale = c(1, .8, .8, 1, 1, .8),
      mainbar.y.label = "GO term Intersections",
      nsets = length(socs),
      sets = c("Apis", "Melipona", "All Advanced Eusocial",
               "Lasioglossum", "Ceratina", "Social Corbiculates", "All Social",
               "Megachile", "Habropoda", "Dufourea", "All Solitary"),
      sets.bar.color = c("orange", "orange", "orange", "green", "green",
                         "green", "green", "blue", "blue", "blue", "blue"),
      sets.x.label = "Number GO terms from PSG",
      set.metadata = list(data = rvg.meta, plots = list(list(type = "matrix_rows", 
                                                             column = "sociality", assign = 10, 
                                                             colors = c(Advanced_Eusocial = "orange", 
                                                                        Social = "green",
                                                                        Solitary = "blue"))),
                          list(type = "text", column = "sociality", assign = 10, 
                               colors = c(Advanced_Eusocial = "orange", 
                                          Social = "green",
                                          Solitary = "blue"),
                               query.name = "Sociality")),
      order.by = c("freq"), keep.order = T,
      queries = list(list(query = elements,
                          params = list("Class", "Canon"), color = "yellow", active = T,
                          query.name = "Canon Immune PSG GO terms"),
                     list(query = elements,
                          params = list("Class", "Background"), color = "sky blue", active = T,
                          query.name = "Non-Canon Immune PSG GO terms")))

#this doesn't work the way I wanted it to....
#idea .... do I condense the runs per sociality and then have a row per gene class ? 
#so there would be a social canon, social non-canon, social background, solitary canon etc
#this would then be overlaps between gene classes as well as socialities

rvg.comb <- go.base
key <- c("Advanced_Eusocial.Canon", "Advanced_Eusocial.Non", "Advanced_Eusocial.Background",
         "Social.Canon", "Social.Non", "Social.Background",
         "Solitary.Canon", "Solitary.Non", "Solitary.Background")

for (i in 1:length(key)){
  s <- str_split(key[i], "\\.")[[1]][1]
  c <- str_split(key[i], "\\.")[[1]][2]
  tmp <- rvg.complete[rvg.complete$Eliminated == "False" & 
                        rvg.complete$Sociality == paste(paste(s)) &
                        rvg.complete$Class == paste(c),]
  tmp <- tmp$TermID
  rvg.comb[,i+1] <- ifelse(rvg.comb$GOterm %in% tmp, paste(1), paste(0))
  names(rvg.comb)[i+1] <- paste(key[i])
}
head(rvg.comb)



names(rvg.comb) <- c("GOterms", "Advanced Eusocial Canon Immune", "Advanced Eusocial Non-Canon Immune", "Advanced Eusocial Background",
                     "Social Canon Immune", "Social Non-Canon Immune", "Social Background",
                     "Solitary Canon Immune", "Solitary Non-Canon Immune", "Solitary Background")
for (i in 2:ncol(rvg.comb)){
  rvg.comb[,i] <- as.numeric(rvg.comb[,i])
}

upset(rvg.comb,
      sets = c("Advanced Eusocial Canon Immune", "Advanced Eusocial Non-Canon Immune", "Advanced Eusocial Background",
               "Social Canon Immune", "Social Non-Canon Immune", "Social Background",
               "Solitary Canon Immune", "Solitary Non-Canon Immune", "Solitary Background"),
      keep.order = T)


#And also colour the set size by sociality ....
rvg.comb$Imm <- ifelse(rvg.comb$GOterm %in% go.immune$V1, paste(TRUE), paste(FALSE))
head(rvg.comb[rvg.comb$Imm == TRUE,])
rvg.comb$Imm <- as.character(rvg.comb$Imm)


upset(rvg.comb, 
      query.legend = "top",
      mb.ratio = c(0.5, 0.5),
      nsets = length(socs),
      sets = c("Advanced Eusocial Canon Immune", "Advanced Eusocial Non-Canon Immune", "Advanced Eusocial Background",
                       "Social Canon Immune", "Social Non-Canon Immune", "Social Background",
                       "Solitary Canon Immune", "Solitary Non-Canon Immune", "Solitary Background"),
      sets.bar.color = c("orange", "orange", "orange", "green", "green",
                         "green", "blue", "blue", "blue"),
      order.by = c("freq"), keep.order = T,
      queries = list(list(query = elements,
                          params = list("Imm", "TRUE"), color = "red", active = F,
                          query.name = "Immune GOterms")))


#and add sociality as metdata item .. 
key <- c("Advanced Eusocial Canon Immune", "Advanced Eusocial Non-Canon Immune", "Advanced Eusocial Background",
"Social Canon Immune", "Social Non-Canon Immune", "Social Background",
"Solitary Canon Immune", "Solitary Non-Canon Immune", "Solitary Background")
key2 <- c("Advanced_Eusocial", "Advanced_Eusocial", "Advanced_Eusocial", 
          "Social", "Social", "Social", "Solitary", "Solitary", "Solitary")
rvg.meta2 <- cbind(key,key2)

rvg.meta2 <- as.data.frame(rvg.meta2)
rvg.meta2

#try and plot

upset(rvg.comb, 
      query.legend = "top",
      mb.ratio = c(0.55, 0.45),
      nsets = length(socs),
      sets = c("Advanced Eusocial Background", "Advanced Eusocial Non-Canon Immune", "Advanced Eusocial Canon Immune",
               "Social Background", "Social Non-Canon Immune", "Social Canon Immune",
               "Solitary Background", "Solitary Non-Canon Immune", "Solitary Canon Immune"),
      sets.bar.color = c("orange", "orange", "orange", "green", "green",
                         "green", "blue", "blue", "blue"),
      sets.x.label = "Number GO terms from PSG",
      set.metadata = list(data = rvg.meta2, plots = list(list(type = "matrix_rows", 
                                                             column = "key2", assign = 10, 
                                                             colors = c(Advanced_Eusocial = "orange", 
                                                                        Social = "green",
                                                                        Solitary = "blue"))),
                          list(type = "text", column = "key2", assign = 10, 
                               colors = c(Advanced_Eusocial = "orange", 
                                          Social = "green",
                                          Solitary = "blue"),
                               query.name = "key2")),
      order.by = c("freq"), keep.order = T,
      queries = list(list(query = elements,
                          params = list("Imm", "TRUE"), color = "red", active = T,
                          query.name = "Immune GO terms")))


  str(rvg.meta)
str(rvg.meta2)



