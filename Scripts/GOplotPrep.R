#9th April 2021
#I have new reduced GO terms I need to read in and massage into a format that can
#allow me to make UpSetR plots or Chord Diagrams
#So here we go

library("dplyr")
library("stringr")
library("UpSetR")

#Read in
filelist <- list.files(path = "output/ReViGo/", pattern = "*csv")

#read files into list
rvg.list <- lapply(filelist, function(x){
  read.table(file = paste("output/ReViGo/", x, sep = ""), header = T, sep = "\t", quote ="")
})

head(rvg.list[[1]])


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


rvg.complete$Sociality[rvg.complete$SocOrigin %in% adv] <- "Advanced_Eusocial"
rvg.complete$Sociality[rvg.complete$SocOrigin %in% soc] <- "Social"
rvg.complete$Sociality[rvg.complete$SocOrigin %in% sol] <- "Solitary"


for (i in 1:nrow(rvg.complete)){
  rvg.complete$Lineage[i] <- posh[(grep(paste(rvg.complete$SocOrigin[i]), socs))]
}

tail(rvg.complete)


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


#Background
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
      keep.order = T,
      queries = list(list(query = elements,
                          params = list("Imm", "TRUE"), color = "red", active = T,
                          query.name = "Immune GOterms")))

