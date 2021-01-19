movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                   header = T, sep = ";")
head(movies)

require(ggplot2); require(plyr); require(gridExtra); require(grid);
library("plyr")
library("dplyr")
library("UpSetR")

between <- function(row, min, max){
  newData <- (row["ReleaseDate" < max] & (row["ReleaseDate"] > min))
}

plot1 <- function(mydata, x){
  myplot <- (ggplot(mydata, aes_string(x = x, fill = "color"))
             + geom_histogram() + scale_fill_identity()
             + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))
}

plot2 <- function(mydata, x, y){
  myplot <- (ggplot(data = mydata, aes_string(x = x, y = y, colour = "color"), alpha = 0.5)
             + geom_point() + scale_fill_identity()
             + theme_bw() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")))
}

attributeplots <- list(gridrows = 55,
                       plots = list(list(plot = plot1, x = "ReleaseDate", queries = F),
                                    list(plot = plot1, x = "ReleaseDate", queries = T),
                                    list(plot = plot2, x = "ReleaseDate", y = "AvgRating", queries = F),
                                    list(plot = plot2, x = "ReleaseDate", y = "AvgRating", queries = T)),
                       ncols = 3)


options(device = "quartz")
upset(movies, nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(T,F))
head(movies)



