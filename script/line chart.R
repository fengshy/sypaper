rm(list = ls())
gc()

library(abind)
library(ggplot2)
dirname <- dir("G:/1.dyland/output data/historical/area/", full.names = TRUE)

carea <- list()
for(i in 1:length(dirname)){
  load(dirname[i])
  carea[[i]] <- area[109:165,]
}

Area <- abind(carea, along = 2)
Area <- data.frame(Area)
Area$year <- 1958:2014




