rm(list = ls())
gc()
library(ncdf4)
library(stringr)
library(plyr)
library(sp)
library(raster)
library(abind)
library(magrittr)
library(maptools)
library(fields)
setwd("H:/CMIP6_data/")


#-------------------------------------------------------------------------------
#1.detect the models
Files = list()
Files[[1]] <- dir("hfls/historical/", full.names = TRUE)
Files[[2]] <- dir("hfss/historical/", full.names = TRUE)
Files[[3]] <- dir("hurs/historical/", full.names = TRUE)
Files[[4]] <- dir("tasmax/historical/", full.names = TRUE)
Files[[5]] <- dir("tasmin/historical/", full.names = TRUE)
Files[[6]] <- dir("sfcWind/historical/", full.names = TRUE)
Files[[7]] <- dir("pr/historical/", full.names = TRUE)

filess <- list()
for(i in 1:length(Files)){
  nums <- which(str_detect(Files[[i]], paste('_r1i1p1f1_', sep = '')))
  filess[[i]] <- Files[[i]][nums]
}
save(filess, file = "G:/1.dyland/output data/historical/other/filess_historical.RData")

{
  name1 <- unique(str_replace_all(filess[[1]], '.+Amon_|_historical.+', ""))
  name2 <- unique(str_replace_all(filess[[2]], '.+Amon_|_historical.+', ""))
  name3 <- unique(str_replace_all(filess[[3]], '.+Amon_|_historical.+', ""))
  name4 <- unique(str_replace_all(filess[[4]], '.+Amon_|_historical.+', ""))
  name5 <- unique(str_replace_all(filess[[5]], '.+Amon_|_historical.+', ""))
  name6 <- unique(str_replace_all(filess[[6]], '.+Amon_|_historical.+', ""))
  name7 <- unique(str_replace_all(filess[[7]], '.+Amon_|_historical.+', ""))
}
model <- intersect(intersect(intersect(intersect(name1,name2),intersect(name3,name4)),intersect(name5, name6)),name7)
save(model, file = "G:/1.dyland/output data/historical/other/model_historical.RData")
#---------------------------------------------------------------------delete up

load("G:/1.dyland/output data/historical/other/filess_historical.RData")
load("G:/1.dyland/output data/historical/other/model_historical.RData")
model <- model[-10]
#2.read the variable

historical <- function(x){
  llply(1:length(model), function(i){
    Data <- bind_GCM(files = filess[[x]], name = model[i])
    variable <- Data$variable
    return(variable)
  }, .progress = 'text')
}


variable <- historical(x = 1)


#3. cal pet of historical
{
  load("G:/1.dyland/output data/historical/variable/hfls.RData")
  load("G:/1.dyland/output data/historical/variable/hfss.RData")
  load("G:/1.dyland/output data/historical/variable/hurs.RData")
  load("G:/1.dyland/output data/historical/variable/sfcWind.RData")
  load("G:/1.dyland/output data/historical/variable/tasmax.RData")
  load("G:/1.dyland/output data/historical/variable/tasmin.RData")
}
load("G:/1.dyland/output data/historical/other/dem_historical.RData")
pet <- function(i){
  PET <- llply(1:dim(hfls[[i]])[3],function(x){
    Penman_future(Tmax = tasmax[[i]][,,x] - 273.16, Tmin = tasmin[[i]][,,x] - 273.16,
                  RH = hurs[[i]][,,x], WS = sfcWind[[i]][,,x],
                  SH = hfss[[i]][,,x]*24*60*60/1e6, LH = hfls[[i]][,,x]*24*60*60/1e6, Z = dem[[i]] )
  }, .progress = 'text')
}

for(i in 1:length(model)){
  Pet <- pet(i)
  Pet <- abind(Pet, along = 3)
  attr(Pet, 'dimnames') <- NULL
  path = paste("G:/1.dyland/output data/historical/pet/",model[i],"_historical.RData")
  save(Pet, file = path)
  rm(Pet)
  print(i)
}



#4.read in nc

dirname <- dir(path = "G:/1.dyland/output data/historical/pet/", full.names = TRUE)
for(i in 1:length(dirname)){
  load(dirname[i])
  lonlat <- bind_lonlat(files = filess[[1]], name = model[i])
  lon <- lonlat[[1]]$lon
  lat <- lonlat[[1]]$lat
  time <- bind_time(files = filess[[1]], name = model[i])
  nums <- which(str_detect(filess[[1]], paste(model[i],sep = '')))
  nc <- nc_open(filess[[1]][nums[1]])
  lon <- ncdim_def(name = 'lon', "degrees_east", vals = lon)
  lat <- ncdim_def(name = 'lat', 'degrees_north', vals = lat)
  time <- ncdim_def(name = 'time',calendar = nc$dim$time$calendar,units = nc$dim$time$units,vals = time)
  pet <- ncvar_def(name = 'pet',units = 'mm day-1', dim = list(lon, lat, time),missval = 1e+20, longname = 'potential evaporation',prec = 'float')
  ncnew <- nc_create(filename = paste("G:/1.dyland/output data/historical/pet_nc/","pet_Amon_",model[i],"_historical_r1i1p1f1.nc",sep = ''),vars = list(pet))
  ncvar_put(nc = ncnew, varid = pet,vals = Pet)
  print(i)
}

#5. cal AI of historical
load("G:/1.dyland/output data/historical/variable/pr.RData")
dirnames <- dir("G:/1.dyland/output data/historical/pet_nc/", full.names = TRUE)

for(i in 1:length(dirnames)){
  nc <- nc_open(dirnames[i])
  pet <- ncvar_get(nc, varid = 'pet')
  cal <- rep(1:(dim(pet)[3]/12), each = 12)
  petmm <- apply(pet, c(1,2), function(x){
    aggregate(x = x, by = list(cal), mean, na.rm = TRUE)$x
  })
  prmm <- apply(pr[[i]],c(1,2), function(x){
    aggregate(x = x, by = list(cal), mean, na.rm = TRUE)$x
  })
  AI <- prmm*24*3600/petmm
  path <- paste("G:/1.dyland/output data/historical/AI/", model[i],"_historical.RData")
  save(AI,file = path)
  print(i)
}


#6. cal area of historical

dirname <- dir("G:/1.dyland/output data/historical/AI/", full.names = TRUE)
dirnames <- dir("G:/1.dyland/output data/historical/pet_nc/", full.names = TRUE)

for(i in 1:length(model)){
  load(dirname[i])
  AI <- aperm(AI, c(2,3,1))
  nc <- nc_open(dirnames[i])
  lon = nc$dim$lon$vals
  lat = nc$dim$lat$vals
  lat = seq(min(lat),max(lat),length.out = length(lat))
  lonlat <- expand.grid(lon,lat)
  Area <- F_area(lonlat)
  area <- matrix(NA, dim(AI)[3])
  for(j in 1:dim(AI)[3]){
    area[j] <- sum(Area[AI[,,j] < 0.65], na.rm = TRUE)
  }
  path <- paste("G:/1.dyland/output data/historical/area/",model[i],"area_historical.RData")
  save(area, file = path)
  print(i)
}




















