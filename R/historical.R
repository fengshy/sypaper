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


#1.Interp 10 model of hfls hfss hurs tasmax tasmin sfcWind tas
setwd("H:/CMIP6_data/")
Files = list()
Files[[1]] <- dir("hfls/historical/", full.names = TRUE)
Files[[2]] <- dir("hfss/historical/", full.names = TRUE)
Files[[3]] <- dir("hurs/historical/", full.names = TRUE)
Files[[4]] <- dir("tasmax/historical/", full.names = TRUE)
Files[[5]] <- dir("tasmin/historical/", full.names = TRUE)
Files[[6]] <- dir("sfcWind/historical/", full.names = TRUE)
Files[[7]] <- dir("tas/historical/", full.names = TRUE)
Files[[8]] <- dir("pr/historical/", full.names = TRUE)

{
  name1 <- unique(str_replace_all(Files[[1]], '.+Amon_|_historical.+', ""))
  name2 <- unique(str_replace_all(Files[[2]], '.+Amon_|_historical.+', ""))
  name3 <- unique(str_replace_all(Files[[3]], '.+Amon_|_historical.+', ""))
  name4 <- unique(str_replace_all(Files[[4]], '.+Amon_|_historical.+', ""))
  name5 <- unique(str_replace_all(Files[[5]], '.+Amon_|_historical.+', ""))
  name6 <- unique(str_replace_all(Files[[6]], '.+Amon_|_historical.+', ""))
}
model <- intersect(intersect(intersect(name1,name2),intersect(name3,name4)),intersect(name5, name6))
rm(name1)
rm(name2)
rm(name3)
rm(name4)
rm(name5)
rm(name6)
model <- model[-6][-10][-11][-12][-15]
model <- model[-10]
model <- model[-9]
save(model, file = "G:/1.dyland/output data/CMIP6 interp/historical-model.RData")


nc3 <- function(file, varid){
  var_value <- llply(varid, function(x){
    ncvar_get(nc = nc_open(file), x)
  })
  names(var_value) <- varid
  return(var_value)
}

bind_GCMs <- function(files, name){
  variable <- str_sub(string = files[1], start = 1, end = 2)
  nums <- which(str_detect(files, paste(name, '_historical_r1i1p1f1_', sep = '')))
  data <- llply(nums, function(i){
    data <- nc3(file = files[i], varid = c('lon', 'lat', 'time', variable))
    data
  })
  var_value <- llply(data, function(x) x[[4]]) %>% abind(., along = 3)
  attr(var_value, 'dimnames') <- NULL
  return(list(lon = data[[1]]$lon, lat = data[[1]]$lat, date = data[[1]]$time, variable = var_value))
  rm(data); gc()
}


historical <- function(x){
  llply(1:10, function(i){
    Data <- bind_GCMs(files = Files[[x]], name = model[i])
    lon <- Data$lon
    lat <- Data$lat
    variable <- Data$variable
    loc <- c(which(lon > 180), which(lon < 180))
    lon <- lon[loc]
    lon[lon > 180] <- lon[lon > 180] - 360
    variable <- variable[loc,,]

    newdata <- list()
    for(j in 1:dim(variable)[3]){
      newdata[[j]] <- interp.surface.grid(obj = list(x = lon, y = lat, z = variable[,,i]),
                                          grid.list = list(x = seq(-179.75,179.75,0.5),
                                                           y = seq(-89.75,89.75,0.5)))$z
    }
    variable <- abind(newdata, along = 3)
    attr(variable, 'dimnames') <- NULL
    return(variable)
  }, .progress = 'text')
}

hfls <- historical(x = 1)
hfss <- historical(x = 2)
hurs <- historical(x = 3)
tasmax <- historical(x = 4)
tasmin <- historical(x = 5)
sfcWind <- historical(x = 6)
tas <- historical(x = 7)
pr <- historical(x = 8)

save(hfls, file = "G:/1.dyland/output data/CMIP6 interp/hfls.RData")
save(hfss, file = "G:/1.dyland/output data/CMIP6 interp/hfss.RData")
save(hurs, file = "G:/1.dyland/output data/CMIP6 interp/hurs.RData")
save(sfcWind, file = "G:/1.dyland/output data/CMIP6 interp/sfcWind.RData")
save(tasmax, file = "G:/1.dyland/output data/CMIP6 interp/tasmax.RData")
save(tasmin, file = "G:/1.dyland/output data/CMIP6 interp/tasmin.RData")
save(tas, file = "G:/1.dyland/output data/CMIP6 interp/tas.RData")
save(pr, file = "G:/1.dyland/output data/CMIP6 interp/pr.RData")



#2. cal PET of 10 model
{
  load("G:/1.dyland/output data/CMIP6 interp/hfls.RData")
  load("G:/1.dyland/output data/CMIP6 interp/hfss.RData")
  load("G:/1.dyland/output data/CMIP6 interp/hurs.RData")
  load("G:/1.dyland/output data/CMIP6 interp/sfcWind.RData")
  load("G:/1.dyland/output data/CMIP6 interp/tasmax.RData")
  load("G:/1.dyland/output data/CMIP6 interp/tasmin.RData")
  load("G:/1.dyland/output data/Z_0.5_land.RData")

}



pet <- function(i){
  PET <- llply(1:dim(hfls[[i]])[3],function(x){
    Penman_future(Tmax = tasmax[[i]][,,x] - 273.16, Tmin = tasmin[[i]][,,x] - 273.16,
                  RH = hurs[[i]][,,x], WS = sfcWind[[i]][,,x],
                  SH = hfss[[i]][,,x]*24*60*60/1e6, LH = hfls[[i]][,,x]*24*60*60/1e6, Z = Z )
  }, .progress = 'text')
}

model
pet10 <- pet(i = 10)
pet10 <- abind(pet10, along = 3)
attr(pet10, 'dimnames') <- NULL
save(pet10, file = "G:/1.dyland/output data/PET_CMIP6_his/pet_Amon_FGOALS-g3_historical_r1i1p1f1_gn_185001-185912.RData")



#3.save as nc
load("G:/1.dyland/output data/PET_CMIP6_his/pet_Amon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412.RData")
hfss_nc <- nc_open("H:/CMIP6_data/hfss/historical/hfss_Amon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412.nc")
lon = seq(-179.75,179.75,0.5)
lat = seq(-89.75,89.75,0.5)
time = ncvar_get(hfls_nc, varid = 'time')
lon <- ncdim_def(name = 'lon', "degrees_east", vals = lon)
lat <- ncdim_def(name = 'lat', 'degrees_north', vals = lat)
time <- ncdim_def(name = 'time',calendar = "proleptic_gregorian",units = "days since 1850-01-01",vals = time)
pet <- ncvar_def(name = 'pet',units = 'mm day-1', dim = list(lon, lat, time),missval = 1e+20, longname = 'ACCESS-CM2 PET',prec = 'float')
ncnew <- nc_create(filename = "G:/1.dyland/output data/PET-CMIP6-nc/pet_Amon_ACCESS-CM2_historical_r1i1p1f1_gn_185001-201412.nc",vars = list(pet))
ncvar_put(nc = ncnew, varid = pet,vals = PET1)

nc_close(ncnew)





