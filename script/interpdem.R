
#load file and model
load("G:/1.dyland/output data/Z_0.5_land.RData")
lon = seq(-179.75,179.75,0.5)
lat = seq(-89.75,89.75,0.5)
loc = c(which(lon > 0), which(lon < 0))
lon <- lon[loc]
lon[lon < 0] <- lon[lon < 0] +360
Z <- Z[loc,]
lon = seq(0.25, 359.75, 0.5)
setwd("H:/CMIP6_data/")
dem = list()
for(i in 1:length(model)){
  nums <- which(str_detect(filess[[1]], paste(model[i],sep = '')))
  nc <- nc_open(filess[[1]][nums[1]])
  dem[[i]] <- interp.surface.grid(obj = list(x = lon ,y = lat ,z = Z ),
                                  grid.list = list(x = nc$dim$lon$vals,
                                                   y = nc$dim$lat$vals))$z
  print(i)
}
save(dem, file = "G:/1.dyland/output data/")
