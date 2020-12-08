F_area <- function(lonlat){
  grid <- SpatialPixelsDataFrame(points = lonlat, data = data.frame(id = 1:nrow(lonlat)),
                                 proj4string = CRS('+proj=longlat +ellps=WGS84'))
  area <- raster::area(raster(grid))# area weights
  area <- area@data@values
  return(area)
}