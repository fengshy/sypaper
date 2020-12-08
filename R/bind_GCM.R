bind_GCM <- function(files, name){
  variable <- strsplit(files[1], '/',fixed = TRUE)[[1]][1]
  nums <- which(str_detect(files, paste('_',name,'_', sep = '')))
  data <- llply(nums, function(i){
    data <- nc3(file = files[i], varid = c('lon', 'lat', variable))
    data
  })
  var_value <- llply(data, function(x) x[[3]]) %>% abind(., along = 3)
  attr(var_value, 'dimnames') <- NULL
  return(list(lon = data[[1]]$lon, lat = data[[1]]$lat, variable = var_value))
  rm(data); gc()
}
# read lon lat and variable from a or many nc
