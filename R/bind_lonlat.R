bind_lonlat <- function(files,name){
  nums <- which(str_detect(files, paste(name, sep = '')))
  data <- llply(nums, function(i){
    data <- nc3(file = files[i], varid = c('lon', 'lat'))
    data
  })
}
