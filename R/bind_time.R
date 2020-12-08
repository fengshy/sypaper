bind_time <- function(files,name){
  nums <- which(str_detect(files, paste('_',name,'_', sep = '')))
  data <- llply(nums, function(i){
    data <- nc3(file = files[i], varid = 'time')
    data
  })
  timess <- llply(data, function(x) x[[1]]) %>% unlist()
  date = nc_date(fid = files[nums[1]],unit = NULL,to_char = FALSE)
}
