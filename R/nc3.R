nc3 <- function(file, varid){
  var_value <- llply(varid, function(x){
    ncvar_get(nc = nc_open(file), x)
  })
  names(var_value) <- varid
  return(var_value)
}
