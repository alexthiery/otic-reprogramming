# cbind dataframes by rowname
cbind_rownames = function(...){
  test <- list(...)
  
  first.names = rownames(test[[1]])

  if(!all(unlist(lapply(test, function(x){
    all(rownames(x) %in% first.names)
  })))){
    stop("rownames in supplied dataframes are not equal")
  }else{
    
    test <- lapply(test, function(x){
      x[match(rownames(x), first.names),]
      x
    })
    do.call(cbind, test)
  }
}


