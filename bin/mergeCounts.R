#!/usr/bin/env Rscript

files = list.files('*.txt', full.names=T)

dataset = do.call(cbind, lapply(files, function(f){
	df = read.table(f, row.names=1)
	df[seq(nrow(df)-5),, drop=F]
}))

# ss11 ss15 WTCHG_528508_201189.txt -> 201189
# ss89 TSS_P2_E2_13-1234.txt -> P2E2
colnames(dataset) <- unlist(lapply(basename(files), function(x) {
  if (grepl("ss11", x) | grepl("ss15", x)) {
    strsplit(tools::file_path_sans_ext(x), split = "_")[[1]][[3]]
  } else{
    spl1 = strsplit(tools::file_path_sans_ext(x), split = "_")[[1]]
    paste0(spl1[2], strsplit(spl1[3], split = "-")[[1]][1])
  }
}))

write.table(x=dataset, file=paste0('readcounts.csv'), sep='\t', row.names=TRUE, quote=FALSE, col.names=NA)
