#!/usr/bin/env Rscript
files = list.files(pattern = '.txt', full.names=T)

assayData = do.call(cbind, lapply(files, function(f){
	df = read.table(f, row.names=1)
	df[seq(nrow(df)-5),, drop=F]
}))

# ss11 ss15 WTCHG_528508_201189.txt -> 201189
# ss89 TSS_P2_E2_13-1234.txt -> P2E2
phenoData <- setNames(plyr::ldply(basename(files), .fun = function(x) {
  if (grepl("ss11", x)) {
    c(strsplit(tools::file_path_sans_ext(x), split = "_")[[1]][[3]], "ss11", 1, "sc")
  } else  if (grepl("ss15", x)) {
    c(strsplit(tools::file_path_sans_ext(x), split = "_")[[1]][[3]], "ss15", 1, "sc")
  } else {
    spl1 = strsplit(tools::file_path_sans_ext(x), split = "_")[[1]]
    c(paste0(spl1[2], strsplit(spl1[3], split = "-")[[1]][1]), "ss8", 1, "sc")
  }
}), c("cell_ID", "timepoint", "replicate_id", "treatment"))

rownames(phenoData) <- phenoData$cell_ID
phenoData$cell_ID <- NULL
colnames(assayData) <- rownames(phenoData) 

write.table(x=assayData, file=paste0('assayData.csv'), sep='\t', row.names=TRUE, quote=FALSE, col.names=NA)
write.table(x=phenoData, file=paste0('phenoData.csv'), sep='\t', row.names=TRUE, quote=FALSE, col.names=NA)
