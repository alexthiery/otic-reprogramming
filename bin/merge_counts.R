#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set run location
if(length(commandArgs(trailingOnly = TRUE)) == 0){
  cat('No command line arguments provided, user defaults paths are set for running interactively in Rstudio on docker\n')
  opt$runtype = "user"
} else {
  if(is.null(opt$runtype)){
    stop("--runtype must be either 'user' or 'nextflow'")
  }
  if(tolower(opt$runtype) != "user" & tolower(opt$runtype) != "nextflow"){
    stop("--runtype must be either 'user' or 'nextflow'")
  }
}

# Set paths and load data
{
  if (opt$runtype == "user"){
    output_path = "./output/merge_counts/output/"
    input_files <- list.files("./testData/process_counts/", pattern = "*.txt", full.names = T)
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    output_path = "output/"
    input_files <- list.files("./", pattern = "*.txt", full.names = T)
  }
  dir.create(output_path, recursive = T)
  
  library(plyr)
}

assayData = do.call(cbind, lapply(input_files, function(f){
	df = read.table(f, row.names=1)
	df[seq(nrow(df)-5),, drop=F]
}))

# ss11 ss15 WTCHG_528508_201189.txt -> 201189
# ss89 TSS_P2_E2_13-1234.txt -> P2E2
phenoData <- setNames(ldply(basename(input_files), .fun = function(x) {
  if (grepl("ss11", x)) {
    c(strsplit(tools::file_path_sans_ext(x), split = "_")[[1]][[3]], 11, 1, "sc")
  } else  if (grepl("ss15", x)) {
    c(strsplit(tools::file_path_sans_ext(x), split = "_")[[1]][[3]], 15, 1, "sc")
  } else {
    spl1 = strsplit(tools::file_path_sans_ext(x), split = "_")[[1]]
    c(paste0(spl1[2], strsplit(spl1[3], split = "-")[[1]][1]), 8, 1, "sc")
  }
}), c("cell_ID", "timepoint", "replicate_id", "treatment"))

rownames(phenoData) <- phenoData$cell_ID
phenoData$cell_ID <- NULL
colnames(assayData) <- rownames(phenoData) 

# make dataframe for gfp counts
gfpData <- assayData["GFP",]

# remove gfp counts from assayData
assayData <- assayData[-which(rownames(assayData) == "GFP"),]

write.table(x=gfpData, file=paste0(output_path, 'gfpData.csv'), sep='\t', row.names=TRUE, quote=FALSE, col.names=NA)
write.table(x=assayData, file=paste0(output_path, 'assayData.csv'), sep='\t', row.names=TRUE, quote=FALSE, col.names=NA)
write.table(x=phenoData, file=paste0(output_path, 'phenoData.csv'), sep='\t', row.names=TRUE, quote=FALSE, col.names=NA)
