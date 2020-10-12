#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
spec = matrix(c(
  'runtype', 'l', 2, "character"
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
    output_path = "./output/sox8_oe/plots/"
    rds_path = "./output/sox8_oe/rds_files/"
    input_files <- list.files("./output/results-lmx1a/featureCounts", pattern = "*.txt$", full.names = T)
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    sapply(list.files(opt$custom_functions, full.names = T), source)
    output_path = "plots/"
    rds_path = "rds_files/"
    input_files <- list.files("./", pattern = "*.txt$", full.names = T)
  }
  
  dir.create(output_path, recursive = T)
  dir.create(rds_path, recursive = T)
  
  library(biomaRt)
  library(tidyverse)
  library(readr)
  library(RColorBrewer)
  library(VennDiagram)
  
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
  library(DESeq2)
  library(apeglm)
  library(openxlsx)
}