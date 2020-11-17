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
    output_path = "./output/merge_smartseq_bulk/output/"
    
    assay <- read.csv("./results/scRNAseq_alignment/merged_counts/output/assayData.csv", sep = '\t',  row.names = 1, check.names = FALSE)
    pheno <- read.csv("./results/scRNAseq_alignment/merged_counts/output/phenoData.csv", sep = '\t',  row.names = 1, check.names = FALSE)
    
    # read in count data
    read_counts <- as.data.frame(read.table("./output/results-sox8_oe/featureCounts/merged_gene_counts.txt", header = T, stringsAsFactors = F, row.names = 1))
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    output_path = "output/"
    
    assay <- read.csv("./assayData.csv", sep = '\t',  row.names = 1, check.names = FALSE)
    pheno <- read.csv("./phenoData.csv", sep = '\t',  row.names = 1, check.names = FALSE)
    
    # read in count data
    read_counts <- as.data.frame(read.table("./featureCounts/merged_gene_counts.txt", header = T, stringsAsFactors = F, row.names = 1))
  }
  dir.create(output_path, recursive = T)
  
  library(plyr)
}


# make dictionary of sample names and readcount files
# samples are renamed from original sample sheet for downstream analysis. Original sample IDs can be found under Sample_ReadMe.txt
samples_IDs <- data.frame(Sample = c("sox81", "sox82", "sox83", "control1", "control2", "control3"),
                          ID = c("WTCHG_706842_201108", "WTCHG_706842_202120", "WTCHG_706842_203132", "WTCHG_706842_205156",
                                 "WTCHG_706842_207180", "WTCHG_706842_208192"),
                          stringsAsFactors = F)


# remove gene name col
read_counts <- read_counts[,-1]

# rename columns
colnames(read_counts) <- unlist(lapply(colnames(read_counts), function(x) samples_IDs$Sample[grepl(gsub("_1Aligned.*", "", x ), samples_IDs$ID)]))
bulk_pheno <- data.frame(timepoint = rep(0, 6), replicate_id = rep(1, 6), treatment = rep("bulk", 6), row.names = colnames(read_counts))

if(setequal(rownames(assay), rownames(read_counts))){
  read_counts <- read_counts[match(rownames(assay), rownames(read_counts)),]
  combined_assay <- cbind(read_counts, assay)
  combined_pheno <- rbind(bulk_pheno, pheno)
} else {stop("Bulk and single cell data contain different gene names. Check that samples are aligned to the same genome version.")}

write.table(x=combined_assay, file=paste0(output_path, 'assayData_bulk.csv'), sep='\t', row.names=TRUE, quote=FALSE, col.names=NA)
write.table(x=combined_pheno, file=paste0(output_path, 'phenoData_bulk.csv'), sep='\t', row.names=TRUE, quote=FALSE, col.names=NA)





