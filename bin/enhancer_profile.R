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
    output_path = "./output/enhancer_profile/output/"
    # import putative enhancer peaks (ATAC peaks with K27ac - K27me3)
    shared.peaks <- read.delim("./../results/peak_intersect/awk/awk_ATAC.annotated.txt", sep = "\t")
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    output_path = "output/"
    # import putative enhancer peaks (ATAC peaks with K27ac - K27me3)
    shared.peaks <- read.delim("./awk_ATAC.annotated.txt", sep = "\t")
  }
  
  dir.create(output_path, recursive = T)
  
  library(ChIPpeakAnno)
  library(rtracklayer)
}


peaks <- GRanges(seqnames=shared.peaks[,2],
                 ranges=IRanges(start=shared.peaks[,3], 
                                end=shared.peaks[,4], 
                                names=shared.peaks[,1]))


peaks.recentered <- peaks.center <- peaks
start(peaks.center) <- start(peaks) + floor(width(peaks)/2)
width(peaks.center) <- 1
start(peaks.recentered) <- start(peaks.center) - 2000
end(peaks.recentered) <- end(peaks.center) + 2000


ATAC.bw <- import("./../output/results-atac/bwa/mergedLibrary/bigwig/ss8_R1.mLb.clN.bigWig", format="BigWig", which=peaks.recentered, as="RleList")
H3K27Ac.bw <- import("./../output/results-chip/bwa/mergedLibrary/bigwig/ss8-K27Ac_R1.bigWig", format="BigWig", which=peaks.recentered, as="RleList")
H3K27me3.bw <- import("./../output/results-chip/bwa/mergedLibrary/bigwig/ss8-K27me3_R1.bigWig", format="BigWig", which=peaks.recentered, as="RleList")
input.bw <- import("./../output/results-chip/bwa/mergedLibrary/bigwig/ss8-input_R1.bigWig", format="BigWig", which=peaks.recentered, as="RleList")

bw <- list(ATAC = ATAC.bw, H3K27Ac = H3K27Ac.bw, H3K27me3 = H3K27me3.bw, Input = input.bw)

sig <- featureAlignedSignal(bw, peaks.recentered, 
                            upstream=2000, downstream=2000) 


# plot profile around ATAC peaks


png(paste0(output_path, "metaprofile.png"), width=20, height=17, units = 'cm', res = 200)
featureAlignedDistribution(sig, peaks.recentered, upstream=2000, downstream=2000, type="l")
graphics.off()


png(paste0(output_path, "heatmap.png"), width=15, height=15, units = 'cm', res = 200)
featureAlignedHeatmap(sig, peaks.recentered, upstream=2000, downstream=2000, upper.extreme=2.5)
graphics.off()
