#!/usr/bin/env Rscript

library(ChIPpeakAnno)
library(rtracklayer)

output_path = "./output/"
dir.create(output_path, recursive = T)

# import putative enhancer peaks (ATAC peaks with K27ac - K27me3)
shared.peaks <- read.delim("./awk_ATAC.annotated.txt", sep = "\t")

peaks <- GRanges(seqnames=shared.peaks[,2],
                 ranges=IRanges(start=shared.peaks[,3], 
                                end=shared.peaks[,4], 
                                names=shared.peaks[,1]))

# find centre of ATAC peak and get coordinates for +/-2kb
peaks.recentered <- peaks.center <- peaks
start(peaks.center) <- start(peaks) + floor(width(peaks)/2)
width(peaks.center) <- 1
start(peaks.recentered) <- start(peaks.center) - 2000
end(peaks.recentered) <- end(peaks.center) + 2000


# import bigwig files and select regions corresponding to ATAC (putative enhancer) peaks
bigwig_files <- list.files('./', pattern = 'bigWig', full.names = T)
ATAC.bw <- import(bigwig_files[grepl("ss8_R1", bigwig_files)], format="BigWig", which=peaks.recentered, as="RleList")
H3K27Ac.bw <- import(bigwig_files[grepl("H3K27Ac", bigwig_files)], format="BigWig", which=peaks.recentered, as="RleList")
H3K27me3.bw <- import(bigwig_files[grepl("H3K27me3", bigwig_files)], format="BigWig", which=peaks.recentered, as="RleList")
input.bw <- import(bigwig_files[grepl("input", bigwig_files)], format="BigWig", which=peaks.recentered, as="RleList")

# make list of bigwig files
bw <- list(ATAC = ATAC.bw, H3K27Ac = H3K27Ac.bw, H3K27me3 = H3K27me3.bw, Input = input.bw)

# extract signal for +/-2kb around enhancer peak for visualisation
sig <- featureAlignedSignal(bw, peaks.recentered, 
                            upstream=2000, downstream=2000) 


# plot profile around ATAC peaks
png(paste0(output_path, "metaprofile.png"), width=20, height=17, units = 'cm', res = 200)
featureAlignedDistribution(sig, peaks.recentered, upstream=2000, downstream=2000, type="l")
graphics.off()

# plot heatmap
png(paste0(output_path, "heatmap.png"), width=15, height=15, units = 'cm', res = 200)
featureAlignedHeatmap(sig, peaks.recentered, upstream=2000, downstream=2000, upper.extreme=2.5)
graphics.off()
