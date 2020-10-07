BiocManager::install("ChIPpeakAnno")

library(ChIPpeakAnno)



ATAC.peaks <- toGRanges("./../output/results-atac/bwa/mergedLibrary/macs/narrowPeak/ss8_R1.mLb.clN_peaks.narrowPeak")
ChIP.peaks <- toGRanges("./../output/results-chip/bwa/mergedLibrary/macs/broadPeak/ss8-K27Ac_R1_peaks.broadPeak")

data <- list(ATAC = ATAC.peaks, ChIP = ChIP.peaks)

ol <- findOverlapsOfPeaks(data)

#makeVennDiagram(ol)
features <- ol$peaklist[[length(ol$peaklist)]]
wid <- width(features)
feature.recentered <- feature.center <- features
start(feature.center) <- start(features) + floor(wid/2)
width(feature.center) <- 1
start(feature.recentered) <- start(feature.center) - 5000
end(feature.recentered) <- end(feature.center) + 5000








## here we also suggest importData function in bioconductor trackViewer package 
## to import the coverage.
## compare rtracklayer, it will save you time when handle huge dataset.
library(rtracklayer)
ATAC.bw <- import("./../output/results-atac/bwa/mergedLibrary/bigwig/ss8_R1.mLb.clN.bigWig", format="BigWig", which=feature.recentered, as="RleList")
ChIP.bw <- import("./../output/results-chip/bwa/mergedLibrary/bigwig/ss8-K27Ac_R1.bigWig", format="BigWig", which=feature.recentered, as="RleList")

bw <- list(ATAC = ATAC.bw, ChIP = ChIP.bw)

sig <- featureAlignedSignal(bw, feature.center, 
                            upstream=5000, downstream=5000) 


featureAlignedDistribution(sig, feature.center, 
                           upstream=5000, downstream=5000,
                           type="l")





















if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("trackViewer")

library(trackViewer)

gene = "Pax2"
mget(gene, org.Hs.egSYMBOL2EG)




gr <- GRanges("6", IRanges(17836767-5000, 17837601+5000))

K27Ac <- importScore('./../output/results-chip/bwa/mergedLibrary/bigwig/ss8-K27Ac_R1.bigWig', format="BigWig",
                    ranges=gr)


K27Ac$dat <- coverageGR(K27Ac$dat)

viewTracks(trackList(K27Ac), gr=gr, autoOptimizeStyle=TRUE, newpage=FALSE)



dt <- DataTrack(range=K27Ac$dat[strand(K27Ac$dat)=="-"] , 
                genome="hg19", type="hist", name="K27Ac", 
                window=-1, chromosome="chr11", 
                fill.histogram="black", col.histogram="NA",
                background.title="white",
                col.frame="white", col.axis="black",
                col="black", col.title="black")
plotTracks(dt, from=122929275, to=122930122, strand="-")


















BiocManager::install("EnrichedHeatmap")
library(EnrichedHeatmap)
library(rtracklayer)
library(GenomicRanges)


ATAC.bed <- import("./../output/results-atac/bwa/mergedLibrary/macs/narrowPeak/ss8_R1.mLb.clN_summits.bed", format = "bed")
# ATAC.bed <- import("./../testData/atac.bed", format = "bed")
# H3K27ac.bed <- import("./../testData/K27Ac.bed", format = "bed")

H3K27ac.bw <- import("./../output/results-chip/bwa/mergedLibrary/bigwig/ss8-K27Ac_R1.bigWig", format = "BigWig")
ATAC.bw <- import("./../output/results-atac/bwa/mergedLibrary/bigwig/ss8_R1.mLb.clN.bigWig", format = "BigWig")



ATAC.10kb <- resize(ATAC.bed, width = 10000, fix = "center")
ATAC.10kb.center <- resize(ATAC.10kb, width =1, fix = "center")

H3K27ac.mat <- normalizeToMatrix(H3K27ac.bw, ATAC.10kb.center, value_column = "score",
                             mean_mode="w0", w=100, extend = 5000)



ATAC.mat <- normalizeToMatrix(ATAC.bw, ATAC.10kb.center, value_column = "score",
                              mean_mode="w0", w=100, extend = 5000)

## check data range
quantile(H3K27ac.mat, probs = c(0.005, 0.5,0.90))
quantile(ATAC.mat, probs = c(0.005, 0.5,0.90))

## mapping colors
library(circlize)
## from the quantile, I choose the color mapping range
col_fun_ATAC<- circlize::colorRamp2(c(0, 20), c("white", "red"))
col_fun_H3K27ac<- circlize::colorRamp2(c(0, 100), c("white", "red"))


EnrichedHeatmap(ATAC.mat, axis_name_rot = 0, name = "ATAC",
                column_title = "ATAC", use_raster = TRUE, col = col_fun_ATAC,
                top_annotation = HeatmapAnnotation(lines = anno_enriched())) +
  EnrichedHeatmap(H3K27ac.mat, axis_name_rot = 0, name = "H3K27ac",
                  column_title = "H3K27ac", use_raster = TRUE, col = col_fun_H3K27ac,
                  top_annotation = HeatmapAnnotation(lines = anno_enriched()))







# generate a meta profile plot with ChIP-seq data

ATAC_mean<- data.frame(avg = colMeans(ATAC.mat), 
                       CI_lower = apply(ATAC.mat, 2, Hmisc::smean.cl.normal)[2,],
                       CI_upper = apply(ATAC.mat, 2, Hmisc::smean.cl.normal)[3,]) %>% 
  mutate(factor = "ATAC", pos = colnames(ATAC.mat)) 

H3K27ac_mean<- data.frame(avg = colMeans(H3K27ac.mat), 
                          CI_lower = apply(H3K27ac.mat, 2, Hmisc::smean.cl.normal)[2,],
                          CI_upper = apply(H3K27ac.mat, 2, Hmisc::smean.cl.normal)[3,]) %>%
  mutate(factor = "H3K27ac", pos = colnames(H3K27ac.mat))

# install.packages("tidyverse")
library(tidyverse)
combine_both<- bind_rows(ATAC_mean, H3K27ac_mean)

## change position to factor and order it 
combine_both$pos<- factor(combine_both$pos, levels= ATAC_mean$pos)

## without confidence interval

ggplot(combine_both, aes(x = pos,y = avg, group = factor)) + geom_line(aes(color = factor)) + 
  theme_bw(base_size = 14) +
  theme(axis.ticks.x = element_blank()) +
  scale_x_discrete(breaks = c("u1", "d1", "d50"), labels =c ("-5kb", "center", "5kb")) +
  xlab(NULL) + 
  ylab("RPKM")+
  ggtitle("ChIP-seq signal")

## take some touch up to finalize the figure
ggplot(combine_both, aes(x = pos,y = avg, group = factor)) + geom_line(aes(color = factor)) + 
  geom_ribbon(aes(ymin= CI_lower, ymax=CI_upper), alpha=0.2, col = "#8B7E66") +
  theme_bw(base_size = 14) +
  theme(axis.ticks.x = element_blank()) +
  scale_x_discrete(breaks = c("u1", "d1", "d50"), labels =c ("-5kb", "center", "5kb")) +
  xlab(NULL) + 
  ylab("RPKM")+
  ggtitle("ChIP-seq signal")


