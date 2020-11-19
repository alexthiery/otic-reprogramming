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
    output_path = "./output/NF-downstream_analysis/sox8_dea/output/"
    input_file <- "./output/NF-sox8_alignment/featurecounts.merged.counts.tsv"
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    output_path = "output/"
    input_file <- "./featurecounts.merged.counts.tsv"
  }
  
  dir.create(output_path, recursive = T)
  
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

# read in count data and rename columns
read_counts <- read.delim(input_file, stringsAsFactors = FALSE)
colnames(read_counts)[1] <- "gene_id"

# if gene name exists then take gene name, else take ensembl ID and make new name column
read_counts <- read_counts %>% mutate(gene_name = ifelse(!is.na(gene_name), gene_name, gene_id))

# make duplicated gene names unique using "_"
read_counts$gene_name <- make.unique(read_counts$gene_name, sep = "_")

# make gene annotations dataframe
gene_annotations <- read_counts %>% dplyr::select(gene_id, gene_name)

# write CSV for output list
write.csv(read_counts, paste0(output_path, "read_counts.csv"), row.names = F)

# make rownames gene_id and remove ID and names column before making deseq object
rownames(read_counts) <- read_counts$gene_id
read_counts[,1:2] <- NULL

### Add sample group to metadata
col_data <- as.data.frame(sapply(colnames(read_counts), function(x){ifelse(grepl("sox8_oe", x), "Sox8_OE", "Control")}))
colnames(col_data) <- "Group"

### Make deseq object and make Control group the reference level
deseq <- DESeqDataSetFromMatrix(read_counts, design = ~ Group, colData = col_data)
deseq$Group <- droplevels(deseq$Group)
deseq$Group <- relevel(deseq$Group, ref = "Control")

# set plot colours
plot_colours <- list(Group = c(Sox8_OE = "#f55f20", Control = "#957dad"))


### Filter genes which have fewer than 10 readcounts
deseq <- deseq[rowSums(counts(deseq)) >= 10, ]

### Run deseq test - size factors for normalisation during this step are calculated using median of ratios method
deseq <- DESeq(deseq)

### Plot dispersion estimates - dispersion should decrease as counts increase
png(paste0(output_path, "dispersion_est.png"), height = 20, width = 25, units = "cm", res = 400)
plotDispEsts(deseq)
graphics.off()

# LFC shrinkage uses information from all genes to generate more accurate estimates. Specifically, the distribution of
# LFC estimates for all genes is used (as a prior) to shrink the LFC estimates of genes with little information or high
# dispersion toward more likely (lower) LFC estimates.
res <- lfcShrink(deseq, coef="Group_Sox8_OE_vs_Control", type="apeglm")

res$gene_name <- gene_annotations$gene_name[match(rownames(res), gene_annotations$gene_id)]


# plot MA with cutoff for significant genes = padj < 0.05
png(paste0(output_path, "MA_plot.png"), height = 20, width = 25, units = "cm", res = 400)
DESeq2::plotMA(res, alpha = 0.05)
graphics.off()


# Plot volcano plot with padj < 0.05 and abs(fold change) > 1.5 (remove annotation column first)
volc_dat <- as.data.frame(res[,-6])

volc_dat$sig <- apply(volc_dat, 1, function(x) {
  if(!is.na(x["padj"]) & x["padj"]<0.05 & x["log2FoldChange"] > 1.5){
    "upregulated"
  } else if(!is.na(x["padj"]) & x["padj"]<0.05 & x["log2FoldChange"] < -1.5){
    "downregulated"
  } else {"not sig"}
}
)

volc_dat <- volc_dat[order(abs(volc_dat$padj)),]

# add gene name to volcano data
volc_dat$gene <- gene_annotations$gene_name[match(rownames(volc_dat), gene_annotations$gene_id)]

# select genes to add as labels on volcano plot
otic_genes <- c("SOHO-1", "LMX1A", "SOX8", "HOMER2", "DLX3", "ZNF385C", "GATA6", "Six2", "JUN", "PROX1", "HMX1")

downreg <- volc_dat %>%
  dplyr::filter(log2FoldChange < 1.5) %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(gene = as.character(gene)) %>%
  dplyr::filter(!stringr::str_detect(gene, "ENS"))

downreg <- downreg[1:10,"gene"]

labels <- volc_dat[volc_dat$gene %in% c(otic_genes, downreg, "SNAI1"),]

png(paste0(output_path, "volcano.png"), width = 22, height = 16, units = "cm", res = 200)
ggplot(volc_dat, aes(log2FoldChange, -log10(padj))) +
  geom_point(shape=21, aes(colour = sig, fill = sig), size = 0.7) +
  scale_fill_manual(breaks = c("not sig", "downregulated", "upregulated"),
                    values = alpha(c(plot_colours$Group[2], "#c1c1c1", plot_colours$Group[1]), 0.3)) +
  scale_color_manual(breaks = c("not sig", "downregulated", "upregulated"),
                     values= c(plot_colours$Group[2], "#c1c1c1", plot_colours$Group[1])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "top", legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_text_repel(data=labels, size = 3.5, aes(label=gene), segment.color = "black") +
  xlab('log2FC (Sox8_OE - Control)') +
  theme(legend.position = "none")
graphics.off()

################################################################################
# make ordered dataframe for raw counts, normalised counts, and differential expression output
################################################################################

# raw counts dataframe
raw_counts <- as.data.frame(counts(deseq))
colnames(raw_counts) <- paste0("counts_", colnames(raw_counts))
raw_counts$gene_id <- rownames(raw_counts)

# normalised counts dataframe
norm_counts <- as.data.frame(counts(deseq, normalized=TRUE))
colnames(norm_counts) <- paste0("norm_size.adj_", colnames(norm_counts))
norm_counts$gene_id <- rownames(norm_counts)

# differential expression statistics dataframe
DE_res <- as.data.frame(res)
DE_res$gene_id <- rownames(DE_res)

# merge raw_counts, norm_counts and DE_res together into a single dataframe
all_dat <- merge(raw_counts, norm_counts, by = 'gene_id')
all_dat <- merge(all_dat, DE_res, by = 'gene_id')

# move position of gene names column
all_dat <- all_dat[,c(1, ncol(all_dat), 2:{ncol(all_dat)-1})]

# Find which genes are up and downregulated following differential expression analysis
res_up <- all_dat[which(all_dat$padj < 0.05 & all_dat$log2FoldChange > 1.5), ]
res_up <- res_up[order(-res_up$log2FoldChange),]

res_down <- all_dat[which(all_dat$padj < 0.05 & all_dat$log2FoldChange < -1.5), ]
res_down <- res_down[order(res_down$log2FoldChange),]

nrow(res_up)
nrow(res_down)
# 511 genes DE with padj 0.05 & abs(logFC) > 1.5 (399 upregulated, 112 downregulated)

# Write DE data as a csv
res_de <- rbind(res_up, res_down) %>% arrange(-log2FoldChange)

# Write all data as a csv
cat("This table shows the differential expression results for genes with absolute log2FC > 1.5 and adjusted p-value < 0.05 when comparing Sox8 overexpression and control samples (Sox8 - Control)
Reads are aligned to Galgal6 \n
Statistics:
Normalised count: read counts adjusted for library size
pvalue: unadjusted pvalue for differential expression test between Sox8 overexpression and control samples
padj: pvalue for differential expression test between Sox8 overexpression and control samples - adjusted for multiple testing (Benjamini and Hochberg) \n \n",
    file = paste0(output_path, "Supplementary_1.csv"))
write.table(res_de, paste0(output_path, "Supplementary_1.csv"), append=TRUE, row.names = F, na = 'NA', sep=",")


# non-DE genes
res_remain <- all_dat[!rownames(all_dat) %in% rownames(res_up) & !rownames(all_dat) %in% rownames(res_down),]
res_remain <- res_remain[order(-res_remain$log2FoldChange),]

# Make a single dataframe with ordered rows
all_dat <- rbind(res_up, res_down, res_remain)

# Write all data as a csv
cat("This table shows the differential expression results for all genes when comparing Sox8 overexpression and control samples (Sox8 - Control)
Reads are aligned to Galgal6 \n
Statistics:
Normalised count: read counts adjusted for library size
pvalue: unadjusted pvalue for differential expression test between Sox8 overexpression and control samples
padj: pvalue for differential expression test between Sox8 overexpression and control samples - adjusted for multiple testing (Benjamini and Hochberg) \n \n",
    file = paste0(output_path, "Supplementary_2.csv"))
write.table(all_dat, paste0(output_path, "Supplementary_2.csv"), append=TRUE, row.names = F, na = 'NA', sep=",")



####################################################################################
### Plot sample-sample distances and PC plots to show relationship between samples
####################################################################################

# To prevent the highest expressed genes from dominating when clustering we need to rlog (regularised log) transform the data
rld <- rlog(deseq, blind=FALSE)

# Plot sample distance heatmap
sample_dists <- dist(t(assay(rld)))

### check this
sampleDistMatrix <- as.matrix(sample_dists)
rownames(sampleDistMatrix) <- paste(colnames(rld))
colnames(sampleDistMatrix) <- paste(colnames(rld))
colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

png(paste0(output_path, "SampleDist.png"), height = 12, width = 15, units = "cm", res = 400)
pheatmap(sampleDistMatrix, color = colours)
graphics.off()

# Plot sample PCA
png(paste0(output_path, "SamplePCA.png"), height = 12, width = 12, units = "cm", res = 200)
plotPCA(rld, intgroup = "Group") +
  scale_color_manual(values=plot_colours$Group) +
  theme(aspect.ratio=1,
        panel.background = element_rect(fill = "white", colour = "black"))
graphics.off()


# subset genes with padj < 0.05 and abs(LFC) > 1.5
res_sub <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1.5), ]
res_sub <- res_sub[order(-res_sub$log2FoldChange),]

# plot heatmap of DE genes
png(paste0(output_path, "sox8_oe_hm.png"), height = 30, width = 21, units = "cm", res = 200)
pheatmap(assay(rld)[rownames(res_sub),], cluster_rows=T, show_rownames=FALSE,
         show_colnames = F, cluster_cols=T, annotation_col=as.data.frame(colData(deseq)["Group"]),
         annotation_colors = plot_colours, scale = "row", treeheight_row = 0, treeheight_col = 25, cellheight = 1.5, cellwidth = 75)
graphics.off()

#########
# Get biomart GO annotations for TFs
#########

ensembl = useMart("ensembl",dataset="ggallus_gene_ensembl")
TF_subset <- getBM(attributes=c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"),
                   filters = 'ensembl_gene_id',
                   values = rownames(res_sub),
                   mart = ensembl)

# subset genes based on transcription factor GO terms
TF_subset <- TF_subset$ensembl_gene_id[TF_subset$go_id %in% c('GO:0003700', 'GO:0043565', 'GO:0000981')]

res_sub_TF <- res_sub[rownames(res_sub) %in% TF_subset,]


##############################################################
## Save CSV for differentially expressed transcription factors
##############################################################

# subset TFs from all_dat

all_dat_TF <- all_dat[all_dat$gene_id %in% rownames(res_sub_TF),]

cat("This table shows differentially expressed (absolute FC > 1.5 and padj (FDR) < 0.05) transcription factors between Sox8 overexpression and control samples (Sox8 - Control)
Reads are aligned to Galgal6 \n
Statistics:
Normalised count: read counts adjusted for library size
pvalue: unadjusted pvalue for differential expression test between Sox8 overexpression and control samples
padj: pvalue for differential expression test between Sox8 overexpression and control samples - adjusted for multiple testing (Benjamini and Hochberg) \n \n",
    file = paste0(output_path, "Supplementary_3.csv"))
write.table(all_dat_TF, paste0(output_path, "Supplementary_3.csv"), append=TRUE, row.names = F, na = 'NA', sep=",")

##############################################################
# Plot heatmap for differentially expressed transcription factors
##############################################################

rld.plot <- assay(rld)
rownames(rld.plot) <- gene_annotations$gene_name[match(rownames(rld.plot), gene_annotations$gene_id)]


plot <- pheatmap(rld.plot[res_sub_TF$gene_name,], cluster_rows=T, show_rownames=T,
                 show_colnames = F, cluster_cols=T, treeheight_row = 30, treeheight_col = 30,
                 annotation_col=as.data.frame(col_data["Group"]), annotation_colors = plot_colours,
                 scale = "row", main = "Sox8OE enriched TFs (logFC > 1.5, padj = 0.05)", border_color = NA,
                 cellheight = 10, cellwidth = 75)

cat(rownames(rld.plot[res_sub_TF$gene_name,][plot$tree_row[["order"]],]), file=paste0(output_path, "sox8_DE_TFs.txt"))

# plot DE TFs
png(paste0(output_path, "sox8_oe_TFs_hm.png"), height = 20, width = 25, units = "cm", res = 200)
plot
graphics.off()

