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
    input_files <- list.files("./output/results-sox8_oe/featureCounts", pattern = "*.txt$", full.names = T)
    
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

# samples WTCHG_706842_265195 and WTCHG_706842_266183 were not imported for downstream analysis due to low number of reads (see multiqc output)

# make dictionary of sample names and readcount files
# samples are renamed from original sample sheet for downstream analysis. Original sample IDs can be found under Sample_ReadMe.txt
samples_IDs <- data.frame(Sample = c("Sox8OE1", "Sox8OE2", "Sox8OE3", "Control1", "Control2", "Control3"),
                          ID = c("WTCHG_706842_201108", "WTCHG_706842_202120", "WTCHG_706842_203132", "WTCHG_706842_205156",
                                 "WTCHG_706842_207180", "WTCHG_706842_208192"),
                          stringsAsFactors = F)

# read in count data and rename columns
read_counts <- as.data.frame(read.table(input_files, header = T, stringsAsFactors = F))
colnames(read_counts) <- unlist(lapply(colnames(read_counts), function(x) ifelse(grepl("WTCHG", x), samples_IDs$Sample[grepl(gsub("_1Aligned.*", "", x ), samples_IDs$ID)], x)))
colnames(read_counts)[1] <- "gene_id"

# remove samples with poor quality
read_counts <- read_counts[,!is.na(colnames(read_counts))]


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
col_data <- as.data.frame(sapply(colnames(read_counts), function(x){ifelse(grepl("Sox8", x), "Sox8", "Control")}))
colnames(col_data) <- "Group"

### Make deseq object and make Control group the reference level
deseq <- DESeqDataSetFromMatrix(read_counts, design = ~ Group, colData = col_data)
deseq$Group <- droplevels(deseq$Group)
deseq$Group <- relevel(deseq$Group, ref = "Control")

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
res <- lfcShrink(deseq, coef="Group_Sox8_vs_Control", type="apeglm")

res$gene_name <- gene_annotations$gene_name[match(rownames(res), gene_annotations$gene_id)]


# plot MA with cutoff for significant genes = padj < 0.05
png(paste0(output_path, "MA_plot.png"), height = 20, width = 25, units = "cm", res = 400)
DESeq2::plotMA(res, alpha = 0.05)
graphics.off()


# Plot volcano plot with padj < 0.05 and abs(fold change) > 1.5 (remove annotation column first)
volc_dat <- as.data.frame(res[,-6])
volc_dat$sig <- apply(volc_dat, 1, function(x) if(is.na(x["padj"]) | x["padj"]>=0.05 | abs(x["log2FoldChange"]) <=1.5){
  "Not sig"
} else{"Differentially expressed (padj <0.05, absolute log2 FC >1.5"}
)

volc_dat <- volc_dat[order(abs(volc_dat$padj)),]

# add gene name to volcano data
volc_dat$gene <- gene_annotations$gene_name[match(rownames(volc_dat), gene_annotations$gene_id)]

# select genes to add as labels on volcano plot
labels <- head(volc_dat, 50)
labels <- labels[!grepl("ENSGAL", labels$gene),]

# Get biomart GO annotations for TFs
ensembl = useMart("ensembl",dataset="ggallus_gene_ensembl")
TF_subset <- getBM(attributes=c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"),
                   filters = 'ensembl_gene_id',
                   values = rownames(res_sub),
                   mart = ensembl)

# subset genes based on transcription factor GO terms
TF_subset <- TF_subset$ensembl_gene_id[TF_subset$go_id %in% c('GO:0003700', 'GO:0043565', 'GO:0000981')]






png(paste0(output_path, "volcano.png"), width = 22, height = 16, units = "cm", res = 200)
ggplot(volc_dat, aes(log2FoldChange, -log10(padj))) +
  geom_point(shape=21, aes(colour = sig, fill = sig), size = 0.7) +
  scale_fill_manual(breaks = c("Not sig", "Differentially expressed (padj <0.05, absolute log2 FC >1.5"),
                    values= alpha(c("gray40", "red"), 0.3)) +
  scale_color_manual(breaks = c("Not sig", "Differentially expressed (padj <0.05, absolute log2 FC >1.5"),
                     values=c("gray40", "red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "top", legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_text_repel(data=labels, size = 2.5, aes(label=gene), segment.color = "gray80") +
  geom_vline(xintercept = 1.5, linetype="dashed",
             color = "gray20", size=0.4) +
  geom_vline(xintercept = -1.5, linetype="dashed",
             color = "gray20", size=0.4) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed",
             color = "gray20", size=0.4) +
  theme(legend.position = "none")
graphics.off()



png(paste0(output_path, "volcano.png"), width = 22, height = 16, units = "cm", res = 200)
ggplot(volc_dat, aes(log2FoldChange, -log10(padj))) +
  geom_point(shape=21, aes(colour = sig, fill = sig), size = 0.7) +
  scale_fill_manual(breaks = c("Not sig", "Differentially expressed (padj <0.05, absolute log2 FC >1.5"),
                    values= alpha(c("gray40", "red"), 0.3)) +
  scale_color_manual(breaks = c("Not sig", "Differentially expressed (padj <0.05, absolute log2 FC >1.5"),
                     values=c("gray40", "red")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_text_repel(data=head(volc_dat, 40), size = 2.5, aes(label=gene), segment.color = "gray80")
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

# 511 genes DE with padj 0.05 & abs(logFC) > 1.5 (398 upregulated, 112 downregulated)

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
    file = paste0(output_path, "Supplementary_1.csv"))
write.table(all_dat, paste0(output_path, "Supplementary_1.csv"), append=TRUE, row.names = F, na = 'NA', sep=",")

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
  theme(aspect.ratio=1,
        panel.background = element_rect(fill = "white", colour = "black"))
graphics.off()


# subset genes with padj < 0.05 and abs(LFC) > 1.5
res_sub <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1.5), ]
res_sub <- res_sub[order(-res_sub$log2FoldChange),]

# plot heatmap of DE genes
png(paste0(output_path, "Sox8OE_HM.png"), height = 30, width = 21, units = "cm", res = 200)
pheatmap(assay(rld)[rownames(res_sub),], cluster_rows=T, show_rownames=FALSE,
         cluster_cols=T, annotation_col=as.data.frame(colData(deseq)["Group"]),
         scale = "row", treeheight_row = 20, treeheight_col = 40)
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
    file = paste0(output_path, "Supplementary_2.csv"))
write.table(all_dat_TF, paste0(output_path, "Supplementary_2.csv"), append=TRUE, row.names = F, na = 'NA', sep=",")

##############################################################
# Plot heatmap for differentially expressed transcription factors
##############################################################

rld.plot <- assay(rld)
rownames(rld.plot) <- gene_annotations$gene_name[match(rownames(rld.plot), gene_annotations$gene_id)]

# plot DE TFs
png(paste0(output_path, "Sox8OE_TFs_HM.png"), height = 25, width = 20, units = "cm", res = 200)
pheatmap(rld.plot[res_sub_TF$gene_name,], cluster_rows=T, show_rownames=T,
         show_colnames = F, cluster_cols=T, treeheight_row = 30, treeheight_col = 30,
         annotation_col=as.data.frame(col_data(deseq)["Group"]), scale = "row",
         main = "Sox8OE enriched TFs (logFC > 1.5, padj = 0.05)",
         border_color = NA)
graphics.off()


###########################################################################
## Compare DE data with DE TFs from Chen et al. (2017) Development
###########################################################################

# Download supplementary file 4 from Chen et al. (2017)
# OEP vs PPR 5/6ss, PPR 5/6ss vs PPR 8/9ss, PPR 8/9ss vs PPR 11/12ss. Dataset is already filtered for transcription factors

temp <- tempfile()
download.file("http://www.biologists.com/DEV_Movies/DEV148494/TableS4.xlsx", temp)

# read xlsx file
otic_enr <- read.xlsx(temp, startRow = 15)
unlink(temp)

# assign column names
colnames(otic_enr)[1:13] <- paste0(colnames(otic_enr)[1:13], c(rep('_normalised_count', 4), rep('_foldChange', 3),
                                                               rep('_pval', 3), rep('_padj', 3)))

# assign row names
rownames(otic_enr) <- otic_enr[,17]
otic_enr[,17] <- NULL

# remove genes from Chen dataset which are not DE in at least one of the stages (not sure why these are in the supplementary file)
otic_enr <- otic_enr[apply(otic_enr, 1, function(x) any(!is.na(x[c("5-6ss_foldChange", "8-9ss_foldChange", "11-12ss_foldChange")]))),]


# compare Sox8OE data vs otic enriched
# subset genes wich are 1.5FC between either PPR vs 4/5ss, 4/5ss vs 8/9ss, 8/9ss vs 11/12ss
otic_enr <- otic_enr[otic_enr$`5-6ss_foldChange` > 1.5 |
                       otic_enr$`8-9ss_foldChange` > 1.5 |
                       otic_enr$`11-12ss_foldChange` > 1.5,]


#plot venn diagram comparing OOPE and Sox8OE
venn.diagram(list(otic_enr=rownames(otic_enr), sox8OE=res_sub_TF$gene_name),
             category.names = c("logFC > 1.5 in any of PPR vs 5/6ss, 5/6ss vs 8/9ss, 8/9ss vs 11/12ss \n(Chen et al. 2017)", "Sox8OE enriched TFs\n(logFC > 1.5, padj = 0.05)"),
             filename = paste0(output_path, "OticEnr.vs.Sox8OE.png"),
             output = TRUE,
             imagetype = "png",
             height = 1100,
             width = 1100,
             resolution = 600,
             compression = "lzw",
             lwd = 1,
             col=c("#440154ff", '#21908dff'),
             fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
             cex = 0.5,
             fontfamily = "sans",
             cat.cex = 0.25,
             cat.pos = c(0, 0),
             cat.dist = c(0.03, 0.03),
             cat.fontfamily = "sans",
             cat.col = c("#440154ff", '#21908dff'),
             ext.percent = 0
)

# make csv file of genes in each part of venn diagram
# when identifying shared genes between current study and previous studies which are alligned using different genome version, match genes using gene name.
# this is important as some of the Ensembl IDs from past genome versions have been depracated and therefore are absent from our data.
venn.genes <- list("Otic enriched" = rownames(otic_enr)[!rownames(otic_enr) %in% res_sub_TF$gene_name],
                   "Sox8OE" = res_sub_TF$gene_name[!res_sub_TF$gene_name %in% rownames(otic_enr)],
                   "Shared" = res_sub_TF$gene_name[res_sub_TF$gene_name %in% rownames(otic_enr)])

venn.genes.df <- t(plyr::ldply(venn.genes, rbind))
colnames(venn.genes.df) <- venn.genes.df[1,]
venn.genes.df <- venn.genes.df[-1,]


cat("This table provides a list of genes from each part of the venn diagram in Fig.x.
Differentially expressed transcription factors between Sox8 overexpression and control samples were cross compared with genes in supplementary table 4 of Chen et al. (2017) Development
Genes from supplementary table 4 of Chen et al. (2017) Development were filtered and kept if they were found to be differentially expressed (absolute FC > 1.5) between either: PPR vs 5/6ss otic; 5/6ss otic vs 8/9ss otic; 8/9ss otic vs 11/12ss otic \n \n",
    file = paste0(output_path, "Supplementary_3.csv"))
write.table(venn.genes.df, paste0(output_path, "Supplementary_3.csv"), append=TRUE, row.names = F, na = '', sep=",")



# plot heatmap
rld.plot <- assay(rld)
rownames(rld.plot) <- gene_annotations$gene_name[match(rownames(rld.plot), gene_annotations$gene_id)]

png(paste0(output_path, "anyOticvsSox8OE.png"),height = 8, width = 21, units = "cm", res = 200)
pheatmap(rld.plot[venn.genes$Shared,], cluster_rows=T, show_rownames=T,
         show_colnames = F, cluster_cols=T, treeheight_row = 30, treeheight_col = 30,
         annotation_col=as.data.frame(col_data(deseq)["Group"]), scale = "row",
         main = "Shared Otic and Sox8OE enriched TFs \n(logFC > 1.5, padj = 0.05)", cellwidth = 50, cellheight = 10,
         border_color = NA)
graphics.off()


# plot all transcription factors DE in Chen et al. 2017 abs(1.5 FC) using our data - do not filter genes which are not DE in the Sox8OE

png(paste0(output_path, "otic_enr.heatmap.png"),height = 50, width = 21, units = "cm", res = 200)
pheatmap(rld.plot[rownames(otic_enr)[rownames(otic_enr) %in% rownames(rld.plot)],], cluster_rows=T, show_rownames=T,
         show_colnames = F, cluster_cols=T, treeheight_row = 30, treeheight_col = 30,
         annotation_col=as.data.frame(col_data(deseq)["Group"]), scale = "row",
         main = "Otic enriched TFs - not necessarily DE between Sox8 and control \n(logFC > 1.5, padj = 0.05)",
         border_color = NA)
graphics.off()


# This heatmap reveals that although many genes are not statistically DE in the Sox8OE - they are clearly upregulated in two of the three Sox8 samples.
# A possible explanation for this is that the Sox8OE does not necessarily switch on the otic program at the same rate in different samples and different cells.
# There may be modules of genes which are switched on at different points of otic specification. This variation could explain why these
# genes are not statistically differentially expressed.

################################################################################
# save CSV norm counts and Sox8OE DEA for genes from Chen et al. 2017
################################################################################

all_dat_Chen_DE <- all_dat[all_dat$gene_name %in% rownames(otic_enr),]


cat("This table shows genes subset from supplementary table 4 of Chen et al. (2017) Development
These genes were found to be differentially expressed (absolute FC > 1.5) between either: PPR vs 5/6ss otic; 5/6ss otic vs 8/9ss otic; 8/9ss otic vs 11/12ss otic
The data presented in this table are from Sox8 overexpression and control samples
Genes presented in this table are not necessarily differentially expressed between Sox8 overexpression and control samples \n
Statistics:
Normalised count: read counts adjusted for library size
pvalue: unadjusted pvalue for differential expression test between Sox8 overexpression and control samples
padj: pvalue for differential expression test between Sox8 overexpression and control samples - adjusted for multiple testing (Benjamini and Hochberg) \n \n",
    file = paste0(output_path, "Supplementary_4.csv"))
write.table(all_dat_Chen_DE, paste0(output_path, "Supplementary_4.csv"), append=TRUE, row.names = F, na = 'NA', sep=",")




###########

# comparing with NC RNAseq

# Load data from Williams et al. 2019 Developmental Cell
williams_data <- tempfile()
download.file("https://ndownloader.figshare.com/files/12750314", williams_data)
load(williams_data)




#########
# Get biomart GO annotations for TFs
#########

ensembl = useMart("ensembl",dataset="ggallus_gene_ensembl")
TF_subset <- getBM(attributes=c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003"),
                   filters = 'ensembl_gene_id',
                   values = de_5to6$EnsemblID,
                   mart = ensembl)

# subset genes based on transcription factor GO terms
TF_subset <- TF_subset$ensembl_gene_id[TF_subset$go_id %in% c('GO:0003700', 'GO:0043565', 'GO:0000981')]




# Filter transcription factors from NC 5-6ss which are abs(log2FC > 1.5) & padj < 0.05 in citrine+ cells vs citrine- cells
de_5to6 <- de_5to6[abs(de_5to6$log2FoldChange) > 1.5 &
                     !is.na(de_5to6$log2FoldChange) &
                     de_5to6$padj < 0.05 &
                     !is.na(de_5to6$padj) &
                     de_5to6$EnsemblID %in% TF_subset,]

# Filter transcription factors from NC 8-10ss which are abs(log2FC > 1.5) & padj < 0.05 in citrine+ cells vs citrine- cells
de_8to10 <- de_8to10[abs(de_8to10$log2FoldChange) > 1.5 &
                       !is.na(de_8to10$log2FoldChange) &
                       de_8to10$padj < 0.05 &
                       !is.na(de_8to10$padj) &
                       de_8to10$EnsemblID %in% TF_subset,]


# List of TFs at 5-6ss OR 8-10ss
NC_enriched_TFs <- unique(c(de_5to6$EnsemblID, de_8to10$EnsemblID))

# Transcription factors which are differentially expressed in Sox8OE vs control samples and in differentially expressed in NC
SOX8_NC_shared <- rownames(res_sub)[rownames(res_sub) %in% NC_enriched_TFs]


# Plot heatmap

NC_TF_DE_dat <- assay(rld)
NC_TF_DE_dat <- NC_TF_DE_dat[rownames(NC_TF_DE_dat) %in% SOX8_NC_shared,]
rownames(NC_TF_DE_dat) <- gene_annotations$gene_name[match(rownames(NC_TF_DE_dat), gene_annotations$gene_id)]

png(paste0(output_path, "SOX8_NC_shared.heatmap.ens.png"),height = 10, width = 21, units = "cm", res = 200)
pheatmap(NC_TF_DE_dat, cluster_rows=T, show_rownames=T,
         show_colnames = F, cluster_cols=F, treeheight_row = 30, treeheight_col = 30,
         annotation_col=as.data.frame(col_data(deseq)["Group"]), scale = "row",
         main = "Shared Sox8OE enriched and NC enriched TFs \n(logFC > 1.5, padj = 0.05)",
         border_color = NA)
graphics.off()


# plot all transcription factors DE in Williams et al. 2019 using our data - do not filter genes which are NOT DE in the Sox8OE

NC_TF_dat <- assay(rld)
NC_TF_dat <- NC_TF_dat[rownames(NC_TF_dat) %in% NC_enriched_TFs,]
rownames(NC_TF_dat) <- gene_annotations$gene_name[match(rownames(NC_TF_dat), gene_annotations$gene_id)]

png(paste0(output_path, "NC_enriched_TFs.heatmap.png"),height = 25, width = 21, units = "cm", res = 200)
pheatmap(NC_TF_dat, cluster_rows=T, show_rownames=T,
         show_colnames = F, cluster_cols=F, treeheight_row = 30, treeheight_col = 30,
         annotation_col=as.data.frame(col_data(deseq)["Group"]), scale = "row",
         main = "NC enriched TFs - not necessarily DE between Sox8 and control \n(logFC > 1.5, padj = 0.05)",
         border_color = NA)
graphics.off()


################################################################################
# save CSV of norm counts and Sox8OE DEA for TFs from Williams et al. 2019
################################################################################

all_dat_NC_TF <- all_dat[all_dat$gene_id %in% NC_enriched_TFs,]

cat("This table shows genes subset from Williams et al. (2019) Developmental Cell bulk RNAseq data
Bulk RNAseq was carried out on citrine positive and citrine negative cells following electroporation of the NC specific enhancer NC1
Sequencing was done at two stages: 5-6ss and 8-10ss
Genes are filtered for transcription factors differentially expressed between citrine positive and negative samples (absolute FC > 1.5 and padj < 0.05) at either 5-6ss or 8-10ss
The data presented in this table are from Sox8 overexpression and control samples
Genes presented in this table are not necessarily differentially expressed between Sox8 overexpression and control samples \n
Statistics:
Normalised count: read counts adjusted for library size
pvalue: unadjusted pvalue for differential expression test between Sox8 overexpression and control samples
padj: pvalue for differential expression test between Sox8 overexpression and control samples - adjusted for multiple testing (Benjamini and Hochberg) \n \n",
    file = paste0(output_path, "Supplementary_5.csv"))
write.table(all_dat_NC_TF, paste0(output_path, "Supplementary_5.csv"), append=TRUE, row.names = F, na = 'NA', sep=",")









