#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
spec = matrix(c(
  'runtype', 'l', 2, "character",
  'cores'   , 'c', 2, "integer",
  'custom_functions', 'm', 2, "character"
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
  if(tolower(opt$runtype) == "nextflow"){
    if(is.null(opt$custom_functions) | opt$custom_functions == "null"){
      stop("--custom_functions path must be specified in process params config")
    }
  }
}

# Set paths and load data
{
  if (opt$runtype == "user"){
    sapply(list.files('./bin/custom_functions/', full.names = T), source)
    
    output_path = "./output/antler"
    plot_path = "./output/antler/plots/"
    rds_path = "./output/antler/rds_files/"
    merged_counts_path = './output/merged_counts/'
    genome_annotations_path = './output/extract_gtf_annotations'
    gfp_counts = './output/merged_counts/'
    
    ncores = 8
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')
    
    sapply(list.files(opt$custom_functions, full.names = T), source)
    output_path = "./output/"
    plot_path = "./output/plots/"
    rds_path = "./output/rds_files/"
    merged_counts_path = './'
    genome_annotations_path = './'
    gfp_counts = './'
    
    ncores = opt$cores
  }
  
  dir.create(output_path, recursive = T)
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
  
  # loadload required packages
  library(Antler)
  library(velocyto.R)
  library(stringr)
  library(monocle)
  library(plyr)
  library(dplyr)
  library(ggsignif)
  library(cowplot)
  library(rstatix)
}

# set pipeline params
seed=1
perp=5
eta=200

#' Stage colors
stage_cols = setNames(c("#BBBDC1", "#6B98E9", "#05080D"), c('8', '11', '15'))

#' # Load and hygienize dataset
m = Antler$new(plot_folder=plot_path, num_cores=ncores)

# load in phenoData and assayData from ../dataset -> assayData is count matrix; phenoData is metaData (i.e. replicated, conditions, samples etc)
m$loadDataset(folderpath=merged_counts_path)

pData(m$expressionSet)$stage_colors = stage_cols[as.character(pData(m$expressionSet)$timepoint)]

m$plotReadcountStats(data_status="Raw", by="timepoint", category="timepoint", basename="preQC", reads_name="read", cat_colors=unname(stage_cols))

# set gene names
# read in annotations file
gtf_annotations = read.csv(list.files(genome_annotations_path, full.names = T), stringsAsFactors = F)

# add missing annotations to annotations file
extra_annotations = c('FOXI3' = 'ENSGALG00000037457', 'ATN1' = 'ENSGALG00000014554', 'TBX10' = 'ENSGALG00000038767',
                      'COL11A1' = 'ENSGALG00000005180', 'GRHL2' = 'ENSGALG00000037687')

# add extra annotations to annotations csv file
gtf_annotations[,2] <- apply(gtf_annotations, 1, function(x) ifelse(x[1] %in% extra_annotations, names(extra_annotations)[extra_annotations %in% x], x[2]))

# label extra MT genes
MT_genes = c('ND3', 'CYTB', 'COII', 'ATP8', 'ND4', 'ND4L')
gtf_annotations[gtf_annotations[,2] %in% MT_genes,2] <- paste0('MT-', gtf_annotations[gtf_annotations[,2] %in% MT_genes,2])

write.csv(gtf_annotations, paste0(output_path, 'new_annotations.csv'), row.names = F)

# set gene annotations
m$setCurrentGeneNames(geneID_mapping_file=paste0(output_path, 'new_annotations.csv'))



#' Store known genes 
apriori_genes = c(
  'DACH1', 'DLX3', 'DLX5', 'DLX6', 'EYA1', 'EYA2', 'FOXG1', 'FOXI3', 'GATA3', 'GBX2', 'HESX1', 'IRX1', 'IRX2', 'IRX3', 'LMX1A', 'LMX1B', 'PAX2', 'SALL4', 'SIX4', 'SOHO-1', 'SOX10', 'SOX2', 'SOX8', 'SOX9', 'SIX1', 'TBX2', # Known_otic_genes
  'ATN1', 'BACH2', 'CNOT4', 'CXCL14', 'DACH2', 'ETS1', 'FEZ1', 'FOXP3', 'HIPK1', 'HOMER2', 'IRX4', 'IRX5', 'KLF7', 'KREMEN1', 'LDB1', 'MSI1', 'PDLIM4', 'PLAG1', 'PNOC', 'PKNOX2', 'RERE', 'SMOC1', 'SOX13', 'TCF7L2', 'TEAD3', 'ZBTB16', 'ZFHX3', 'ZNF384', 'ZNF385C', # New_otic_TFs
  'PRDM1', 'FOXI1', 'NKX2-6', 'NR2F2', 'PDLIM1', 'PHOX2B', 'SALL1', 'TBX10', 'TFAP2E', 'TLX1', 'VGLL2', # Epibranchial_genes
  'ZNF423', 'CXCR4', 'MAFB', 'MYC', 'ZEB2', # Neural_Genes
  'CD151', 'ETS2', 'FGFR4', 'OTX2', 'Pax3', 'PAX6', 'PAX7', 'SIX3', # Non_otic_placode_genes
  'FOXD3', 'ID2', 'ID4', 'MSX1', 'TFAP2A', 'TFAP2B', 'TFAP2C', # Neural_Crest_Genes
  'GATA2', # Epidermis_genes
  'COL11A1', 'DTX4', 'GRHL2', 'NELL1', 'OTOL1', # Disease_associated_genes
  'ARID3A', 'BMP4', 'CREBBP', 'ETV4', 'ETV5', 'EYA4', 'FOXP4', 'FSTL4', 'HOXA2', 'JAG1', 'LFNG', 'LZTS1', 'MAFA', 'MEIS1', 'MYB', 'MYCN', 'NFKB1', 'NOTCH1', 'SPRY1', 'SPRY2', 'SSTR5', # chen et al. 2017
  'TWIST1', 'ASL1', 'MAF' # other_genes
)

m$favorite_genes <- unique(sort(apriori_genes))

#' Remove outliers genes and cells
m$removeOutliers( lowread_thres = 5e5,   # select cells with more than 500000 reads 
                  genesmin = 1000,       # select cells expressing more than 1k genes
                  cellmin = 3,           # select genes expressed in more than 3 cells)
                  data_status='Raw')

annotations = list(
  "blank"=c('241112', '250184', '265102', '272111', '248185', '274173'),
  "bulk"=c('225110', '251172', '273103', '280110', '235161', '246161'),
  "human"=c('233111', '249196', '257101', '264112', '233185', '247173')
)

m$excludeCellsFromIds(which(m$getCellsNames() %in% unlist(annotations)))

m$plotReadcountStats(data_status="Raw", by="timepoint", category="timepoint", basename="postQC", reads_name="read", cat_colors=unname(stage_cols))

#' Remove cells having more than 6% of mitochondrial read counts
m$removeGenesFromRatio(
  candidate_genes=grep('^MT-', m$getGeneNames(), value=T),
  threshold = 0.06
)

#' Normalize readcounts to CPM
m$normalize(method="Count-Per-Million")

#' Remove genes with an average level of 10 CPM or less, among positive cells
m$excludeUnexpressedGenes(min.cells=3, data_status='Normalized')
m$removeLowlyExpressedGenes(expression_threshold=1, selection_theshold=10, data_status='Normalized')


#' # Transcriptomic features

#' Identify modules composed of correlated genes

# "fastCor" uses the tcrossprod function which by default relies on BLAS to speed up computation. This produces inconsistent result in precision depending on the version of BLAS being used.
# see "matprod" in https://stat.ethz.ch/R-manual/R-devel/library/base/html/options.html
corr.mat = fastCor(t(m$getReadcounts(data_status='Normalized')), method="spearman")


#' ## Gene modules identification
m$identifyGeneModules(
  method="TopCorr_DR",
  corr=corr.mat,
  corr_t = 0.3,
  topcorr_mod_min_cell=0, # default
  topcorr_mod_consistency_thres=0.4, # default
  topcorr_mod_skewness_thres=-Inf, # default
  topcorr_min_cell_level=5,
  topcorr_num_max_final_gms=40,
  data_status='Normalized'
)

names(m$topCorr_DR$genemodules) <- paste0("GM ", seq(length(m$topCorr_DR$genemodules)))

# identify cell clusters
m$identifyCellClusters(method='hclust', used_genes="topCorr_DR.genemodules", data_status='Normalized')


m$plotGeneModules(
  basename='AllCells',
  displayed.gms = 'topCorr_DR.genemodules',
  displayed.geneset=NA,
  use.dendrogram='hclust',
  display.clusters=NULL,
  file_settings=list(list(type='pdf', width=10, height=10)),
  data_status='Normalized',
  gene_transformations='logscaled',
  extra_colors=cbind(
    pData(m$expressionSet)$stage_colors
  ),
  pretty.params=list("size_factor"=0.5, "ngenes_per_lines" = 8, "side.height.fraction"=.3)
)

m$writeGeneModules(basename='AllCells_allGms', gms='topCorr_DR.genemodules')

############################################
#' Remove cells with low summary readcounts
############################################
m2 = m$copy()
m2$excludeCellsFromIds(m$getReadcounts('Normalized')[unlist(m$topCorr_DR$genemodules),] %>% colSums %>% {as.numeric(scale(log(.), center=TRUE, scale=T))} %>% {which(. < -1.5)})
m2$excludeUnexpressedGenes(min.cells=1, data_status="Normalized", verbose=TRUE)

#' ## Manual feature selection

#' We select gene modules containing at least one gene known to be involved in differentiation process
bait_genes = c("HOXA2", "PAX6", "SOX2", "MSX1", "Pax3", "SALL1", "ETS1", "TWIST1", "HOMER2", "LMX1A", "VGLL2", "EYA2", "PRDM1", "FOXI3", "NELL1", "DLX5", "SOX8", "SOX10", "SOHO-1", "IRX4", "DLX6")

m2$dR$genemodules = Filter(function(x){any(bait_genes %in% x)}, m2$topCorr_DR$genemodules)


# cluster into 5 clusters 
clust.colors <- c('#da70d6', '#c71585', '#b0c4de', '#afeeee', '#5f9ea0')
#' Plot final clustering of all cells
m2$identifyCellClusters(method='hclust', clust_name="Mansel", used_genes="dR.genemodules", data_status='Normalized', numclusters=5)

gfp_counts = read.table(file=paste0(gfp_counts, 'gfpData.csv'), header=TRUE, check.names=FALSE)

m2$plotGeneModules(
  basename='AllCellsManualGMselection',
  displayed.gms = 'dR.genemodules',
  displayed.geneset=NA,
  use.dendrogram='Mansel',
  file_settings=list(list(type='pdf', width=10, height=10)),
  data_status='Normalized',
  gene_transformations=c('log', 'logscaled'),
  extra_colors=cbind(
    pData(m2$expressionSet)$stage_colors,
    m2$cellClusters$Mansel$cell_ids %>% clust.colors[.],
    "PAX2_log"=m2$getReadcounts(data_status='Normalized')['PAX2',] %>%
      {log10(1+.)} %>%
      {as.integer(1+100*./max(.))} %>%
      colorRampPalette(c("white", "black"))(n=100)[.],
    "GFP_log"= as.numeric(gfp_counts[m2$getCellsNames()]) %>%
      {log10(1+.)} %>%
      {as.integer(1+100*./max(.))} %>%
      colorRampPalette(c("white", "darkgreen"))(n=100)[.]
  ),
  pretty.params=list("size_factor"=1, "ngenes_per_lines" = 6, "side.height.fraction"=1),
  extra_legend=list("text"=names(stage_cols), "colors"=unname(stage_cols))
)

##################################################################################################################################
# Plot tsne prior to removing Pax2- cells

tsne_path = paste0(plot_path, 'allcells.tsne/') 
dir.create(tsne_path)

# here you can assign cluster colours for the tsne >> change this so that colours are directly selected by cluster number
tsne_plot(m2, m2$dR$genemodules, "allcells_clusters", seed=seed,
          cols=clust.colors[m2$cellClusters$Mansel$cell_ids], perplexity=perp, eta=eta, plot_folder = tsne_path)

tsne_plot(m2, m2$dR$genemodules, "allcells_stage", seed=seed,
          cols=pData(m2$expressionSet)$stage_colors, perplexity=perp, eta=eta, plot_folder = tsne_path)


########################################################################################################################
# Plot expression of genes of interest on tsne before filtering

# plot tsne for gradient expression of select genes in gene_list
gene_list = c('SOX2', 'SOX10', 'SOX8', 'PAX7', 'PAX2', 'LMX1A', 'SOX21', 'SIX1')
for(gn in gene_list){
  path = paste0(tsne_path, gn)
  tsne_plot(m2, m2$dR$genemodules,basename = paste0("allcells.", gn), seed=seed,
            cols=colorRampPalette(c("grey", "darkmagenta"))(n=101)[as.integer(1+100*log10(1+m2$getReadcounts(data_status='Normalized')[gn,]) / max(log10(1+m2$getReadcounts(data_status='Normalized')[gn,])))],
            perplexity=perp, pca=FALSE, eta=eta, plot_folder = tsne_path, main = gn)
}

###############################################################
# DOTPLOTS

# gene list for dotplot
gene_list = c("DLX6", "HOMER2", "FOXI3", "TFAP2E", "ZNF385C", "SIX1", "PAX2", "DLX3", "NELL1", "FGF8", "VGLL2", "EYA1", "SOHO-1", "LMX1A", "SOX8", "ZBTB16", "DLX5", "TFAP2A", # Placodes
              "SOX10", "WNT1", "MSX2", "BMP5", "PAX7", "TFAP2B", "LMO4", "ETS1", "MSX1", "SOX9", # NC
              "ZEB2", "HOXA2", "SOX2", "RFX4", "PAX6", "WNT4", "SOX21", # Neural
              "SIM1", "PITX2", "TWIST1") # Mesoderm

# get cell branch information for dotplot
cell_cluster_data = data.frame(cluster = m2$cellClusters$Mansel$cell_ids) %>%
  tibble::rownames_to_column('cellname') %>%
  dplyr::mutate(celltype = case_when(
    cluster == "1" ~ "OEP",
    cluster == "2" ~ "Late Placodal",
    cluster == "3" ~ 'Neural',
    cluster == "4" ~ 'Mesodermal',
    cluster == "5" ~ 'Neural Crest'
  ))

# gather data for dotplot
dotplot_data <- data.frame(t(m2$getReadcounts('Normalized')[gene_list, ]), check.names=F) %>%
  tibble::rownames_to_column('cellname') %>% 
  tidyr::gather(genename, value, -cellname) %>%
  dplyr::left_join(cell_cluster_data, by="cellname") %>%
  dplyr::group_by(genename, celltype) %>%
  # calculate percentage of cells in each cluster expressing gene
  dplyr::mutate('Proportion of Cells Expressing' = sum(value > 0)/n()) %>%
  # scale data
  dplyr::group_by(genename) %>%
  dplyr::mutate(value = (value - mean(value, na.rm=TRUE)) / sd(value, na.rm=TRUE)) %>%  
  # calculate mean expression
  dplyr::group_by(genename, celltype) %>%
  dplyr::mutate('Scaled Average Expression'=mean(value)) %>%
  dplyr::distinct(genename, celltype, .keep_all=TRUE) %>%
  dplyr::ungroup() %>%
  # make factor levels to order genes in dotplot
  dplyr::mutate(genename = factor(genename, levels = gene_list)) %>%
  # make factor levels to order cells in dotplott
  dplyr::mutate(celltype = factor(celltype, levels = rev(c("OEP", "Late Placodal", "Neural Crest", "Neural", "Mesodermal"))))

png(paste0(plot_path, "all_cells_dotplot.png"), width = 28, height = 12, units = "cm", res = 200)
ggplot(dotplot_data, aes(x=genename, y=celltype, size=`Proportion of Cells Expressing`, color=`Scaled Average Expression`)) +
  geom_count() +
  scale_size_area(max_size=5) +
  scale_x_discrete(position = "top") + xlab("") + ylab("") +
  scale_color_gradient(low = "grey90", high = "blue") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 0, size=7))
graphics.off()


# save rds file of all cells
saveRDS(m2, paste0(rds_path, 'm2.rds'))

########################################################################################################################
#' ## OEP derivative isolation

#' Blue cell cluster is composed of non-oep derived populations (it is also mostly comprising PAX2 negative)
#' We exclude these cells from the analysis (cluster ids, red: 1, blue: 2, green: 3, purple: 4...)
m_oep = m2$copy()
m_oep$excludeCellFromClusterIds(cluster_ids=c(3:5), used_clusters='Mansel', data_status='Normalized')

#' Some genes may not be expressed any more in the remaining cells
m_oep$excludeUnexpressedGenes(min.cells=1, data_status="Normalized", verbose=TRUE)
m_oep$removeLowlyExpressedGenes(expression_threshold=1, selection_theshold=10, data_status='Normalized')


# "fastCor" uses the tcrossprod function which by default relies on BLAS to speed up computation. This produces inconsistent result in precision depending on the version of BLAS being used.
# see "matprod" in https://stat.ethz.ch/R-manual/R-devel/library/base/html/options.html
corr.mat2 = fastCor(t(m_oep$getReadcounts(data_status='Normalized')), method="spearman")

m_oep$identifyGeneModules(
  method="TopCorr_DR",
  corr=corr.mat2,
  corr_t = 0.3,
  topcorr_mod_min_cell=0, # default
  topcorr_mod_consistency_thres=0.4, # default
  topcorr_mod_skewness_thres=-Inf, # default
  topcorr_min_cell_level=5,
  data_status='Normalized'
)

names(m_oep$topCorr_DR$genemodules) <- paste0("GM ", seq(length(m_oep$topCorr_DR$genemodules)))

m_oep$writeGeneModules(basename='OEP_allGms', gms='topCorr_DR.genemodules')

m_oep$identifyCellClusters(method='hclust', used_genes="topCorr_DR.genemodules", data_status='Normalized')

m_oep$plotGeneModules(
  basename='OEPs',
  displayed.gms = 'topCorr_DR.genemodules',
  displayed.geneset=NA,
  use.dendrogram='hclust',
  display.clusters=NULL,
  file_settings=list(list(type='pdf', width=10, height=10)),
  data_status='Normalized',
  gene_transformations='logscaled',
  extra_colors=cbind(
    pData(m_oep$expressionSet)$stage_colors
  ),
  pretty.params=list("size_factor"=1, "ngenes_per_lines" = 8, "side.height.fraction"=.3)
)

# add gene modules txt


#' Manual feature selection
bait_genes = c("HOMER2", "LMX1A", "SOHO-1", "SOX10", "VGLL2", "FOXI3", 'ZNF385C', 'NELL1', "CXCL14", "EYA4")

m_oep$topCorr_DR$genemodules.selected = Filter(function(x){any(bait_genes %in% x)}, m_oep$topCorr_DR$genemodules)

# cluster into 5 clusters 
m_oep$identifyCellClusters(method='hclust', clust_name="Mansel", used_genes="topCorr_DR.genemodules.selected", data_status='Normalized', numclusters=5)

clust.colors <- c('#ffa07a', '#f55f20', '#dda0dd', '#48d1cc', '#b2ffe5')

m_oep$plotGeneModules(
  basename='OEP_GMselection',
  displayed.gms = 'topCorr_DR.genemodules.selected',
  displayed.geneset=NA,
  use.dendrogram='Mansel',
  file_settings=list(list(type='pdf', width=10, height=5)),
  data_status='Normalized',
  gene_transformations=c('log', 'logscaled'),
  extra_colors=cbind(
    pData(m_oep$expressionSet)$stage_colors,
    m_oep$cellClusters$Mansel$cell_ids %>% clust.colors[.]
  ),
  pretty.params=list("size_factor"=1, "ngenes_per_lines" = 6, "side.height.fraction"=1),
  extra_legend=list("text"=names(stage_cols), "colors"=unname(stage_cols))
)



# save rds file of oep cells
saveRDS(m_oep, paste0(rds_path, 'm_oep.rds'))

########################################################################################################################
# Plot tSNE for oep data
tsne_path = paste0(plot_path, 'OEP_subset_tsne/')
dir.create(tsne_path)


tsne_plot(m_oep, m_oep$topCorr_DR$genemodules.selected, plot_folder = tsne_path, basename = "OEP_Clusters", seed=seed,
          cols=clust.colors[m_oep$cellClusters[['Mansel']]$cell_ids], perplexity=perp, pca=FALSE, eta=eta)

tsne_plot(m_oep, m_oep$topCorr_DR$genemodules.selected, plot_folder = tsne_path, basename = "OEP_samples", seed=1,
          cols=pData(m_oep$expressionSet)$stage_colors, perplexity=perp, pca=FALSE, eta=eta)

# Plot expression of genes of interest on OEP tsne
gene_list = c('PAX2', 'LMX1A', 'SOX8', 'TFAP2A', 'FOXI3', 'SIX1', 'ZBTB16', 'SOHO-1', 'FOXG1', 'NELL1', 'PDLIM1', 'VGLL2')

for(gn in gene_list){
  path = paste0(tsne_path, gn)
  tsne_plot(m_oep, m_oep$topCorr_DR$genemodules.selected, basename = paste0("OEP.subset.", gn), seed=seed,
            cols=colorRampPalette(c("grey", "darkmagenta"))(n=101)[as.integer(1+100*log10(1+m2$getReadcounts(data_status='Normalized')[gn,]) / max(log10(1+m2$getReadcounts(data_status='Normalized')[gn,])))],
            perplexity=perp, pca=FALSE, eta=eta, plot_folder = tsne_path, main = gn)
}



########################################################################################################################
# Plot tSNE co-expression plots

gene_pairs <- list(c("PAX2", "LMX1A"), c("PAX2", "SOX8"), c("FOXI3", "LMX1A"))
lapply(gene_pairs, function(x) {plot_tsne_coexpression(m_oep, m_oep$topCorr_DR$genemodules.selected,
                                                       gene1 = x[1], gene2 = x[2], plot_folder = tsne_path, seed=seed, perplexity=perp, pca=FALSE, eta=eta)})


##################################################################
#' ## Monocle 2

monocle.input_dims = unlist(m_oep$topCorr_DR$genemodules.selected)

df_pheno = pData(m_oep$expressionSet)
df_feature = cbind(fData(m_oep$expressionSet), 'gene_short_name'=fData(m_oep$expressionSet)$current_gene_names) # "gene_short_name" may be required by monocle
rownames(df_feature) <- df_feature$gene_short_name

HSMM <- monocle::newCellDataSet(
  as.matrix(m_oep$getReadcounts(data_status='Normalized')),
  phenoData = new("AnnotatedDataFrame", data = df_pheno),
  featureData = new('AnnotatedDataFrame', data = df_feature),
  lowerDetectionLimit = .1,
  expressionFamily=VGAM::tobit()
)

HSMM <- estimateSizeFactors(HSMM)

HSMM <- setOrderingFilter(HSMM, monocle.input_dims)

HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')

HSMM <- orderCells(HSMM)

#' Order cells from earliest "State" (ie DDRTree branch)
HSMM <- orderCells(HSMM, root_state = which.max(table(pData(HSMM)$State, pData(HSMM)$timepoint)[, "8"]))

#' Plot cell samples over projected coordinates
curr.plot.folder = paste0(plot_path, "monocle_tsne/")
dir.create(curr.plot.folder)

pdf(paste0(curr.plot.folder, 'Monocle_DDRTree_samples.pdf'))
z_order = sample(seq(m_oep$getNumberOfCells()))
plot(t(reducedDimS(HSMM))[z_order,], col=pData(m_oep$expressionSet)$stage_colors[z_order], pch=16, main=names(d), xaxt='n', ann=FALSE, yaxt='n', asp=1, cex=1.5)
dev.off()

#' Plot cell clusters over projected coordinates
#' 
#' problem here! it's not respecting the new colors
#' col=pData(m_oep$expressionSet)$stage_colors[z_order]
#' clust.colors = c('#FFA500', '#FF7F50', '#CC99CC', '#E78AC3', '#66C2A5', '#98FB98', '#E5C494', '#B3B3B3', RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))
pdf(paste0(curr.plot.folder, 'Monocle_DDRTree_Clusters.pdf'))
plot(t(reducedDimS(HSMM)), col=clust.colors[m_oep$cellClusters[['Mansel']]$cell_ids], pch=16, main=names(d), xaxt='n', ann=FALSE, yaxt='n', asp=1)
dev.off()


########################################################################
# plot gradient gene expression on monocle embeddings

gene_list = c('PAX2', 'SOX8', 'TFAP2E', 'LMX1A', 'FOXI3')
for(gn in gene_list){
  print(gn)
  pdf(paste0(curr.plot.folder, "monocle.gradient.", gn, '.pdf'))
  plot(t(reducedDimS(HSMM)), pch=16, main=gn, xlab="", ylab="", xaxt='n', yaxt='n', asp=1,
       col=colorRampPalette(c("grey", "red"))(n=101)[as.integer(1+100*log10(1+m_oep$getReadcounts(data_status='Normalized')[gn,]) / max(log10(1+m_oep$getReadcounts(data_status='Normalized')[gn,])))],
  )
  dev.off()
}

########################################################################
# plot gradient gene co-expression on monocle embeddings

gene_pairs <- list(c("FOXI3", "PAX2"), c("FOXI3", "SOX8"), c("FOXI3", "LMX1A"), c("TFAP2E", "SOX8"), c("TFAP2E", "LMX1A"))
lapply(gene_pairs, function(x) {monocle_coexpression_plot(m_oep, m_oep$topCorr_DR$genemodules.selected, monocle_obj = HSMM, gene1 = x[1], gene2 = x[2], plot_folder = curr.plot.folder)})





#' <a href="./suppl_files/monocle_plots.zip">Download all gene pattern plots</a>
#'

#' Plot annotated trajectories

p1 = plot_cell_trajectory(HSMM, color_by = "cells_samples")
p2 = plot_cell_trajectory(HSMM, color_by = "timepoint")
p3 = plot_cell_trajectory(HSMM, color_by = "Pseudotime")

pdf(paste0(plot_path, "Monocle_DDRTree_trajectories.pdf"), width=15, height=8)
gridExtra::grid.arrange(grobs=list(p1, p2, p3), layout_matrix=matrix(seq(3), ncol=3, byrow=T))
graphics.off()

#' <a href="./suppl_files/Monocle_DDRTree_trajectories.pdf">Download PDF</a>
#' <p align="center"><img src="./suppl_files/Monocle_DDRTree_trajectories.png" width="100%"></p>
#'

#' State subplots
pdf(paste0(plot_path, 'Monocle_DDRTree_State_facet.pdf'), width=7, height=4)
plot_cell_trajectory(HSMM, color_by = "State") + facet_wrap(~State, nrow = 1)
graphics.off()




########################################################################
#' Generate a monocle projection plot for each known genes

monocle_plot_folder = paste0(plot_path, 'Monocle_plots/')
dir.create(monocle_plot_folder, showWarnings = FALSE, recursive = TRUE)

genes_sel = sort(intersect(
  getDispersedGenes(m_oep$getReadcounts('Normalized'), -1),
  getHighGenes(m_oep$getReadcounts('Normalized'), mean_threshold=5)
))

gene_level = m_oep$getReadcounts("Normalized")[genes_sel %>% .[. %in% m_oep$getGeneNames()],]
gene_level.2 = t(apply(log(.1+gene_level), 1, function(x){
  pc_95 = quantile(x, .95)
  if(pc_95==0){
    pc_95=1
  }
  x <- x/pc_95
  x[x>1] <- 1
  as.integer(cut(x, breaks=10))
}))

for(n in rownames(gene_level.2)){
  print(n)
  pdf(paste0(monocle_plot_folder, 'Monocle_DDRTree_projection_', n, '.pdf'))
  plot(t(reducedDimS(HSMM)), pch=16, main=n, xlab="", ylab="", xaxt='n', yaxt='n', asp=1,
       col=colorRampPalette(c("#0464DF", "#FFE800"))(n = 10)[gene_level.2[n,]]
  )
  dev.off()
}

system(paste0("zip -rj ", plot_path, "/monocle_plots.zip ", monocle_plot_folder))
unlink(monocle_plot_folder, recursive=TRUE, force=TRUE)




###############################################################
# DOTPLOTS

# gene list for dotplot

gene_list = c("SOX10", "SOHO-1", "SOX8", "LMX1A", "PAX2", "HOMER2", "TFAP2E", "OTX2", "FOXI3", "NELL1", "PDLIM1", "VGLL2")

# get cell branch information for dotplot
cell_branch_data = pData(HSMM)[, "State", drop=F] %>%
  tibble::rownames_to_column('cellname') %>%
  dplyr::rename(branch = State) %>%
  dplyr::mutate(celltype = case_when(
    branch == "1" ~ "otic",
    branch == "2" ~ "epibranchial",
    branch == "3" ~ 'OEP'
  ))

# gather data for dotplot
dotplot_data <- data.frame(t(m_oep$getReadcounts('Normalized')[gene_list, ]), check.names=F) %>%
  tibble::rownames_to_column('cellname') %>% 
  tidyr::gather(genename, value, -cellname) %>%
  dplyr::left_join(cell_branch_data, by="cellname") %>%
  dplyr::group_by(genename, celltype) %>%
  # calculate percentage of cells in each cluster expressing gene
  dplyr::mutate('proportion of cells expressing' = sum(value > 0)/n()) %>%
  # scale data
  dplyr::group_by(genename) %>%
  dplyr::mutate(value = (value - mean(value, na.rm=TRUE)) / sd(value, na.rm=TRUE)) %>%  
  # calculate mean expression
  dplyr::group_by(genename, celltype) %>%
  dplyr::mutate('scaled average expression' = mean(value)) %>%
  dplyr::distinct(genename, celltype, .keep_all=TRUE) %>%
  dplyr::ungroup() %>%
  # make factor levels to order plot
  dplyr::mutate(genename = factor(genename, levels = gene_list))



png(paste0(plot_path, "m_oep_dotplot.png"), width = 15, height = 8, units = "cm", res = 200)
ggplot(dotplot_data, aes(x=genename, y=celltype, size=`proportion of cells expressing`, color=`scaled average expression`)) +
  geom_count() +
  scale_size_area(max_size=5) +
  scale_x_discrete(position = "top") + xlab("") + ylab("") +
  scale_color_gradient(low = "grey90", high = "blue") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 0, size=7))
graphics.off()







###############################################################
# COEXPRESSION ANALYSIS


# plot gradient gene expression on monocle embeddings
curr.plot.folder = paste0(plot_path, 'coexpression_plots/')
dir.create(curr.plot.folder)

# list of genes used for coexpression analysis
otic = c("SOX10", "SOX8", "HOMER2", 'LMX1A')
epi = c("NELL1", "FOXI3", "PDLIM1", "TFAP2E")


# plot gradient expression for genes used for coexpression analysis
for(gn in c(otic, epi)){
  print(gn)
  pdf(paste0(curr.plot.folder, gn, '.pdf'))
  plot(t(reducedDimS(HSMM)), pch=16, main=gn, xlab="", ylab="", xaxt='n', yaxt='n', asp=1,
       col=colorRampPalette(c("grey", "red"))(n=101)[as.integer(1+100*log10(1+m_oep$getReadcounts(data_status='Normalized')[gn,]) / max(log10(1+m_oep$getReadcounts(data_status='Normalized')[gn,])))],
  )
  dev.off()
}

# generate gene pairs for coexpression
comb <- t(combn(c(otic, epi), 2))

# use data cell branch data from dotplots
coexpression_data = data.frame()
for(pair in 1:nrow(comb)){
  if(sum(as.character(comb[pair,]) %in% otic) == 2){comparison = "o-o"}else if(sum(as.character(comb[pair,]) %in% otic) == 1){comparison = "o-e"}else{comparison = "e-e"}
  newdat = ldply(unique(cell_branch_data$celltype), function(y) {
    c(sum(colSums(m_oep$getReadcounts(data_status='Normalized')[c(comb[pair,1], comb[pair,2]), cell_branch_data[cell_branch_data$celltype == y, 'cellname']] > 0) == 2)/sum(cell_branch_data$celltype == y),
      comparison, y)
  })
  newdat[,1] <- as.numeric(newdat[,1])
  coexpression_data = rbind(coexpression_data, newdat)
}

# rename dataframe columns
colnames(coexpression_data) = c("value", "comparison", "branch")

# t-test between co-expression in OEP lineage and epibranchial/otic branches
coexpression_data %>%
  group_by(comparison) %>%
  pairwise_t_test(
    value ~ branch,
    ref.group = "OEP",
    p.adjust.method = "bonferroni"
  )

# calculate mean value per group and SD for plotting bar plot
plot_dat <- coexpression_data %>%
  dplyr::group_by(branch, comparison) %>%
  dplyr::summarise(
    "proportion cells co-expressing" = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  group_split(comparison)

# rename dataframes split by comparison
names(plot_dat) <- lapply(plot_dat, function(x) as.character(unique(x[["comparison"]])))

# bar plots
oe_plot <- ggplot(plot_dat$`o-e`, aes(x=branch,y=`proportion cells co-expressing`, fill = c("black", "orange", "blue"))) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=`proportion cells co-expressing`-sd, ymax=`proportion cells co-expressing`+sd), width=.2,
                position=position_dodge(.9)) +
  geom_signif(comparisons=list(c("OEP", "otic")), annotations = "***",
              y_position = 0.83, tip_length = 0.02, vjust=0.4) +
  geom_signif(comparisons=list(c("OEP", "epibranchial")), annotations = "***",
              y_position = 0.8, tip_length = 0.02, vjust=0.4) +
  ylim(c(0, 0.85)) +
  theme_classic() +
  theme(legend.position="none") +
  ggtitle("O-E") +
  theme(plot.title = element_text(hjust = 0.5))

ee_plot <- ggplot(plot_dat$`e-e`, aes(x=branch,y=`proportion cells co-expressing`, fill = c("black", "orange", "blue"))) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=`proportion cells co-expressing`-sd, ymax=`proportion cells co-expressing`+sd), width=.2,
                position=position_dodge(.9)) +
  geom_signif(comparisons=list(c("OEP", "otic")), annotations = "***",
              y_position = 0.83, tip_length = 0.02, vjust=0.4) +
  geom_signif(comparisons=list(c("OEP", "epibranchial")), annotations = "*",
              y_position = 0.8, tip_length = 0.02, vjust=0.4) +
  ylim(c(0, 0.85)) +
  theme_classic() +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        legend.position="none") +
  ggtitle("E-E") +
  theme(plot.title = element_text(hjust = 0.5))


oo_plot <- ggplot(plot_dat$`o-o`, aes(x=branch,y=`proportion cells co-expressing`, fill = c("black", "orange", "blue"))) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=`proportion cells co-expressing`-sd, ymax=`proportion cells co-expressing`+sd), width=.2,
                position=position_dodge(.9)) +
  geom_signif(comparisons=list(c("OEP", "otic")), annotations = "ns",
              y_position = 0.83, tip_length = 0.02, vjust=0.4) +
  geom_signif(comparisons=list(c("OEP", "epibranchial")), annotations = "**",
              y_position = 0.8, tip_length = 0.02, vjust=0.4) +
  ylim(c(0, 0.85)) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),
        legend.position="none") +
  ggtitle("O-O") +
  theme(plot.title = element_text(hjust = 0.5))


png(paste0(curr.plot.folder, "coexpression_test.png"), width = 20, height = 12, units = "cm", res = 200)
plot_grid(oe_plot, ee_plot, oo_plot, align = "hv", nrow = 1, rel_widths = c(1,1,1))
graphics.off()






########################################################################
# plot multiple monocle pseudotime projections in a single plots
# in order to separate the plots change separate plots to T

pseudotime_multiplot(data = HSMM, gene_list = c("PAX2", "FOXI3", "SOX8"), separate_plots = F, out_path = curr.plot.folder,
                     basename = "pseudotime.PAX2_FOXI3_SOX8")

pseudotime_multiplot(data = HSMM, gene_list = c("PAX2", "FOXI3", "LMX1A"), separate_plots = F, out_path = curr.plot.folder,
                     basename = "pseudotime.PAX2_FOXI3_LMX1A")

pseudotime_multiplot(data = HSMM, gene_list = c("PAX2", "TFAP2E", "SOX8"), separate_plots = F, out_path = curr.plot.folder,
                     basename = "pseudotime.PAX2_TFAP2E_SOX8")




# plot gradient gene expression on monocle embeddings
mon_path = paste0(plot_path, 'monocle_grad_expression/')
dir.create(mon_path)

gene_list = c('PAX2', 'SOX8', 'TFAP2E')
for(gn in gene_list){
  print(gn)
  pdf(paste0(mon_path, gn, '.pdf'))
  plot(t(reducedDimS(HSMM)), pch=16, main=gn, xlab="", ylab="", xaxt='n', yaxt='n', asp=1,
       col=colorRampPalette(c("grey", "red"))(n=101)[as.integer(1+100*log10(1+m_oep$getReadcounts(data_status='Normalized')[gn,]) / max(log10(1+m_oep$getReadcounts(data_status='Normalized')[gn,])))],
  )
  dev.off()
}





#' Plot some genes along pseudotime

some_genes=c("HOMER2", "LMX1A", "SOHO-1", "PRDM12", "FOXI3",  "TFAP2E", "VGLL2", "PDLIM1")

# if cells are missing then check pData(HSMM) to see if the correct cells are being excluded in the state column
branch1 = my_plot_genes_in_pseudotime(HSMM[some_genes, which(pData(HSMM)$State != 2)], color_by = "timepoint", relative_expr=FALSE)
branch2 = my_plot_genes_in_pseudotime(HSMM[some_genes, which(pData(HSMM)$State != 1)], color_by = "timepoint", relative_expr=FALSE)


smooth_curves = rbind.data.frame(
  cbind(branch1, "branch"="Otic"),
  cbind(branch2, "branch"="Epibranchial")
)

smooth_curves$timepoint = factor(smooth_curves$timepoint, levels=sort(unique(smooth_curves$timepoint)))

min_expr= .1

q <- ggplot(aes(Pseudotime, expression), data = smooth_curves)
q <- q + geom_point(aes_string(color = "timepoint"), size = I(.5),
                    position = position_jitter(NULL, NULL))
q <- q + geom_line(aes(x = Pseudotime, y = expectation),
                   data = smooth_curves, size=1)
q <- q + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x)))
q <- q + facet_grid(~feature_label~branch)#, nrow = NULL,
# ncol = 2, scales = "fixed")
if (min_expr < 1) {
  q <- q + expand_limits(y = c(min_expr, 1))
}
q <- q + ylab("Absolute Expression")
q <- q + xlab("Pseudotime")
# q <- q + monocle:::monocle_theme_opts()
q <- q + scale_colour_manual(values=c("#BBBDC1", "#6B98E9", "#05080D"))# breaks=c(8.5, 11, 15))
q <- q + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

pdf(paste0(plot_path, 'Monocle_DDRTree_some_genes_along_PT.pdf'), width=7, height=7)
print(q)
graphics.off()


############################################################
# BEAM

#' Use BEAM to filter the branches (from pre-filter gene list)

genes_sel = intersect(
  getDispersedGenes(m_oep$getReadcounts('Normalized'), -1),
  getHighGenes(m_oep$getReadcounts('Normalized'), mean_threshold=5)
)

# genes_sel = m_oep$getGeneNames()

branch_point_id = 1

BEAM_res <- BEAM(HSMM[genes_sel, ], branch_point = branch_point_id, cores = m_oep$num_cores)

BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

pdf(paste0(plot_path, 'Monocle_Beam.pdf'), width=7, height=40)
beam_hm = plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res, qval < .05)),],
                                      branch_point = branch_point_id,
                                      num_clusters = 20,
                                      cores = 1,
                                      use_gene_short_name = T,
                                      show_rownames = T,
                                      return_heatmap=T,
                                      branch_colors=RColorBrewer::brewer.pal(8, "Set2")[c(4,1,6)],
                                      branch_labels=c('Otic', "Epibranchial"))
graphics.off()

# save beam score to file
write.csv(BEAM_res %>% dplyr::arrange(pval), paste0(plot_path, 'beam_scores.csv'), row.names=F)


#' BEAM plot of TFs

TF_sel <- genes_to_TFs(m_oep, genes_sel)
pdf(paste0(plot_path, 'Monocle_Beam_TFs.pdf'), width=7, height=30)
beam_hm = plot_genes_branched_heatmap(HSMM[TF_sel,],
                                      branch_point = branch_point_id,
                                      # num_clusters = 4,
                                      cluster_rows=FALSE,
                                      cores = 1,
                                      use_gene_short_name = T,
                                      show_rownames = T,
                                      return_heatmap=T,
                                      branch_colors=RColorBrewer::brewer.pal(8, "Set2")[c(4,1,6)],
                                      branch_labels=c('Epibranchial', 'Otic')
)
graphics.off()



#' BEAM plot of the original known genes
pdf(paste0(plot_path, 'Monocle_Beam_knownGenes.pdf'), width=7, height=10)
beam_hm = plot_genes_branched_heatmap(HSMM[m_oep$favorite_genes,],
                                      branch_point = branch_point_id,
                                      num_clusters = 4,
                                      cores = 1,
                                      use_gene_short_name = T,
                                      show_rownames = T,
                                      return_heatmap=T,
                                      branch_colors=RColorBrewer::brewer.pal(8, "Set2")[c(4,1,6)],
                                      branch_labels=c('Epibranchial', 'Otic'))
graphics.off()


#' BEAM plot of the selected genes

beam_sel = c("FOXI3","HOMER2","PAX2","LMX1A","ZBTB16","SOHO-1","ZNF385C","SOX8","SOX10","PDLIM1","VGLL2","TFAP2E","GBX2","OTX2","DLX5","PRDM1","PRDM12","PDLIM4","EYA1","EYA2","ETV4", "NELL1") # SIX1

pdf(paste0(plot_path, 'Monocle_Beam_selGenes.pdf'), width=7, height=5)
beam_hm = plot_genes_branched_heatmap(HSMM[beam_sel,],
                                      branch_point = branch_point_id,
                                      # num_clusters = 4,
                                      cluster_rows=FALSE,
                                      cores = 1,
                                      use_gene_short_name = T,
                                      show_rownames = T,
                                      return_heatmap=T,
                                      branch_colors=RColorBrewer::brewer.pal(8, "Set2")[c(4,1,6)],
                                      branch_labels=c('Epibranchial', 'Otic')
)
graphics.off()








#####################################################################################
######                    Velocity - read and clean loom data                  ######
#####################################################################################

curr.plot.folder = paste0(plot_path, "velocyto/")
dir.create(curr.plot.folder)

# read in loom data, with ensembl ID as rownames instead of gene name
velocyto_dat <- custom_read_loom(list.files(paste0(input_path, 'velocyto'), full.names = T))

# change cell names in velocyto dat to match antler cell names
velocyto_dat <- lapply(velocyto_dat, function(x) {
  colnames(x) <- gsub(".*:", "", colnames(x))
  colnames(x) <- gsub("\\..*", "", colnames(x))
  colnames(x) <- unname(sapply(colnames(x), function(y) ifelse(
    grepl("ss8-TSS", y),
    sapply(strsplit(y, split = "_"), function(z){paste0(z[[2]], z[[3]])}),
    strsplit(y, split = "_")[[1]][[3]])))
  x
})


#####################################################################################
# run velocity on m_oep cells

# get gene annotations from anlter object
antler.gene.names <- m_oep$expressionSet@featureData@data

# get cell names in remaining dataset and associated cluster colours
clust.colors = c('#FFA500', '#FF7F50', '#CC99CC', '#E78AC3', '#66C2A5', '#98FB98', '#E5C494', '#B3B3B3', RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))
cell.colors <- clust.colors[m_oep$cellClusters$Mansel$cell_ids]
names(cell.colors) <- names(m_oep$cellClusters$Mansel$cell_ids)

# keep only cells in m_oep
m_oep_velocyto_dat <- lapply(velocyto_dat,function(x) {
  x[,colnames(x) %in% names(cell.colors)]
})

# keep only genes in cleaned antler dataset and rename genes based on antler names
m_oep_velocyto_dat <- lapply(m_oep_velocyto_dat, function(x){
  x <- x[rownames(x) %in% fData(m_oep$expressionSet)$ensembl_gene_id,]
  rownames(x) <- fData(m_oep$expressionSet)$current_gene_names[match(rownames(x), fData(m_oep$expressionSet)$ensembl_gene_id)]
  x
})

# exonic read (spliced) expression matrix
emat <- m_oep_velocyto_dat$spliced
# intronic read (unspliced) expression matrix
nmat <- m_oep_velocyto_dat$unspliced
# spanning read (intron+exon) expression matrix
smat <- m_oep_velocyto_dat$spanning

# calculate cell velocity
rvel <- gene.relative.velocity.estimates(emat,nmat,smat=smat, kCells = 5, fit.quantile = 0.05, diagonal.quantiles = TRUE)

# get tsne embeddings for m_oep cells
tsne.embeddings = tsne_embeddings(m_oep, m_oep$topCorr_DR$genemodules.selected, seed=seed, perplexity=perp, pca=FALSE, eta=eta)

# plot cell velocity on embeddings from tsne for m_oep cells
pdf(paste0(curr.plot.folder, 'OEP.subset.velocity_inc.spanning.pdf'))
show.velocity.on.embedding.cor(tsne.embeddings, rvel, n=100, scale='sqrt', cell.colors=ac(cell.colors, alpha=0.4),
                               cex=1, arrow.scale=6, arrow.lwd=1)
dev.off()

########################################################################
# velocity on monocle plots

# plot cell velocity on monocle embeddings
pdf(paste0(curr.plot.folder, 'Monocle.OEP.subset.velocity_inc.spanning.pdf'), width = 7, height = 7)
show.velocity.on.embedding.cor(t(reducedDimS(HSMM))[z_order,], rvel, n=100, scale='sqrt', cell.colors=ac(cell.colors, alpha=0.4),
                               cex=1, arrow.scale=1, arrow.lwd=0.5)
dev.off()







# ###########

# Dotplot for Williams et al. 2019 neural crest genes

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
                   mart = ensembl,
                   useCache = FALSE)

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

# make dotplot
# gene list for dotplot
gene_list = fData(m2$expressionSet) %>% filter(ensembl_gene_id %in% NC_enriched_TFs) %>% dplyr::pull(external_gene_name)

# get cell cluster information for dotplot
NC_cells = data.frame(cluster = m2$cellClusters$Mansel$cell_ids) %>%
  tibble::rownames_to_column('cellname') %>%
  filter(cluster == '5') %>%
  dplyr::mutate(celltype = 'neural crest') %>%
  dplyr::select(-cluster)

# get cell branch information for dotplot
otic_NC_cells = pData(HSMM)[, "State", drop=F] %>%
  tibble::rownames_to_column('cellname') %>%
  filter(State == '1') %>%
  dplyr::mutate(celltype = 'otic') %>%
  dplyr::select(-State) %>%
  dplyr::bind_rows(NC_cells)


# gather data for dotplot
dotplot_data <- data.frame(t(m2$getReadcounts('Normalized')[gene_list, otic_NC_cells[['cellname']]]), check.names=F) %>%
  tibble::rownames_to_column('cellname') %>% 
  tidyr::gather(genename, value, -cellname) %>%
  dplyr::left_join(otic_NC_cells, by="cellname") %>%
  dplyr::group_by(genename, celltype) %>%
  # calculate percentage of cells in each cluster expressing gene
  dplyr::mutate('proportion of cells expressing' = sum(value > 0)/n()) %>%
  # scale data
  dplyr::group_by(genename) %>%
  dplyr::mutate(value = (value - mean(value, na.rm=TRUE)) / sd(value, na.rm=TRUE)) %>%  
  # calculate mean expression
  dplyr::group_by(genename, celltype) %>%
  dplyr::mutate('scaled average expression'=mean(value)) %>%
  dplyr::distinct(genename, celltype, .keep_all=TRUE) %>%
  dplyr::ungroup() %>%
  # make factor levels to order cells in dotplott
  dplyr::mutate(celltype = factor(celltype, levels = c("neural crest", "otic")))

# order genes for dotplot based on average expression levels in placodal population
gene_order <- dotplot_data %>%
  filter(celltype == 'otic') %>%
  dplyr::arrange(-`scaled average expression`) %>%
  dplyr::pull(genename)

# add gene order to dotplot_data
dotplot_data <- dotplot_data %>%
  dplyr::mutate(genename = factor(genename, levels = gene_order))


png(paste0(plot_path, "Williams_NC_dotplot.png"), width = 30, height = 8, units = "cm", res = 200)
ggplot(dotplot_data, aes(x=genename, y=celltype, size=`proportion of cells expressing`, color=`scaled average expression`)) +
  geom_count() +
  scale_size_area(max_size=5) +
  scale_x_discrete(position = "top") + xlab("") + ylab("") +
  scale_color_gradient(low = "grey90", high = "blue") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 0, size=7)) +
  theme(legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.4, "lines"))
graphics.off()



# ###########

# Dotplot for SOX8 differentially expressed TFs across NC and otic cells
sox8_TFs <- scan('./output/antler/sox8_DE_TFs.txt', character(), quote = "")

# remove genes not expressed in m2 object
sox8_TFs <- sox8_TFs[sox8_TFs %in% rownames(m2$getReadcounts('Normalized'))]

# get cell cluster information for dotplot
NC_cells = data.frame(cluster = m2$cellClusters$Mansel$cell_ids) %>%
  tibble::rownames_to_column('cellname') %>%
  filter(cluster == '5') %>%
  dplyr::mutate(celltype = 'neural crest') %>%
  dplyr::select(-cluster)

# get cell branch information for dotplot
otic_NC_cells = pData(HSMM)[, "State", drop=F] %>%
  tibble::rownames_to_column('cellname') %>%
  filter(State == '1') %>%
  dplyr::mutate(celltype = 'otic') %>%
  dplyr::select(-State) %>%
  dplyr::bind_rows(NC_cells)


# gather data for dotplot
dotplot_data <- data.frame(t(m2$getReadcounts('Normalized')[sox8_TFs, otic_NC_cells[['cellname']]]), check.names=F) %>%
  tibble::rownames_to_column('cellname') %>% 
  tidyr::gather(genename, value, -cellname) %>%
  dplyr::left_join(otic_NC_cells, by="cellname") %>%
  dplyr::group_by(genename, celltype) %>%
  # calculate percentage of cells in each cluster expressing gene
  dplyr::mutate('proportion of cells expressing' = sum(value > 0)/n()) %>%
  # scale data
  dplyr::group_by(genename) %>%
  dplyr::mutate(value = (value - mean(value, na.rm=TRUE)) / sd(value, na.rm=TRUE)) %>%  
  # calculate mean expression
  dplyr::group_by(genename, celltype) %>%
  dplyr::mutate('scaled average expression'=mean(value)) %>%
  dplyr::distinct(genename, celltype, .keep_all=TRUE) %>%
  dplyr::ungroup() %>%
  # make factor levels to order cells in dotplot
  dplyr::mutate(celltype = factor(celltype, levels = c("neural crest", "otic"))) %>%
  # make factor levels to order genes in dotplot
  dplyr::mutate(genename = factor(genename, levels = sox8_TFs))

png(paste0(plot_path, "sox8_TFs_dotplot.png"), width = 30, height = 8, units = "cm", res = 200)
ggplot(dotplot_data, aes(x=genename, y=celltype, size=`proportion of cells expressing`, color=`scaled average expression`)) +
  geom_count() +
  scale_size_area(max_size=5) +
  scale_x_discrete(position = "top") + xlab("") + ylab("") +
  scale_color_gradient(low = "grey90", high = "blue") +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 0, size=7)) +
  theme(legend.title = element_text(size = 8),
        legend.text  = element_text(size = 8),
        legend.key.size = unit(0.4, "lines"))
graphics.off()
