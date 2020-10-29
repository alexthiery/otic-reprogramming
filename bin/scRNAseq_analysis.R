custom_functions = "./bin/custom_functions/"
input_path = "./output/"
output_path = "./output/antler"
plot_path = "./output/antler/plots/"
rds_path = "./output/antler/rds_files/"
dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)


# loadload required packages
library(Antler)

# load custom functions
sapply(list.files(custom_functions, full.names = T), source)

# set pipeline params
seed=1
perp=5
eta=200

#' Stage colors
stage_cols = setNames(c("#BBBDC1", "#6B98E9", "#05080D"), c('8', '11', '15'))

#' # Load and hygienize dataset
m = Antler$new(plot_folder=plot_path, num_cores=6)

# load in phenoData and assayData from ../dataset -> assayData is count matrix; phenoData is metaData (i.e. replicated, conditions, samples etc)
m$loadDataset(folderpath= paste0(input_path, "merged_counts/"))

pData(m$expressionSet)$stage_colors = stage_cols[as.character(pData(m$expressionSet)$timepoint)]

m$plotReadcountStats(data_status="Raw", by="timepoint", category="timepoint", basename="preQC", reads_name="read", cat_colors=unname(stage_cols))

# assigning ensembl names
m$setCurrentGeneNames(geneID_mapping_file=system.file("extdata", "Annotations/biomart_ensemblid_genename_ggallus.csv", package="Antler"))

#' Store known genes 
apriori_genes = c(
  'DACH1', 'DLX3', 'DLX5', 'DLX6', 'EYA1', 'EYA2', 'FOXG1', 'FOXI3', 'GATA3', 'GBX2', 'HESX1', 'IRX1', 'IRX2', 'IRX3', 'LMX1A', 'LMX1B', 'Pax-2', 'SALL4', 'SIX4', 'SOHO1', 'SOX10', 'SOX2', 'SOX8', 'SOX9', 'Six1', 'TBX2', # Known_otic_genes
  'ATN1', 'BACH2', 'CNOT4', 'CXCL14', 'DACH2', 'ETS1', 'FEZ1', 'FOXP1', 'HIPK1', 'HOMER2', 'IRX4', 'IRX5', 'KLF7', 'KREMEN1', 'LDB1', 'MSI1', 'PDLIM4', 'PLAG1', 'PNOC', 'PREP2', 'RERE', 'SMOC1', 'SOX13', 'TCF7L2', 'TEAD3', 'ZBTB16', 'ZFHX3', 'ZNF384', 'ZNF385C', # New_otic_TFs
  'BLIMP1', 'FOXI1', 'NKX2-6', 'NR2F2', 'PDLIM1', 'PHOX2B', 'SALL1', 'TBX10', 'TFAP2E', 'TLX1', 'VGLL2', # Epibranchial_genes
  'ZNF423', 'CXCR4', 'MAFB', 'MYC', 'Sip1', # Neural_Genes
  'CD151', 'ETS2', 'FGFR4', 'OTX2', 'PAX3', 'PAX6', 'Pax-7', 'SIX3', # Non_otic_placode_genes
  'FOXD3', 'ID2', 'ID4', 'MSX1', 'TFAP2A', 'TFAP2B', 'TFAP2C', # Neural_Crest_Genes
  'GATA2', # Epidermis_genes
  'COL11A1', 'DTX4', 'GRHL2', 'NELL1', 'OTOG', 'OTOL1', # Disease_associated_genes
  'ARID3A', 'BMP4', 'CREBBP', 'ETV4', 'ETV5', 'EYA4', 'FOXP4', 'FSTL4', 'HOXA2', 'JAG1', 'LFNG', 'LZTS1', 'MAFA', 'MEIS1', 'MYB', 'MYCN', 'NFKB1', 'NOTCH1', 'SPRY1', 'SPRY2', 'SSTR5', # chen et al. 2017
  'TWIST1', 'ENSGALG00000001876', 'ASL1', 'ENSGALG00000002558', 'ENSGALG00000011695', 'ENSGALG00000012644', 'ENSGALG00000015112', 'MAF' # other_genes
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

# save correlation matrix
saveRDS(corr.mat, paste0(rds_path, 'corr.mat.rds'))
# load cor mat if needed
# corr.mat = readRDS(paste0(rds_path, 'corr.mat.rds'))

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
    pData(m$expressionSet)$stage_colors,
    m$getReadcounts(data_status='Normalized')['Pax-2',] %>%
      ifelse(.==0, NA, .) %>% {cut(log10(1+.), breaks=100)} %>%
      colorRampPalette(c("white", "black"))(n=100)[.] %>% ifelse(is.na(.), 'red', .) %>%
      {matrix(rep(., 2), ncol=2, dimnames=list(list(), list('Pax-2', 'Red: Null')))},
    "Poorly characterized"=m$getReadcounts('Normalized')[unlist(m$topCorr_DR$genemodules),] %>% colSums %>% {as.numeric(scale(log(.), center=TRUE, scale=T))} %>% {ifelse(. < -1.5, "black", "white")}
  ),
  pretty.params=list("size_factor"=0.5, "ngenes_per_lines" = 8, "side.height.fraction"=.3)
)

m$writeGeneModules(basename='AllCells_allGms', gms='topCorr_DR.genemodules')

saveRDS(m, paste0(rds_path, 'm.rds'))
# m = readRDS(paste0(rds_path, 'm.rds'))

############################################
#' Remove cells with low summary readcounts
############################################
m2 = m$copy()
m2$excludeCellsFromIds(m$getReadcounts('Normalized')[unlist(m$topCorr_DR$genemodules),] %>% colSums %>% {as.numeric(scale(log(.), center=TRUE, scale=T))} %>% {which(. < -1.5)})
m2$excludeUnexpressedGenes(min.cells=1, data_status="Normalized", verbose=TRUE)

#' ## Manual feature selection

#' We select gene modules containing at least one gene known to be involved in differentiation process

bait_genes = c("HOXA2", "PAX6", "SOX2", "MSX1", "PAX3", "SALL1", "ETS1", "TWIST1", "HOMER2", "LMX1A", "VGLL2", "EYA2", "BLIMP1", "FOXI3", "NELL1", "DLX5", "SOX8", "SOX10", "SOHO1", "IRX4", "DLX6") #, "CXCR4")

m2$dR$genemodules = Filter(function(x){any(bait_genes %in% x)}, m2$topCorr_DR$genemodules)

#' Plot final clustering of all cells
m2$identifyCellClusters(method='hclust', clust_name="Mansel", used_genes="dR.genemodules", data_status='Normalized', numclusters=2)

gfp_counts = read.table(file=paste0(input_path, 'merged_counts/gfpData.csv'), header=TRUE, check.names=FALSE)

m2$plotGeneModules(
  basename='AllCellsManualGMselection',
  displayed.gms = 'dR.genemodules',
  displayed.geneset=NA,
  use.dendrogram='Mansel',
  display.clusters='Mansel',
  file_settings=list(list(type='pdf', width=10, height=10)),
  data_status='Normalized',
  gene_transformations=c('log', 'logscaled'),
  extra_colors=cbind(
    pData(m2$expressionSet)$stage_colors,
    "Pax-2_log"=m2$getReadcounts(data_status='Normalized')['Pax-2',] %>%
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


#' GFP distribution per timepoint X cluster id

gfp_time.df = cbind(
  "GFP"=unlist(gfp_counts[m2$getCellsNames()]),
  "clust"=factor(m2$cellClusters[['Mansel']]$cell_ids, levels=c(1,2)),
  pData(m2$expressionSet)[, c('timepoint','replicate_id')]
)
gfp_time.df$timepoint = factor(gfp_time.df$timepoint, levels=sort(unique(gfp_time.df$timepoint)))

p = gfp_time.df %>%
  ggplot() +
  geom_violin(aes(x=timepoint, y=log10(1+GFP), fill=clust, colour=clust)) + 
  scale_fill_manual(values=getClusterColors(v=2)[1:2], aesthetics = "fill") + 
  scale_colour_manual(values=getClusterColors(v=2)[1:2], aesthetics = "colour")

pdf(paste0(plot_path, 'GFP_stat.pdf'), width=4, height=4)
print(p)
graphics.off()

saveRDS(m2, paste0(rds_path, 'm2.rds'))
# m2 = readRDS(paste0(rds_path, 'm2.rds'))

##################################################################################################################################
# Plot tsne prior to removing Pax2- cells

tsne_path = paste0(plot_path, 'allcells.tsne/') 
dir.create(tsne_path)

# here you can assign cluster colours for the tsne >> change this so that colours are directly selected by cluster number
clust.colors = getClusterColors(v=2)[1:2]
tsne_plot(m2, m2$dR$genemodules, "allcells_clusters", seed=seed,
              cols=clust.colors[m2$cellClusters$Mansel$cell_ids], perplexity=perp, eta=eta, plot_folder = tsne_path)

tsne_plot(m2, m2$dR$genemodules, "allcells_stage", seed=seed,
              cols=pData(m2$expressionSet)$stage_colors, perplexity=perp, eta=eta, plot_folder = tsne_path)


########################################################################################################################
# Plot expression of Pax-2 Pax7 and Sox21 on tsne before filtering

# plot tsne for gradient expression of select genes in gene_list


gene_list = c('SOX2', 'SOX10', 'SOX8', 'Pax-7', 'Pax-2', 'LMX1A', 'SOX21', 'Six1')
for(gn in gene_list){
  path = paste0(tsne_path, gn)
  tsne_plot(m2, m2$dR$genemodules,basename = paste0("allcells.", gn), seed=seed,
                cols=colorRampPalette(c("grey", "darkmagenta"))(n=101)[as.integer(1+100*log10(1+m2$getReadcounts(data_status='Normalized')[gn,]) / max(log10(1+m2$getReadcounts(data_status='Normalized')[gn,])))],
                perplexity=perp, pca=pca, eta=eta, plot_folder = tsne_path, main = gn)
}

########################################################################################################################
#' ## OEP derivative isolation

#' Blue cell cluster is composed of non-oep derived populations (it is also mostly comprising Pax-2 negative)
#' We exclude these cells from the analysis (cluster ids, red: 1, blue: 2, green: 3, purple: 4...)

# m2 <- readRDS('./output/antler/rds_files/m2.rds')
m_oep = m2$copy()
m_oep$excludeCellFromClusterIds(cluster_ids=c(2), used_clusters='Mansel', data_status='Normalized')

#' Some genes may not be expressed any more in the remaining cells

m_oep$excludeUnexpressedGenes(min.cells=1, data_status="Normalized", verbose=TRUE)
m_oep$removeLowlyExpressedGenes(expression_threshold=1, selection_theshold=10, data_status='Normalized')


# "fastCor" uses the tcrossprod function which by default relies on BLAS to speed up computation. This produces inconsistent result in precision depending on the version of BLAS being used.
# see "matprod" in https://stat.ethz.ch/R-manual/R-devel/library/base/html/options.html
corr.mat2 = fastCor(t(m_oep$getReadcounts(data_status='Normalized')), method="spearman")


#saveRDS(corr.mat2, paste0(rds_path, 'corr.mat2.rds'))
#Load RDS corr.mat2:
# corr.mat2 = readRDS(paste0(rds_path, 'corr.mat2.rds'))

m_oep$identifyGeneModules(
  method="TopCorr_DR",
  corr=corr.mat2,
  corr_t = 0.3,
  topcorr_mod_min_cell=0, # default
  topcorr_mod_consistency_thres=0.4, # default
  topcorr_mod_skewness_thres=-Inf, # default
  topcorr_min_cell_level=5,
  topcorr_num_max_final_gms=100,
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
    pData(m_oep$expressionSet)$stage_colors,
    m_oep$getReadcounts(data_status='Normalized')['Pax-2',] %>%
      ifelse(.==0, NA, .) %>% {cut(log10(1+.), breaks=100)} %>%
      colorRampPalette(c("white", "black"))(n=100)[.] %>% ifelse(is.na(.), 'red', .) %>%
      {matrix(rep(., 2), ncol=2, dimnames=list(list(), list('Pax-2', 'Red: Null')))},
    "Poorly characterized"=m_oep$getReadcounts('Normalized')[unlist(m_oep$topCorr_DR$genemodules),] %>% colSums %>% {as.numeric(scale(log(.), center=TRUE, scale=T))} %>% {ifelse(. < -1.5, "black", "white")}
  ),
  pretty.params=list("size_factor"=1, "ngenes_per_lines" = 8, "side.height.fraction"=.3)
)

#' Manual feature selection
bait_genes = c("HOMER2", "LMX1A", "SOHO1", "SOX10", "VGLL2", "FOXI3", "ZNF385C", "NELL1")
# bait_genes = c("HOMER2", "LMX1A", "SOHO1", "SOX10", "VGLL2", "FOXI3", 'ZNF385C', 'NELL1', "CXCL14", "EYA4")

m_oep$topCorr_DR$genemodules.selected = Filter(function(x){any(bait_genes %in% x)}, m_oep$topCorr_DR$genemodules)

m_oep$identifyCellClusters(method='hclust', clust_name="Mansel", used_genes="topCorr_DR.genemodules.selected", data_status='Normalized', numclusters=5)

clust.colors = c('#FFA500', '#FF7F50', '#CC99CC', '#E78AC3', '#66C2A5', '#98FB98', '#E5C494', '#B3B3B3', RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))

m_oep$plotGeneModules(
  basename='OEP_GMselection',
  displayed.gms = 'topCorr_DR.genemodules.selected',
  displayed.geneset=NA,
  use.dendrogram='Mansel',
  display.clusters='Mansel',
  file_settings=list(list(type='pdf', width=10, height=5)),
  data_status='Normalized',
  gene_transformations=c('log', 'logscaled'),
  extra_colors=cbind(
    pData(m_oep$expressionSet)$stage_colors,
    "Pax-2_log"=m_oep$getReadcounts(data_status='Normalized')['Pax-2',] %>%
      {log10(1+.)} %>%
      {as.integer(1+100*./max(.))} %>%
      colorRampPalette(c("white", "black"))(n=100)[.],
    "GFP_log"= as.numeric(gfp_counts[m_oep$getCellsNames()]) %>%
      {log10(1+.)} %>%
      {as.integer(1+100*./max(.))} %>%
      colorRampPalette(c("white", "darkgreen"))(n=100)[.]
  ),
  pretty.params=list("size_factor"=1, "ngenes_per_lines" = 6, "side.height.fraction"=1),
  extra_legend=list("text"=names(stage_cols), "colors"=unname(stage_cols))
)

saveRDS(m_oep, paste0(rds_path, 'm_oep.rds'))
# m_oep <- readRDS(paste0(rds_path, 'm_oep.rds'))



########################################################################################################################
# Plot tSNE for oep data

tsne_path = paste0(plot_path, 'OEP_subset_tsne/')
dir.create(tsne_path)


tsne_plot(m_oep, m_oep$topCorr_DR$genemodules.selected, plot_folder = tsne_path, basename = "OEP_Clusters", seed=seed,
              cols=clust.colors[m_oep$cellClusters[['Mansel']]$cell_ids], perplexity=perp, pca=pca, eta=eta)

tsne_plot(m_oep, m_oep$topCorr_DR$genemodules.selected, plot_folder = tsne_path, basename = "OEP_samples", seed=1,
              cols=pData(m_oep$expressionSet)$stage_colors, perplexity=perp, pca=pca, eta=eta)

# Plot gradient expression of Pax-2 on OEP tsne
tsne_plot(m_oep, m_oep$topCorr_DR$genemodules.selected, basename = "OEP.subset.Pax-2", seed=seed,
              cols=colorRampPalette(c("grey", "darkmagenta"))(n=101)[as.integer(1+100*log10(1+m_oep$getReadcounts(data_status='Normalized')['Pax-2',]) / max(log10(1+m_oep$getReadcounts(data_status='Normalized')['Pax-2',])))],
              perplexity=perp, pca=pca, eta=eta, plot_folder = tsne_path, main = 'Pax-2')


########################################################################################################################
# Plot tSNE co-expression plots

gene_pairs <- list(c("Pax-2", "LMX1A"), c("Pax-2", "SOX8"), c("FOXI3", "LMX1A"))
lapply(gene_pairs, function(x) {plot_tsne_coexpression(m_oep, m_oep, m_oep$topCorr_DR$genemodules.selected,
              gene1 = x[1], gene2 = x[2], plot_folder = tsne_path, seed=seed, perplexity=perp, pca=pca, eta=eta)})


##################################################################
#' ## Monocle 2

library(monocle)

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

gene_list = c('Pax-2', 'SOX8', 'TFAP2E', 'LMX1A', 'FOXI3')
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

gene_pairs <- list(c("FOXI3", "Pax-2"), c("FOXI3", "SOX8"), c("FOXI3", "LMX1A"), c("TFAP2E", "SOX8"), c("TFAP2E", "LMX1A"))
lapply(gene_pairs, function(x) {monocle_coexpression_plot(m_oep, m_oep$topCorr_DR$genemodules.selected, monocle_obj = HSMM, gene1 = x[1], gene2 = x[2], plot_folder = curr.plot.folder)})

#####################################################################################
######                    Velocity - read and clean loom data                  ######
#####################################################################################

library(velocyto.R)
library(stringr)

velocyto_plot_path = "./output/antler/plots/velocyto"
dir.create(velocyto_plot_path)

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
# run velocity on m2 cells

# get gene annotations from anlter object
antler.gene.names <- m2$expressionSet@featureData@data

# get cell names in remaining dataset and associated cluster colours
clust.colors = getClusterColors(v=2)[1:2]
clust.colors <- clust.colors[m2$cellClusters$Mansel$cell_ids]
names(clust.colors) <- names(m2$cellClusters$Mansel$cell_ids)

# keep only cells in m2
m2_velocyto_dat <- lapply(velocyto_dat,function(x) {
  x[,colnames(x) %in% names(clust.colors)]
})

# keep only genes in cleaned antler dataset and rename genes based on antler names
m2_velocyto_dat <- lapply(m2_velocyto_dat, function(x){
  x <- x[rownames(x) %in% fData(m2$expressionSet)$ensembl_gene_id,]
  rownames(x) <- fData(m2$expressionSet)$current_gene_names[match(rownames(x), fData(m2$expressionSet)$ensembl_gene_id)]
  x
})

# exonic read (spliced) expression matrix
emat <- m2_velocyto_dat$spliced
# intronic read (unspliced) expression matrix
nmat <- m2_velocyto_dat$unspliced

# spanning read (intron+exon) expression matrix
smat <- m2_velocyto_dat$spanning

# calculate cell velocity
rvel <- gene.relative.velocity.estimates(emat,nmat,smat=smat, kCells = 5, fit.quantile = 0.05, diagonal.quantiles = TRUE)

# get tsne embeddings for all(m2) cells
tsne.embeddings = tsne_embeddings(m_oep, m_oep$topCorr_DR$genemodules.selected, seed=seed, perplexity=perp, pca=FALSE, eta=eta)

# plot cell velocity on embeddings from tsne for all(m2) cells
pdf(paste0(velocyto_plot_path, 'allcells.velocity_inc.spanning.pdf'))
show.velocity.on.embedding.cor(tsne.embeddings, rvel, n=100, scale='sqrt', cell.colors=ac(cell.colors, alpha=0.4),
                               cex=1, arrow.scale=6, arrow.lwd=1)
dev.off()



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
pdf(paste0(tsne_path, 'OEP.subset.velocity_inc.spanning.pdf'))
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











########################################################################
# plot multiple monocle pseudotime projections in a single plots
# in order to separate the plots change separate plots to T

pseudotime_multiplot(data = HSMM, gene_list = c("Pax-2", "FOXI3", "SOX8"), separate_plots = F, out_path = curr.plot.folder,
                     basename = "pseudotime.Pax-2_FOXI3_SOX8")

pseudotime_multiplot(data = HSMM, gene_list = c("Pax-2", "FOXI3", "LMX1A"), separate_plots = F, out_path = curr.plot.folder,
                     basename = "pseudotime.Pax-2_FOXI3_LMX1A")

pseudotime_multiplot(data = HSMM, gene_list = c("Pax-2", "TFAP2E", "SOX8"), separate_plots = F, out_path = curr.plot.folder,
                     basename = "pseudotime.Pax-2_TFAP2E_SOX8")


########################################################################


#' <a href="./suppl_files/Monocle_DDRTree_Clusters.pdf">Download PDF</a>
#' <p align="center"><img src="./suppl_files/Monocle_DDRTree_Clusters.png" width="80%"></p>
#'  

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





# plot gradient gene expression on monocle embeddings
mon_path = paste0(plot_path, 'monocle_grad_expression/') 
dir.create(mon_path)

gene_list = c('Pax-2', 'SOX8', 'TFAP2E')
for(gn in gene_list){
  print(gn)
  pdf(paste0(mon_path, gn, '.pdf'))
  plot(t(reducedDimS(HSMM)), pch=16, main=gn, xlab="", ylab="", xaxt='n', yaxt='n', asp=1,
       col=colorRampPalette(c("grey", "red"))(n=101)[as.integer(1+100*log10(1+m_oep$getReadcounts(data_status='Normalized')[gn,]) / max(log10(1+m_oep$getReadcounts(data_status='Normalized')[gn,])))],
  )
  dev.off()
}






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

#' <a href="./suppl_files/Monocle_DDRTree_State_facet.pdf">Download PDF</a>
#' <p align="center"><img src="./suppl_files/Monocle_DDRTree_State_facet.png" width="100%"></p>
#'  

#' Plot some genes along pseudotime

some_genes=c("HOMER2", "LMX1A", "SOHO1", "PRDM12", "FOXI3",  "TFAP2E", "VGLL2", "PDLIM1")

# if cells are missing then check pData(HSMM) to see if the correct cells are being excluded in the state column
branch2 = my_plot_genes_in_pseudotime(HSMM[some_genes, which(pData(HSMM)$State != 1)], color_by = "timepoint", relative_expr=FALSE)
branch3 = my_plot_genes_in_pseudotime(HSMM[some_genes, which(pData(HSMM)$State != 3)], color_by = "timepoint", relative_expr=FALSE)


smooth_curves = rbind.data.frame(
  cbind(branch2, "branch"="Otic"),
  cbind(branch3, "branch"="Epibranchial")
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

#' <a href="./suppl_files/Monocle_DDRTree_some_genes_along_PT.pdf">Download PDF</a>
#' <p align="center"><img src="./suppl_files/Monocle_DDRTree_some_genes_along_PT.png" width="100%"></p>
#'  

#' ## BEAM

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

pdf(paste0(plot_path, 'Monocle_Beam.pdf'), width=7, height=30)
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

#' <a href="./suppl_files/Monocle_Beam.pdf">Download PDF</a>
#' <p align="center"><img src="./suppl_files/Monocle_Beam.png" width="100%"></p>
#'  

# save beam score to file
write.csv(BEAM_res %>% dplyr::arrange(pval), paste0(plot_path, 'beam_scores.csv'), row.names=F)

#' <a href="./suppl_files/beam_scores.csv">Download BEAM scores</a>
#'  

#' BEAM plot of the selected known genes

beam_sel = c("FOXI3","HOMER2","Pax-2","LMX1A","ZBTB16","SOHO1","ZNF385C","SOX8","SOX10","PDLIM1","VGLL2","TFAP2E","GBX2","OTX2","DLX5","BLIMP1","PRDM12","PDLIM4","EYA1","EYA2","ETV4") # SIX1

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
                                      branch_labels=c('Otic', "Epibranchial")
)
graphics.off()


#' <a href="./suppl_files/Monocle_Beam_selGenes.pdf">Download PDF</a>
#' <p align="center"><img src="./suppl_files/Monocle_Beam_selGenes.png" width="100%"></p>
#'  


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
                                      branch_labels=c('Otic', "Epibranchial"))
graphics.off()


#' <a href="./suppl_files/Monocle_Beam_knownGenes.pdf">Download PDF</a>
#' <p align="center"><img src="./suppl_files/Monocle_Beam_knownGenes.png" width="100%"></p>
#'  

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
                                      branch_labels=c('Otic', "Epibranchial")
)
graphics.off()


