# install and load required packages
#renv::init()

# if renv is not installed install using
#install.packages('devtools')
library(devtools)

# changed from Antler_dev back to Antler as Antler_dev would not download
#devtools::install_github("juliendelile/Antler_dev@994c1b5b", auth_token="aafbbda529e912b5b09a6bc0a0582a4ab7e4f165")
# devtools::install_github("juliendelile/Antler")
library(Antler)

#' Extra packages are also required
# install.packages(c('Rtsne', 'gProfileR', 'tidyverse', 'stringr', 'Seurat'))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("monocle")
# devtools::install_github('hadley/multidplyr')

# velocyto.R only runs on UNIX
# BiocManager::install("pcaMethods")
# devtools::install_github("velocyto-team/velocyto.R")

# switch for use paths
user = "Alex"

if(user == "Alex"){
  # path containing functions
  my_funcs = "/Users/alex/dev/repos/OEP_APT/R_files/APT_functions/"
  
  # set path for reading in data
  data_path = '~/dev/data/OEP_APT/'
  
  #' The output path can be changed to any existing directory path
  output_path = '~/dev/output/OEP_APT/plots/'
  dir.create(output_path, recursive = T)
  
  # make new folder called rds for outputting r objects to for fast loading
  rds_path = '~/dev/output/OEP_APT/rds/'
  dir.create(rds_path)
  
} else if (user == "Ailin"){
  my_funcs = "./APT_functions/"
  data_path = '/Users/Ailin/Documents/dev_bioinf/data/'
  
  output_path = '/Users/Ailin/Documents/dev_bioinf/output/plots/'
  dir.create(output_path, recursive = T)
  rds_path = '/Users/Ailin/Documents/dev_bioinf/output/rds/'
  dir.create(rds_path)
}


###

# load APT functions
sapply(list.files(my_funcs, full.names = T), source)

#' Stage colors
stage_cols = setNames(c("#BBBDC1", "#6B98E9", "#05080D"), c('8', '11', '15'))

#' # Load and hygienize dataset
m = Antler$new(plot_folder=output_path, num_cores=4)
# devtools::load_all('.'); test=m$copy(); m=test$copy() # hack for reloading functions

# load in phenoData and assayData from ../dataset -> assayData is count matrix; phenoData is metaData (i.e. replicated, conditions, samples etc)
m$loadDataset(folderpath=paste0(data_path, 'julien_data/galgal6'))

# change timepoint values from 8.5 to 8
for(timepoint in 1:length(pData(m$expressionSet)[['timepoint']])){
  if(pData(m$expressionSet)[['timepoint']][timepoint] == 8.5){
    pData(m$expressionSet)[['timepoint']][timepoint] = 8 
  } else {
  }
}

pData(m$expressionSet)$stage_colors = stage_cols[as.character(pData(m$expressionSet)$timepoint)]


#' Display counts pre-QC
#' <p align="center"><img src="./suppl_files/preQC_read_statistics_all.png" width="80%"></p>
#' <p align="center"><img src="./suppl_files/preQC_statistics_cellNumber_replicate_id_by_timepoint.png" width="80%"></p>
#' <p align="center"><img src="./suppl_files/preQC_read_statistics_replicate_id_by_timepoint.png" width="80%"></p>
#'  

m$plotReadcountStats(data_status="Raw", by="timepoint", category="timepoint", basename="preQC", reads_name="read", cat_colors=unname(stage_cols))

###############################################
# ACCESS ANTLER DATA

# access readcount matrix from Antler object
#m$readcounts_raw

# access columnns from phenodata
#pData(m$expressionSet)[,1,drop=F]
###############################################



#' <p align="center"><img src="./suppl_files/preQC_statistics_counts_timepoint_by_timepoint.png" width="60%"></p>
#' <p align="center"><img src="./suppl_files/preQC_statistics_geneCounts_timepoint_by_timepoint.png" width="60%"></p>
#'  

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

#' Save dataset to file
saveRDS(m, file = paste0(rds_path, 'm_raw.rds'))

#' Read RDS file if needed
#' m=readRDS(paste0(rds_path, 'm_raw.rds'))

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

#' Display counts post-QC

#' <p align="center"><img src="./suppl_files/postQC_statistics_cellNumber_replicate_id_by_timepoint.png" width="60%"></p>
#' <p align="center"><img src="./suppl_files/postQC_statistics_counts_replicate_id_by_timepoint.png" width="60%"></p>
#' <p align="center"><img src="./suppl_files/postQC_statistics_geneCounts_replicate_id_by_timepoint.png" width="60%"></p>
#'  

m$plotReadcountStats(data_status="Raw", by="timepoint", category="timepoint", basename="postQC", reads_name="read", cat_colors=unname(stage_cols))

#' <p align="center"><img src="./suppl_files/postQC_statistics_counts_timepoint_by_timepoint.png" width="60%"></p>
#' <p align="center"><img src="./suppl_files/postQC_statistics_geneCounts_timepoint_by_timepoint.png" width="60%"></p>
#'  

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

#' Save dataset to file
saveRDS(m, file=paste0(rds_path, 'm_clean.rds'))

#' Read m_clean RDS file if needed
#' m=readRDS(paste0(rds_path, 'm_clean.rds'))



#' # Transcriptomic features

#' ## Gene modules identification

#' Identify modules composed of correlated genes

# "fastCor" uses the tcrossprod function which by default relies on BLAS to speed up computation. This produces inconsistent result in precision depending on the version of BLAS being used.
# see "matprod" in https://stat.ethz.ch/R-manual/R-devel/library/base/html/options.html
corr.mat = fastCor(t(m$getReadcounts(data_status='Normalized')), method="spearman")

# save correlation matrix
saveRDS(corr.mat, paste0(rds_path, 'corr.mat.rds'))

# load cor mat if needed
# corr.mat = readRDS(paste0(rds_path, 'corr.mat.rds'))

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
  displayed.geneset=m$favorite_genes,
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
  pretty.params=list("size_factor"=5, "ngenes_per_lines" = 8, "side.height.fraction"=.3)
)

#' <a href="./suppl_files/AllCells_topCorr_DR.genemodules_Normalized_logscaled.pdf">Download PDF</a>

#' <p align="center"><img src="./suppl_files/AllCells_topCorr_DR.genemodules_Normalized_logscaled.png" width="100%"></p>

m$writeGeneModules(basename='AllCells_allGms', gms='topCorr_DR.genemodules')

#' <a href="./suppl_files/AllCells_allGms_topCorr_DR.genemodules.txt">Download predicted genes</a>
#'  

#' Remove low summary readcounts cells (black cells of last side row in previous plot)
m2 = m$copy()
m2$excludeCellsFromIds(m$getReadcounts('Normalized')[unlist(m$topCorr_DR$genemodules),] %>% colSums %>% {as.numeric(scale(log(.), center=TRUE, scale=T))} %>% {which(. < -1.5)})
m2$excludeUnexpressedGenes(min.cells=1, data_status="Normalized", verbose=TRUE)

#' ## Manual feature selection

#' We select gene modules containing at least one gene known to be involved in differentiation process

bait_genes = c("HOXA2", "PAX6", "SOX2", "MSX1", "PAX3", "SALL1", "ETS1", "TWIST1", "HOMER2", "LMX1A", "VGLL2", "EYA2", "BLIMP1", "FOXI3", "NELL1", "DLX5", "SOX8", "SOX10", "SOHO1") #, "CXCR4")

m2$dR$genemodules = Filter(function(x){any(bait_genes %in% x)}, m2$topCorr_DR$genemodules)

#' Plot final clustering of all cells
m2$identifyCellClusters(method='hclust', clust_name="Mansel", used_genes="dR.genemodules", data_status='Normalized', numclusters=2)

gfp_counts = read.table(file=paste0(data_path, 'julien_data/gfp_counts.csv'), header=TRUE, check.names=FALSE)

m2$plotGeneModules(
  basename='AllCellsManualGMselection',
  displayed.gms = 'dR.genemodules',
  displayed.geneset=m2$favorite_genes,
  use.dendrogram='Mansel',
  display.clusters='Mansel',
  file_settings=list(list(type='pdf', width=10, height=10)),
  data_status='Normalized',
  gene_transformations=c('log', 'logscaled'),
  show_cells_colors=FALSE,
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
  pretty.params=list("size_factor"=5, "ngenes_per_lines" = 6, "side.height.fraction"=1),
  extra_legend=list("text"=names(stage_cols), "colors"=unname(stage_cols))
)

#' <a href="./suppl_files/AllCellsManualGMselection_dR.genemodules_Normalized_logscaled.pdf">Download PDF</a>
#' <p align="center"><img src="./suppl_files/AllCellsManualGMselection_dR.genemodules_Normalized_logscaled.png" width="100%"></p>
#'  

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

pdf(paste0(output_path, 'GFP_stat.pdf'), width=4, height=4)
print(p)
graphics.off()

# saveRDS(m2, paste0(rds_path, 'm2.rds'))
m2 = readRDS(paste0(rds_path, 'm2.rds'))

##################################################################################################################################
# Plot tsne prior to removing Pax2- cells

#check gene names
#m2$getGeneNames()[grepl("six", m2$getGeneNames(), ignore.case = T)]

tsne_path = paste0(output_path, 'allcells.tsne/') 
dir.create(tsne_path)

makeTSNEplots <- function(curr_m, gms, basename, cols=clust.colors[curr_m$cellClusters[['hclust']]$cell_ids], seed=1, pca=F, perplexity=12, eta=200, plot_folder=output_path, main = NULL){
  set.seed(seed)
  data.logscaled =  t(scale(t(log(curr_m$getReadcounts(data_status='Normalized')[unlist(gms),]+1)), center=TRUE, scale=TRUE))
  tsne_xy = Rtsne::Rtsne(t(data.logscaled), pca=pca, perplexity=perplexity, eta=eta, max_iter=2000, verbose=T)$Y
  pdf(paste0(plot_folder, '/', basename, '_TSNE.pdf'))
  plot(tsne_xy, col=cols, pch=16, main = main)
  dev.off()
}

getTSNEembeddings <- function(curr_m, gms, seed=1, pca=F, perplexity=12, eta=200){
  set.seed(seed)
  data.logscaled =  t(scale(t(log(curr_m$getReadcounts(data_status='Normalized')[unlist(gms),]+1)), center=TRUE, scale=TRUE))
  tsne_xy = Rtsne::Rtsne(t(data.logscaled), pca=pca, perplexity=perplexity, eta=eta, max_iter=2000, verbose=T)$Y
  rownames(tsne_xy) <- colnames(data.logscaled)
  return(tsne_xy)
}

perp=5
eta=200
pca=FALSE
pname="temp" # paste0("OEP_4clusters_Perp_", perp, "_eta_", eta, "_PCA_", pca)

# here you can assign cluster colours for the tsne >> change this so that colours are directly selected by cluster number
clust.colors = getClusterColors(v=2)[1:2]
makeTSNEplots(m2, m2$dR$genemodules, "allcells_clusters", seed=1,
              cols=clust.colors[m2$cellClusters$Mansel$cell_ids], perplexity=perp, pca=pca, eta=eta, plot_folder = tsne_path)

makeTSNEplots(m2, m2$dR$genemodules, "allcells_stage", seed=1,
              cols=pData(m2$expressionSet)$stage_colors, perplexity=perp, pca=pca, eta=eta, plot_folder = tsne_path)

tsne.embeddings = getTSNEembeddings(m2, m2$dR$genemodules, seed=1, perplexity=perp, pca=pca, eta=eta)

########################################################################################################################
# Plot expression of Pax-2 Pax7 and Sox21 on tsne before filtering

# plot tsne for gradient expression of select genes in gene_list


gene_list = c('SOX2', 'SOX10', 'SOX8', 'Pax-7', 'Pax-2', 'LMX1A', 'SOX21', 'Six1')
for(gn in gene_list){
  path = paste0(tsne_path, gn)
  makeTSNEplots(m2, m2$dR$genemodules,basename = paste0("allcells.", gn), seed=1,
                cols=colorRampPalette(c("grey", "darkmagenta"))(n=101)[as.integer(1+100*log10(1+m2$getReadcounts(data_status='Normalized')[gn,]) / max(log10(1+m2$getReadcounts(data_status='Normalized')[gn,])))],
                perplexity=perp, pca=pca, eta=eta, plot_folder = tsne_path, main = gn)
}

#####################################################################################
######                    Velocity - read and clean loom data                  ######
#####################################################################################

library(velocyto.R)
library(stringr)

# read in loom data, with ensembl ID as rownames instead of gene name
ss8 <- custom.read.loom("~/dev/data/OEP_APT/alignment_200423/velocity_alignment/ss8/onefilepercell_TSS_P1_A10-1234_and_others_BXX77.loom",)
ss11 <- custom.read.loom("~/dev/data/OEP_APT/alignment_200423/velocity_alignment/ss11/onefilepercell_WTCHG_524615_201106_and_others_T7BPH.loom",)
ss15 <- custom.read.loom("~/dev/data/OEP_APT/alignment_200423/velocity_alignment/ss15/onefilepercell_WTCHG_528508_201108_and_others_KBRCF.loom",)


# change column names to match antler cell names
ss8 <- lapply(ss8, function(x) {
  colnames(x) <- gsub(".*:TSS_", "", gsub("-1234\\..*", "", colnames(x)))
  colnames(x) <- ifelse(str_count(colnames(x), "_") > 1, gsub("_", "", gsub('(.*)_\\w+', '\\1', colnames(x))), gsub("_", "", colnames(x)))
  x
})

ss11 <- lapply(ss11, function(x){
  colnames(x) <- gsub("\\..*", "", gsub(".*_","", colnames(ss11$spliced)))
  x
})

ss15 <- lapply(ss15, function(x){
  colnames(x) <- gsub("\\..*", "", gsub(".*_","", colnames(ss15$spliced)))
  x
})

# merge data from different stages (cbind.rownames orders the rownames in each dataframe and then do.call(cbind))
ldat <- list(spliced = cbind.rownames(ss8$spliced, ss11$spliced, ss15$spliced),
             unspliced = cbind.rownames(ss8$unspliced, ss11$unspliced, ss15$unspliced),
             ambiguous = cbind.rownames(ss8$ambiguous, ss11$ambiguous, ss15$ambiguous),
             spanning = cbind.rownames(ss8$spanning, ss11$spanning, ss15$spanning))

#####################################################################################
# run velocity on m2 cells

# get gene annotations from anlter object
antler.gene.names <- m2$expressionSet@featureData@data

# get cell names in remaining dataset and associated cluster colours
clust.colors = getClusterColors(v=2)[1:2]
clust.colors <- clust.colors[m2$cellClusters$Mansel$cell_ids]
names(clust.colors) <- names(m2$cellClusters$Mansel$cell_ids)

# keep onluy cells in m2
ldat <- lapply(ldat,function(x) {
  x[,colnames(x) %in% names(clust.colors)]
})

# keep only genes in cleaned antler dataset and rename genes based on antler names
ldat <- lapply(ldat, function(x){
  x <- x[rownames(x) %in% fData(m2$expressionSet)$ensembl_gene_id,]
  rownames(x) <- fData(m2$expressionSet)$current_gene_names[match(rownames(x), fData(m2$expressionSet)$ensembl_gene_id)]
  x
})

# exonic read (spliced) expression matrix
emat <- ldat$spliced
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$spanning

# calculate cell velocity
rvel <- gene.relative.velocity.estimates(emat,nmat,smat=smat, kCells = 5, fit.quantile = 0.05, diagonal.quantiles = TRUE)

# plot cell velocity on embeddings from tsne for all cells

pdf(paste0(tsne_path, 'allcells.velocity_inc.spanning.pdf'))
show.velocity.on.embedding.cor(tsne.embeddings, rvel, n=100, scale='sqrt', cell.colors=ac(clust.colors, alpha=0.4),
                               cex=1, arrow.scale=4, arrow.lwd=1)
dev.off()



########################################################################################################################
#' ## OEP derivative isolation

#' Blue cell cluster is composed of non-oep derived populations (it is also mostly comprising Pax-2 negative)
#' We exclude these cells from the analysis (cluster ids, red: 1, blue: 2, green: 3, purple: 4...)


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
  displayed.geneset=m_oep$favorite_genes,
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
  pretty.params=list("size_factor"=5, "ngenes_per_lines" = 8, "side.height.fraction"=.3)
)

#' <a href="./suppl_files/OEPs_topCorr_DR.genemodules_Normalized_logscaled.pdf">Download PDF</a>
#' <p align="center"><img src="./suppl_files/OEPs_topCorr_DR.genemodules_Normalized_logscaled.png" width="100%"></p>
#'  

#' Manual feature selection

bait_genes = c("HOMER2", "LMX1A", "SOHO1", "SOX10", "VGLL2", "FOXI3") # 'ZNF385C', 'NELL1' # "CXCL14", "EYA4", "RERE", "SOX13", 

m_oep$topCorr_DR$genemodules.selected = Filter(function(x){any(bait_genes %in% x)}, m_oep$topCorr_DR$genemodules)

m_oep$identifyCellClusters(method='hclust', clust_name="Mansel", used_genes="topCorr_DR.genemodules.selected", data_status='Normalized', numclusters=6)

clust.colors = c('#FFA500', '#FF7F50', '#CC99CC', '#E78AC3', '#66C2A5', '#98FB98', '#E5C494', '#B3B3B3', RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))

#Previous cluster color code below:
#clust.colors = c('#FFD92F', '#FC8D62', '#8DA0CB', '#E78AC3', '#66C2A5', '#A6D854', '#E5C494', '#B3B3B3', RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))
#color test1:
#clust.colors = c('#FFD92F', '#FF7F50', '#CC99CC', '#E78AC3', '#66C2A5', '#98FB98', '#E5C494', '#B3B3B3', RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))
#clust colors test 2:
#clust.colors = c('#FFA500', '#FF7F50', '#DD8AE7', '#E78AC3', '#66C2A5', '#98FB98', '#E5C494', '#B3B3B3', RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))

m_oep$plotGeneModules(
  basename='OEP_GMselection',
  displayed.gms = 'topCorr_DR.genemodules.selected',
  displayed.geneset=m_oep$favorite_genes,
  use.dendrogram='Mansel',
  display.clusters='Mansel',
  file_settings=list(list(type='pdf', width=10, height=5)),
  data_status='Normalized',
  gene_transformations=c('log', 'logscaled'),
  show_cells_colors=FALSE,
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
  pretty.params=list("size_factor"=5, "ngenes_per_lines" = 6, "side.height.fraction"=1),
  cluster_colors = clust.colors,
  extra_legend=list("text"=names(stage_cols), "colors"=unname(stage_cols))
)

#' <a href="./suppl_files/OEP_GMselection_topCorr_DR.genemodules.selected_Normalized_logscaled.pdf">Download PDF</a>

#' <p align="center"><img src="./suppl_files/OEP_GMselection_topCorr_DR.genemodules.selected_Normalized_logscaled.png" width="100%"></p>
#'  

saveRDS(m_oep, paste0(rds_path, 'm_oep.rds'))
# m_oep <- readRDS(paste0(rds_path, 'm_oep.rds'))
#' # 2D Projections & Pseudotime Reconstruction

#' ## tSNE

makeTSNEplots <- function(curr_m, gms, basename, cols=clust.colors[curr_m$cellClusters[['hclust']]$cell_ids], seed=1, pca=F, perplexity=12, eta=200, plot_folder=output_path, main = NULL){
  
  set.seed(seed)
  
  # tsne_xy = Rtsne::Rtsne(t(curr_m$getReadcounts(data_status='Normalized')[unlist(gms),]), pca=F, perplexity=12, max_iter=2000, verbose=T)$Y
  # pdf(paste0(output_path, basename, '.TSNE.pdf'))
  # plot(tsne_xy, col=clust.colors[curr_m$cellClusters[['hclust']]$cell_ids], pch=16)
  # dev.off()
  
  data.logscaled =  t(scale(t(log(curr_m$getReadcounts(data_status='Normalized')[unlist(gms),]+1)), center=TRUE, scale=TRUE))
  
  tsne_xy = Rtsne::Rtsne(t(data.logscaled), pca=pca, perplexity=perplexity, eta=eta, max_iter=2000, verbose=T)$Y
  
  
  pdf(paste0(plot_folder, '/', basename, '_TSNE.pdf'))
  plot(tsne_xy, col=cols, pch=16, main = main)
  dev.off()
}

perp=5
eta=200
pca=FALSE


curr.plot.folder = paste0(output_path, 'OEP_subset_tsne/')
dir.create(curr.plot.folder)

makeTSNEplots(m_oep, m_oep$topCorr_DR$genemodules.selected, plot_folder = curr.plot.folder, basename = "OEP_Clusters", seed=1,
              cols=clust.colors[m_oep$cellClusters[['Mansel']]$cell_ids], perplexity=perp, pca=pca, eta=eta)


#' <a href="./suppl_files/OEP_Clusters_TSNE.pdf">Download PDF</a>
#' <p align="center"><img src="./suppl_files/OEP_Clusters_TSNE.png" width="80%"></p>
#'  

makeTSNEplots(m_oep, m_oep$topCorr_DR$genemodules.selected, plot_folder = curr.plot.folder, basename = "OEP_samples", seed=1,
              cols=pData(m_oep$expressionSet)$stage_colors, perplexity=perp, pca=pca, eta=eta)

tsne.embeddings = getTSNEembeddings(m_oep, m_oep$topCorr_DR$genemodules.selected, seed=1, perplexity=perp, pca=pca, eta=eta)

########################################################################################################################
# PlottSNE co-expression plots

perp=5
eta=200
pca=FALSE

TSNEcoexpression.plots(m_oep, m_oep$topCorr_DR$genemodules.selected, gene1 = "Pax-2", gene2 = "LMX1A",
                       plot_folder = "~/Desktop/", seed=1, perplexity=perp, pca=pca, eta=eta)

TSNEcoexpression.plots(m_oep, m_oep$topCorr_DR$genemodules.selected, gene1 = "Pax-2", gene2 = "SOX8",
                       plot_folder = "~/Desktop/", seed=1, perplexity=perp, pca=pca, eta=eta)

TSNEcoexpression.plots(m_oep, m_oep$topCorr_DR$genemodules.selected, gene1 = "FOXI3", gene2 = "LMX1A",
                       plot_folder = "~/Desktop/", seed=1, perplexity=perp, pca=pca, eta=eta)

########################################################################################################################
# Plot gradient expression of Pax-2 on OEP tsne

makeTSNEplots(m_oep, m_oep$topCorr_DR$genemodules.selected, basename = "OEP.subset.Pax-2", seed=1,
              cols=colorRampPalette(c("grey", "darkmagenta"))(n=101)[as.integer(1+100*log10(1+m_oep$getReadcounts(data_status='Normalized')['Pax-2',]) / max(log10(1+m_oep$getReadcounts(data_status='Normalized')['Pax-2',])))],
              perplexity=perp, pca=pca, eta=eta, plot_folder = curr.plot.folder, main = 'Pax-2')


#####################################################################################
# run velocity on m_oep cells

# get gene annotations from anlter object
antler.gene.names <- m_oep$expressionSet@featureData@data

# get cell names in remaining dataset and associated cluster colours
clust.colors = c('#FFA500', '#FF7F50', '#CC99CC', '#E78AC3', '#66C2A5', '#98FB98', '#E5C494', '#B3B3B3', RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))
cell.colors <- clust.colors[m_oep$cellClusters$Mansel$cell_ids]
names(cell.colors) <- names(m_oep$cellClusters$Mansel$cell_ids)

# keep onluy cells in m_oep
ldat <- lapply(ldat,function(x) {
  x[,colnames(x) %in% names(cell.colors)]
})

# keep only genes in cleaned antler dataset and rename genes based on antler names
ldat <- lapply(ldat, function(x){
  x <- x[rownames(x) %in% fData(m_oep$expressionSet)$ensembl_gene_id,]
  rownames(x) <- fData(m_oep$expressionSet)$current_gene_names[match(rownames(x), fData(m_oep$expressionSet)$ensembl_gene_id)]
  x
})

# exonic read (spliced) expression matrix
emat <- ldat$spliced
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$spanning

# calculate cell velocity
rvel <- gene.relative.velocity.estimates(emat,nmat,smat=smat, kCells = 5, fit.quantile = 0.05, diagonal.quantiles = TRUE)

# plot cell velocity on embeddings from tsne for all cells

pdf(paste0(curr.plot.folder, 'OEP.subset.velocity_inc.spanning.pdf'))
show.velocity.on.embedding.cor(tsne.embeddings, rvel, n=100, scale='sqrt', cell.colors=ac(cell.colors, alpha=0.4),
                               cex=1, arrow.scale=4, arrow.lwd=1)
dev.off()

########################################################################################################################

# get transcription factor subset from gene modules

TF_sub <- modules_to_TFs(Antler_obj = m_oep, gene_modules = m_oep$topCorr_DR$genemodules.selected)

##################################################################
# Plot dotplot for all cells -> to identify identity of unknown clusters

# re-cluster OEP subset into 3 clusters
m_oep$identifyCellClusters(method='hclust', clust_name="Mansel", used_genes="topCorr_DR.genemodules.selected", data_status='Normalized', numclusters=3)

# set cluster colours
clust.colors = c('#FFA500', '#FF7F50', '#CC99CC', '#E78AC3', '#66C2A5', '#98FB98', '#E5C494', '#B3B3B3', RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))

# plot gene modules for OEPs with new clustering
m_oep$plotGeneModules(
  basename='OEP_3_clust',
  displayed.gms = 'topCorr_DR.genemodules.selected',
  displayed.geneset=m_oep$favorite_genes,
  use.dendrogram='Mansel',
  display.clusters='Mansel',
  file_settings=list(list(type='pdf', width=10, height=5)),
  data_status='Normalized',
  gene_transformations=c('log', 'logscaled'),
  #  show_cells_colors=FALSE,
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
  pretty.params=list("size_factor"=5, "ngenes_per_lines" = 6, "side.height.fraction"=1),
  # cluster_colors = clust.colors,
  extra_legend=list("text"=names(stage_cols), "colors"=unname(stage_cols))
)

# reload entire Antler dataset for normalisation in Seurat
Seur_obj <- readRDS(paste0(rds_path, 'm_raw.rds'))

# initialise seurat object and carry out CPM and log normalisation
library(Seurat)
Seur_obj <- CreateSeuratObject(Seur_obj$getReadcounts(data_status='Raw'), project = "Ailin", assay = "raw")
Seur_obj <- NormalizeData(Seur_obj, normalization.method = 'LogNormalize', scale.factor = 1e6)

# filter seurat object by m_oep (ANTLER object) cells and genes
Seur_obj <- Seur_obj[rownames(m_oep$getReadcounts(data_status = 'Raw')),colnames(m_oep$getReadcounts(data_status = 'Raw'))]

# add cluster information (from ANTLER m_oep) to seurat metadata
Seur_obj@meta.data$Antler_clusters <- m_oep$cellClusters$Mansel$cell_ids[match(names(m_oep$cellClusters$Mansel$cell_ids), rownames(Seur_obj@meta.data))]

Seur_obj@meta.data$cluster_ident <- apply(Seur_obj@meta.data, 1, function(x)
  if(x["Antler_clusters"] == 1){"Epibranchial"} else if(x["Antler_clusters"] == 2){"OEP"}
  else if(x["Antler_clusters"] == 3){"Otic"} else {NA})

# select genes and plot dotplot
gene_list <- c('TFAP2E', 'FOXI3', 'OTX2', 'Pax-2', "SOX8", 'LMX1A', 'SOHO1', 'VGLL2')

pdf(paste0(output_path, 'OEP_dotplot.pdf'), height = 8, width = 10)
DotPlot(Seur_obj, features = gene_list, group.by = "cluster_ident") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
graphics.off()

##################################################################




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
curr.plot.folder = paste0(output_path, "monocle_tsne/")
dir.create(curr.plot.folder)

set.seed(1)
pdf(paste0(curr.plot.folder, 'Monocle_DDRTree_samples.pdf'))
z_order = sample(seq(m_oep$getNumberOfCells()))
plot(t(reducedDimS(HSMM))[z_order,], col=pData(m_oep$expressionSet)$stage_colors[z_order], pch=16, main=names(d), xaxt='n', ann=FALSE, yaxt='n', asp=1, cex=1.5)
# plot(t(reducedDimS(HSMM)), col="black", bg=ailin_cols[as.character(pData(m_oep$expressionSet)$timepoint)], pch=21, main=names(d), xaxt='n', ann=FALSE, yaxt='n', asp=1, cex=1.5)
dev.off()

#' <a href="./suppl_files/Monocle_DDRTree_samples.pdf">Download PDF</a>
#' <p align="center"><img src="./suppl_files/Monocle_DDRTree_samples.png" width="80%"></p>
#'  

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
Monocle.coexpression.plots(m_oep, m_oep$topCorr_DR$genemodules.selected, monocle_obj = HSMM, gene1 = "FOXI3", gene2 = "Pax-2", plot_folder = "~/Desktop/")
Monocle.coexpression.plots(m_oep, m_oep$topCorr_DR$genemodules.selected, monocle_obj = HSMM, gene1 = "FOXI3", gene2 = "SOX8", plot_folder = "~/Desktop/")
Monocle.coexpression.plots(m_oep, m_oep$topCorr_DR$genemodules.selected, monocle_obj = HSMM, gene1 = "FOXI3", gene2 = "LMX1A", plot_folder = "~/Desktop/")
Monocle.coexpression.plots(m_oep, m_oep$topCorr_DR$genemodules.selected, monocle_obj = HSMM, gene1 = "TFAP2E", gene2 = "SOX8", plot_folder = "~/Desktop/")
Monocle.coexpression.plots(m_oep, m_oep$topCorr_DR$genemodules.selected, monocle_obj = HSMM, gene1 = "TFAP2E", gene2 = "LMX1A", plot_folder = "~/Desktop/")

########################################################################
# velocity on monocle plots

# get gene annotations from anlter object
antler.gene.names <- m_oep$expressionSet@featureData@data

# merge data from different stages (cbind.rownames orders the rownames in each dataframe and then do.call(cbind))
ldat <- list(spliced = cbind.rownames(ss8$spliced, ss11$spliced, ss15$spliced),
             unspliced = cbind.rownames(ss8$unspliced, ss11$unspliced, ss15$unspliced),
             ambiguous = cbind.rownames(ss8$ambiguous, ss11$ambiguous, ss15$ambiguous),
             spanning = cbind.rownames(ss8$spanning, ss11$spanning, ss15$spanning))

# get cell names in remaining dataset and associated cluster colours
clust.colors = c('#FFD92F', '#FC8D62', '#8DA0CB', '#E78AC3', '#66C2A5', '#A6D854', '#E5C494', '#B3B3B3', RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))
cell.colors <- clust.colors[m_oep$cellClusters$Mansel$cell_ids]
names(cell.colors) <- names(m_oep$cellClusters$Mansel$cell_ids)

# keep onluy cells in m_oep
ldat <- lapply(ldat,function(x) {
  x[,colnames(x) %in% names(cell.colors)]
})

# keep only genes in cleaned antler dataset and rename genes based on antler names
ldat <- lapply(ldat, function(x){
  x <- x[rownames(x) %in% fData(m_oep$expressionSet)$ensembl_gene_id,]
  rownames(x) <- fData(m_oep$expressionSet)$current_gene_names[match(rownames(x), fData(m_oep$expressionSet)$ensembl_gene_id)]
  x
})

# exonic read (spliced) expression matrix
emat <- ldat$spliced
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$spanning

# calculate cell velocity
rvel <- gene.relative.velocity.estimates(emat,nmat,smat=smat, kCells = 5, fit.quantile = 0.05, diagonal.quantiles = TRUE)

# plot cell velocity on embeddings from tsne for all cells
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

monocle_plot_folder = paste0(output_path, 'Monocle_plots/')
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

system(paste0("zip -rj ", output_path, "/monocle_plots.zip ", monocle_plot_folder))
unlink(monocle_plot_folder, recursive=TRUE, force=TRUE)





# plot gradient gene expression on monocle embeddings
mon_path = paste0(output_path, 'monocle_grad_expression/') 
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

pdf(paste0(output_path, "Monocle_DDRTree_trajectories.pdf"), width=7, height=4)
gridExtra::grid.arrange(grobs=list(p1, p2, p3), layout_matrix=matrix(seq(3), ncol=3, byrow=T))
graphics.off()

#' <a href="./suppl_files/Monocle_DDRTree_trajectories.pdf">Download PDF</a>
#' <p align="center"><img src="./suppl_files/Monocle_DDRTree_trajectories.png" width="100%"></p>
#'  

#' State subplots
pdf(paste0(output_path, 'Monocle_DDRTree_State_facet.pdf'), width=7, height=4)
plot_cell_trajectory(HSMM, color_by = "State") + facet_wrap(~State, nrow = 1)
graphics.off()

#' <a href="./suppl_files/Monocle_DDRTree_State_facet.pdf">Download PDF</a>
#' <p align="center"><img src="./suppl_files/Monocle_DDRTree_State_facet.png" width="100%"></p>
#'  

#' Plot some genes along pseudotime

some_genes=c("HOMER2", "LMX1A", "SOHO1", "PRDM12", "FOXI3",  "TFAP2E", "VGLL2", "PDLIM1")

source('../R_files/my_plot_genes_in_pseudotime.R')
# if cells are missing then check pData(HSMM) to see if the correct cells are being excluded in the state column
branch2 = my_plot_genes_in_pseudotime(HSMM[some_genes, which(pData(HSMM)$State != 1)], color_by = "timepoint", relative_expr=FALSE)
branch3 = my_plot_genes_in_pseudotime(HSMM[some_genes, which(pData(HSMM)$State != 2)], color_by = "timepoint", relative_expr=FALSE)


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

pdf(paste0(output_path, 'Monocle_DDRTree_some_genes_along_PT.pdf'), width=7, height=7)
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

BEAM_res <- BEAM(HSMM[genes_sel, ], branch_point = branch_point_id, cores = m$num_cores)

BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

pdf(paste0(output_path, 'Monocle_Beam.pdf'), width=7, height=30)
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
write.csv(BEAM_res %>% dplyr::arrange(pval), paste0(output_path, 'beam_scores.csv'), row.names=F)

#' <a href="./suppl_files/beam_scores.csv">Download BEAM scores</a>
#'  

#' BEAM plot of the selected known genes

beam_sel = c("FOXI3","HOMER2","Pax-2","LMX1A","ZBTB16","SOHO1","ZNF385C","SOX8","SOX10","PDLIM1","VGLL2","TFAP2E","GBX2","OTX2","DLX5","BLIMP1","PRDM12","PDLIM4","EYA1","EYA2","ETV4") # SIX1

pdf(paste0(output_path, 'Monocle_Beam_selGenes.pdf'), width=7, height=5)
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
pdf(paste0(output_path, 'Monocle_Beam_knownGenes.pdf'), width=7, height=10)
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

pdf(paste0(output_path, 'Monocle_Beam_TFs.pdf'), width=7, height=30)
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





#' # Combinatorial tests

#' Identify patterned genes

library(monocle)

clust_ids = m_oep$cellClusters[['Mansel']]$cell_ids
numclusters_final = max(m_oep$cellClusters[['Mansel']]$cell_ids)

# prefilter with dispersion and mean level

genes_sel = intersect(
  getDispersedGenes(m_oep$getReadcounts('Normalized'), -1),
  getHighGenes(m_oep$getReadcounts('Normalized'), mean_threshold=3)
)

HSMM <- monocle::newCellDataSet(
  as.matrix(m_oep$getReadcounts(data_status='Normalized')),
  phenoData = new("AnnotatedDataFrame", data=cbind(pData(m_oep$expressionSet), clust_ids = clust_ids)),
  featureData = new('AnnotatedDataFrame', data = data.frame(row.names=m_oep$getGeneNames(), "gene_short_name"=m_oep$getGeneNames())),
  lowerDetectionLimit = 10,
  expressionFamily=VGAM::tobit()
)

HSMM <- estimateSizeFactors(HSMM)

test0 = monocle::differentialGeneTest(HSMM[genes_sel,],
                                      fullModelFormulaStr="~clust_ids",
                                      reducedModelFormulaStr="~1",
                                      cores=10,
                                      relative_expr=1, 
                                      verbose=T)

qval_threshold = 5e-2
patterned_genes = test0 %>% as_tibble %>% dplyr::arrange(qval) %>% dplyr::filter(qval < qval_threshold) %>% .$gene_short_name %>% as.character
# patterned_genes <- patterned_genes[1:10]

#' Generate all combinations of cluster identities

allcomb = expand.grid(rep(list(0:1), numclusters_final))
colnames(allcomb) <- as.character(seq(numclusters_final))

gene_domain_comb_tbl = as_tibble(allcomb) %>% 
  tidyr::unite(model, remove=F, sep="") %>% # add model name (concatenate all design)
  dplyr::mutate(model=paste0("m", model)) %>%
  dplyr::filter(!model %in% c(
    paste0("m", paste0(rep(1, numclusters_final), collapse="")),
    paste0("m", paste0(rep(0, numclusters_final), collapse=""))
  )) %>% # remove the no/all domain(s) models
  tidyr::crossing(tibble(genename=patterned_genes), ., by=NULL)


all_design_domain_level = gene_domain_comb_tbl %>% dplyr::select(-genename) %>% dplyr::distinct(model, .keep_all=T)

clust_cell_tbl = tibble(clust_ids=as.character(clust_ids), cellname=m_oep$getCellsNames())

all_design_cell_level = all_design_domain_level %>%
  tidyr::gather(clust_ids, part_of, -model) %>%
  dplyr::left_join(clust_cell_tbl, by=c("clust_ids"="clust_ids")) %>%
  dplyr::select(-clust_ids) %>%
  tidyr::spread(model, part_of) %>%
  dplyr::arrange(match(cellname, m_oep$getCellsNames())) %>% # order as in original dataset
  tibble::column_to_rownames("cellname") %>%
  as.data.frame()

# add models for tests
pData(HSMM) <- cbind(pData(HSMM), all_design_cell_level)

#' Run tests

model_test <- function(x){
  
  # print(x)
  test= monocle::differentialGeneTest(HSMM[unlist(x$genelist), ,drop=F],
                                      # test= monocle::differentialGeneTest(HSMM[unlist(x$genelist), ,drop=F],
                                      fullModelFormulaStr=paste0("~", x$model),
                                      reducedModelFormulaStr = "~1",
                                      cores=1,
                                      relative_expr=F, 
                                      verbose=T)
  # print(test)
  
  cbind(
    test,
    model=x$model
  )
}

cl <- multidplyr::create_cluster(10) # too many cores trigger segmentation fault (lack of memory when HSMM is copied on each cluster)
multidplyr::cluster_copy(cl, model_test)
multidplyr::cluster_copy(cl, HSMM)

# multidplyr::cluster_assign_value(cl, "my_diff_test_helper", my_diff_test_helper)
# multidplyr::cluster_assign_value(cl, "my_differentialGeneTest", my_differentialGeneTest)
# multidplyr::cluster_assign_value(cl, "my_compareModels", my_compareModels)
# multidplyr::cluster_assign_value(cl, "my_smartEsApply", my_smartEsApply)

cl %>% multidplyr::cluster_ls()

t0 = Sys.time()
Monocle_run = gene_domain_comb_tbl %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(genelist = list(genename)) %>%
  multidplyr::partition(model, cluster=cl) %>%
  dplyr::do(model_test(.)) %>%
  dplyr::collect()
t1 = Sys.time()

print(t1 - t0)

parallel::stopCluster(cl)

#' Select top model and filter genes

# pop_sign: 1 for positive, 0 for negative
# type_step1: "All", "Progenitor", or "Neuron"
mean_level_domains <- function(x, pop_sign=1, type_step1="All"){
  domain_list = seq(numclusters_final)[which(strsplit(substring(x$model, 2), split="")[[1]]==
                                               pop_sign)]
  cell_list = which(m_oep$cellClusters[['Mansel']]$cell_ids %in% domain_list)
  
  mean(m_oep$getReadcounts('Normalized')[x$genename, cell_list, drop=F])
}

ratio_expr <- function(x, pop_sign=1, type_step1="All"){
  domain_list = seq(numclusters_final)[which(strsplit(substring(x$model, 2), split="")[[1]]==
                                               pop_sign)]
  cell_list = which(m_oep$cellClusters[['Mansel']]$cell_ids %in% domain_list)
  
  mean(1*(m_oep$getReadcounts('Raw')[x$genename, cell_list, drop=F]>0))
}

#' Select top model and calculate associated metrics (fold-change, mean level in both clusters…)

Monocle_run_topModel = Monocle_run %>%
  dplyr::ungroup() %>%
  dplyr::mutate(genename=as.character(gene_short_name),
                gene_short_name=NULL) %>%
  dplyr::group_by(genename) %>%
  dplyr::filter(pval <= 1.00001 * min(pval)) %>% # we add this step before adding metrics to reduce computation time. 
  dplyr::ungroup() %>%
  # dplyr::arrange(num_domain, model) %>%
  dplyr::rowwise() %>%
  dplyr::do({res=dplyr::as_data_frame(.)      # add time correlation
  res$mean_pos = mean_level_domains(res, pop_sign=1, type_step1="All")
  res$mean_neg = mean_level_domains(res, pop_sign=0, type_step1="All")
  res$ratio_expr_pos = ratio_expr(res, pop_sign=1, type_step1="All")
  res$ratio_expr_neg = ratio_expr(res, pop_sign=0, type_step1="All")
  res$log2fc = log2(res$mean_pos) - log2(res$mean_neg)
  res
  }) %>% 
  dplyr::ungroup() %>%
  dplyr::group_by(genename) %>%
  dplyr::filter(log2fc > 0) %>% # some model pairs have exactly the same pval/qval because they are mirrored. we keep the "positive" model
  dplyr::filter(pval == min(pval)) %>% # then we filter the top model 
  dplyr::ungroup()


models_annotations = as_tibble(do.call(rbind, lapply(unique(Monocle_run$model), function(model){
  data.frame(
    model=model,
    model_names_color=paste0(c('Red', 'Blue', 'Green', 'Purple', 'Orange', 'Yellow', 'Brown', 'Pink', 'Grey')[which(strsplit(substring(model, 2), split="")[[1]]==1)], collapse='/'),
    model_names_number=model,
    num_domain=sum(as.integer(strsplit(substring(model, 2), split="")[[1]])),
    first_domain_factor=factor(seq(numclusters_final)[which(strsplit(substring(model, 2), split="")[[1]]==1)[1]], levels=seq(numclusters_final))
  )
})))

#' Save all top models to file

Monocle_run_topModel %>% dplyr::left_join(models_annotations, by="model") %>% dplyr::select(model_names_color, genename, status, pval, qval, mean_pos, mean_neg, ratio_expr_pos, ratio_expr_neg, log2fc, model) %>% dplyr::arrange(model, pval) %>% write.csv(paste0(m$plot_folder, '/combinatorialDEtest_result.csv'), row.names=F)

#' <a href="./suppl_files/combinatorialDEtest_result.csv">Download (combinatorial) differential expression test results</a>
#'  
#' Filter genes

pval_thres = 1e-5
pos_pop_level_threshold = 0 #.5
pos_pop_ratio_threshold=0 #.15
log2fc_thres = 1

level_filtered_genes = Monocle_run_topModel %>% dplyr::filter(mean_pos >= pos_pop_level_threshold) %>% .$genename %>% unique%>% sort
ratio_filtered_genes = Monocle_run_topModel %>% dplyr::filter(ratio_expr_pos >= pos_pop_ratio_threshold) %>% .$genename %>% unique%>% sort
fc_filtered_genes = Monocle_run_topModel %>% dplyr::filter(log2fc > log2fc_thres | mean_neg == 0) %>% .$genename %>% unique%>% sort
selected_genes = Reduce(intersect, list(level_filtered_genes, fc_filtered_genes, ratio_filtered_genes))

Monocle_run_filtered = Monocle_run_topModel %>%
  dplyr::filter(genename %in% selected_genes) %>%
  dplyr::filter(pval <= pval_thres)

Monocle_run_genelist = Monocle_run_filtered %>%
  dplyr::group_by_(.dots="model") %>% 
  dplyr::summarise(genelist = list(genename)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(models_annotations, by="model")

testGMs = Monocle_run_genelist %>%
  # dplyr::arrange_at(c("num_domain", "first_domain_factor")) %>%
  dplyr::arrange_at(c("first_domain_factor", "num_domain")) %>%
  {"names<-"(.$genelist, .$model_names_color)}

m_oep$topCorr_DR$genemodules <- testGMs

output_basename = 'OEP_Final_DEtests'

m_oep$plotGeneModules(
  basename=output_basename,
  displayed.gms = 'topCorr_DR.genemodules',
  displayed.geneset=NA,
  use.dendrogram='Mansel',
  display.clusters='Mansel',
  file_settings=list(list(type='pdf', width=10, height=as.integer(0.02*length(unlist(m_oep$topCorr_DR$genemodules)) + 6))),
  data_status='Normalized',
  gene_transformations=c('logscaled'),
  extra_colors=cbind(
    pData(m_oep$expressionSet)$stage_colors),
  pretty.params=list("size_factor"=1, "ngenes_per_lines" = 8, "side.height.fraction"=.3),
  cluster_colors=clust.colors,
  show_cells_colors=FALSE,
  extra_legend=list("text"=names(stage_cols), "colors"=unname(stage_cols))
)

#' <a href="./suppl_files/OEP_Final_DEtests_topCorr_DR.genemodules_Normalized_logscaled.pdf">Download PDF</a>

#' <p align="center"><img src="./suppl_files/OEP_Final_DEtests_topCorr_DR.genemodules_Normalized_logscaled.png" width="100%"></p>
#'  

m_oep$writeGeneModules(basename=output_basename, gms='topCorr_DR.genemodules')

#' <a href="./suppl_files/OEP_Final_DEtests_topCorr_DR.genemodules.txt">Download modules of differentially expressed genes (by combinatorial categories)</a>
#'  

saveRDS(m_oep, paste0(m$plot_folder, '/m_oep_withDEtests.rds'))
# m_oep = readRDS(paste0(m$plot_folder, '/m_oep_withDEtests.rds'))



#' `r knitr::knit_exit()`
#'  

stop()



devtools::load_all("/camp/lab/briscoej/working/julien/Antler_dev/")

test = m_oep$copy()
m_oep = test$copy()


m_oep$carve(
  used_gms = m_oep$topCorr_DR$genemodules.selected,
  output_folder = paste0(m_oep$plot_folder,'/carvings_test_200_200/'),
  extra_pd = NULL,
  num_generations = 200,
  pop_size = 200,
  cluster_size_limit = NULL,
  cluster_size_limit_ratio = .05,
  expand_neighb_cutoff = 2
)

# Final score(s): -0.682610831832344 0.231087192287043

saveRDS(m_oep, paste0(output_path, '/m_oep_carved_biased.rds'))

test=readRDS(paste0(output_path, '/m_oep_carved_biased.rds'))

devtools::load_all('../Antler_dev/'); test2=test$copy(); test=test2$copy() # hack 

# Smoothed version
test2$readcounts_smoothed = getMagicDiffusionData(
  weighted_adjacency_matrix=test2$cellStateGraph$dist_matrix * igraph::as_adjacency_matrix(test2$cellStateGraph$graph, sparse=F),
  data_orig=test2$getReadcounts(data_status="Normalized"),
  graph_autotune = 3,
  epsilon = 1,
  t=5,
  rescale=.99
)


test2$calculatePTdynamics_No_kSPs(basename="AllComms", data_status='Smoothed', gene_list=some_genes)


test2$plotCellStateGraph(basename="TEST", plotted_gms=NULL, refresh_layout=FALSE, projection_method='gephiForceAtlas2', vertex.size=1000/test2$getNumberOfCells(), data_status='Normalized', genelist=some_genes, plot=c("genelists", "Mansel"))


#' <embed src=`r paste0(m$plot_folder, "AllComms_landmark_per_community.pdf")` type="application/pdf" width="100%" height="500">

#' <embed src=`r paste0(m$plot_folder, "AllComms_PTgraph_tree_plot_community_colors.pdf")` type="application/pdf" width="100%" height="500">

#' From the PT graph, we select two end populations

# End population 1 -> highest Tubb3, End population 2 -> highest Fabp7

select_leaves <- function(top_genes, currm){
  available_leaves = currm$PT %>%
    dplyr::filter(
      pt > 1 &
        name %in% names(which(degree(currm$PT_graph)==1))
    ) %>%
    as_tibble %>%
    dplyr::select("path_id") %>%
    dplyr::distinct() %>%
    unlist()
  
  do.call(rbind, lapply(top_genes, function(end_pop_marker){
    currm$PT %>% 
      dplyr::filter(path_id %in% available_leaves) %>%
      dplyr::group_by(path_id) %>%
      dplyr::filter(pt==max(pt) & gene_name == end_pop_marker) %>%
      dplyr::ungroup() %>% 
      dplyr::filter(mean.smooth.avg==max(mean.smooth.avg)) %>%
      dplyr::select(path_id, landmark_id)
  }))
}

selected_ends = select_leaves(c("VGLL2", "SOHO1"), test2)

print(selected_ends)

test2$plotPTdynamics(basename="AllComms1", gene_list=c("VGLL2", "SOHO1"), selected_paths=selected_ends$path_id, title="OEP bifuraction", subtitle='Test')

test2$plotPTdynamics(basename="AllComms2", gene_list=c("LMX1A", "TFAP2E"), selected_paths=selected_ends$path_id, title="OEP bifuraction", subtitle='Test')

test2$plotPTdynamics(basename="AllComms3", gene_list=c("FOXI3", "DLX5"), selected_paths=selected_ends$path_id, title="OEP bifuraction", subtitle='Test')





m_oep$carve(
  used_gms = m_oep$topCorr_DR$genemodules,
  output_folder = paste0(m_oep$plot_folder,'/carvings_test_200_200_unbiased/'),
  extra_pd = NULL,
  num_generations = 200,
  pop_size = 200,
  cluster_size_limit = NULL,
  cluster_size_limit_ratio = .05,
  expand_neighb_cutoff = 2
)

saveRDS(m_oep, paste0(output_path, '/m_oep_carved_unbiased.rds'))




#' # Extra

#' ### GO term analysis - Biological Processes

gProfileR::gprofiler(m$topCorr_DR$genemodules, organism='ggallus', png_fn=paste0(output_path, '/gprofiler_allGMs_BP.png'), include_graph=T, src_filter="GO:BP")
#'  

#' ### GO term analysis - Cellular Compartment

gProfileR::gprofiler(m$topCorr_DR$genemodules, organism='ggallus', png_fn=paste0(output_path, '/gprofiler_allGMs_CC.png'), include_graph=T, src_filter="GO:CC")
#'  

#' ### GO term analysis - Molecular Function

gProfileR::gprofiler(m$topCorr_DR$genemodules, organism='ggallus', png_fn=paste0(output_path, '/gprofiler_allGMs_MF.png'), include_graph=T, src_filter="GO:MF")
#'  

#' ### GO term analysis - Transcription Factor Binding

gProfileR::gprofiler(m$topCorr_DR$genemodules, organism='ggallus', png_fn=paste0(output_path, '/gprofiler_allGMs_TF.png'), include_graph=T, src_filter="TF")
#'  

#' ### Bulk comparison

#' Adding bulk samples + extra genes on initial heatmap I

mb = Antler$new(plot_folder=output_path, num_cores=4)

mb$loadDataset(folderpath='./dataset/scOEP_with_all_bulks/')

mb$setCurrentGeneNames(geneID_mapping_file=system.file("extdata", "Annotations/biomart_ensemblid_genename_ggallus.csv", package="Antler"))

mb$removeOutliers(lowread_thres = 5e5,  # select cells with more than 500000 reads 
                  genesmin = 1000,      # select cells expressing more than 1k genes
                  cellmin = 3,
                  data_status='Raw')          # select genes expressed in more than 3 cells)

annotations = list(
  "blank"=c('241112', '250184', '265102', '272111', '248185', '274173'),
  "bulk"=c('225110', '251172', '273103', '280110', '235161', '246161'),
  "human"=c('233111', '249196', '257101', '264112', '233185', '247173')
)
mb$excludeCellsFromIds(which(mb$getCellsNames() %in% unlist(annotations)))

# Remove cells having more than 6% of mitochondrial read counts
mb$removeGenesFromRatio(
  candidate_genes=grep('^MT-', mb$getGeneNames(), value=T),
  threshold = 0.06)

mb$normalize(method="Count-Per-Million")

mb$topCorr_DR$genemodules.selected <- m2$topCorr_DR$genemodules.selected

gl=c('FGF3', 'FGF10', 'FGF8', 'GSC', 'ACVR2A', 'EYA2', 'FZD1', 'NEUROD1', 'NEUROD2', 'NEUROD4', 'HES5', 'NOTCH1', 'DLL1', 'JAG2', 'ISL1', 'RPLP1', 'GAPDH')

mb$identifyCellClusters(method='hclust', used_genes="topCorr_DR.genemodules.selected", data_status='Normalized')

mb$plotGeneModules(
  basename='AllCellsAPrioriGMselection_extraGenes',
  displayed.gms = 'topCorr_DR.genemodules.selected',
  displayed.geneset=m$favorite_genes,
  use.dendrogram='hclust',
  display.clusters=NULL,
  file_settings=list(list(type='pdf', width=10, height=10)),
  data_status='Normalized',
  gene_transformations=c('log', 'logscaled'),
  extra_colors=apply(mb$getReadcounts('Normalized')[gl, ], 1, function(x){ unname(log(1+x)/max(log(1+x))) %>% {colorRampPalette(c("#0464DF", "#FFE800"))(n = 1000)[1+as.integer(999*.)]}}),
  pretty.params=list("size_factor"=5, "ngenes_per_lines" = 8, "side.height.fraction"=.3)
)


#' Adding bulk samples + extra genes on initial heatmap II
mb$dR$genemodules <- m2$dR$genemodules

mb$identifyCellClusters(method='hclust', clust_name="withBulk", used_genes="dR.genemodules", data_status='Normalized', numclusters=2)

pData(mb$expressionSet)$cells_colors <- generateCellColors(pData(mb$expressionSet), hue_shift=.3)

mb$plotGeneModules(
  basename='AllCellsManualGMselection_extraGenes',
  displayed.gms = 'dR.genemodules',
  displayed.geneset=m$favorite_genes,
  use.dendrogram='withBulk',
  display.clusters='withBulk',
  file_settings=list(list(type='pdf', width=10, height=10)),
  data_status='Normalized',
  gene_transformations=c('log', 'logscaled'),
  extra_colors=apply(mb$getReadcounts('Normalized')[gl, ], 1, function(x){ unname(log(1+x)/max(log(1+x))) %>% {colorRampPalette(c("#0464DF", "#FFE800"))(n = 1000)[1+as.integer(999*.)]}}),
  pretty.params=list("size_factor"=5, "ngenes_per_lines" = 6, "side.height.fraction"=.3)
)

#' ### Differential expression test OEP-derived (red) and non-OEP-derived (blue) cluster

DEtest_res = m2$runMonocle2DETestOnClusters(clustering_method="Mansel")

DEtest_res.2 = m2$getReadcounts('Normalized') %>%
  as.data.frame %>%
  tibble::rownames_to_column('genename') %>%
  tidyr::gather(cellname, value, -genename) %>%
  dplyr::left_join(
    data.frame(
      cellname=m2$getCellsNames(),
      clust_id=m2$cellClusters[['Mansel']]$cell_ids),
    by=c("cellname"="cellname")
  ) %>%
  dplyr::mutate(clust_id=paste0("c_", clust_id)) %>%
  dplyr::group_by(genename, clust_id) %>%
  dplyr::summarise(mean_val=mean(value)) %>%
  tidyr::spread(clust_id, mean_val) %>%
  dplyr::mutate(
    fold_change=log2(c_1/c_2),
    up_oep=ifelse(fold_change > 0, 1, 0)
  ) %>%
  dplyr::left_join(DEtest_res, by=c("genename"='gene_short_name')) %>%
  dplyr::arrange(up_oep, pval)

write.table(DEtest_res.2, file=paste0(output_path,'/DE_test_red_vs_blue.csv'), row.names=F, sep=";", quote=FALSE)

#' Select somes with pval, fold-change and mean level cutoffs

DE_genes_filtered = DEtest_res.2 %>%
  dplyr::filter(pval < 1e-5 & abs(fold_change) > 2 & max(c_1, c_2) > 10) %>%
  dplyr::arrange(up_oep, pval)

m2$dR$genemodules.selected = DE_genes_filtered%>%
  dplyr::group_by(up_oep) %>%
  dplyr::summarise(gl=list(genename)) %>%
  .$gl

#' Plot selected DE genes

m2$plotGeneModules(
  basename='AllCells_DEgenes',
  displayed.gms = 'dR.genemodules.selected',
  displayed.geneset="naked",
  use.dendrogram='Mansel',
  display.clusters='Mansel',
  file_settings=list(list(type='pdf', width=10, height=30)),
  data_status='Normalized',
  gene_transformations=c('log', 'logscaled'),
  pretty.params=list("size_factor"=3, "ngenes_per_lines" = 6, "side.height.fraction"=.6)
)

# render rmd file
#' 
#' #' `r knitr::knit_exit()`
#' #'  
#' 
#' stop()
#' 
#' # Render md file
#' source('~/working/julien/Antler/R/common_dev.R') # renderSource
#' 
#' # renderSource(source_r = "/camp/lab/briscoej/working/julien/OEPLineages/R_files/OEPLineages_analysis_galgal6.R", output_dir= paste0("/camp/lab/briscoej/working/julien/Muscle/inst/html/OEP_trajectories_", format(Sys.time(), '%y%m%d')), plot_path = "/camp/lab/briscoej/working/julien/OEPLineages/output/", output_repo="/camp/lab/briscoej/working/julien/OEPLineages/")
#' 
#' 
#' renderSource(source_r = "/camp/lab/briscoej/working/julien/OEPLineages/R_files/OEPLineages_analysis_galgal6.R", output_dir= paste0("/camp/lab/briscoej/working/julien/temp/OEP_trajectories_", format(Sys.time(), '%y%m%d')), plot_path = "/camp/lab/briscoej/working/julien/OEPLineages/output/galgal6/", output_repo="/camp/lab/briscoej/working/julien/OEPLineages/")
#' 

