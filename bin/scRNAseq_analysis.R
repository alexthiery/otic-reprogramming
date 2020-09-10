custom_functions = "./bin/custom_functions/"
input_path = "./output/merged_counts/"
plot_path = "./output/antler/plots/"
rds_path = "./output/antler/rds_files/"
dir.create(plot_path, recursive = T)
dir.create(rds_path, recursive = T)

# loadload required packages
library(Antler)

# load custom functions
sapply(list.files(custom_functions, full.names = T), source)

#' Stage colors
stage_cols = setNames(c("#BBBDC1", "#6B98E9", "#05080D"), c('8', '11', '15'))

#' # Load and hygienize dataset
m = Antler$new(plot_folder=plot_path, num_cores=6)

# load in phenoData and assayData from ../dataset -> assayData is count matrix; phenoData is metaData (i.e. replicated, conditions, samples etc)
m$loadDataset(folderpath= input_path)

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

m$writeGeneModules(basename='AllCells_allGms', gms='topCorr_DR.genemodules')

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

