---
layout: page
label: SmartSeq2 scRNAseq
category: Downstream analysis
order: 3
---

## SmartSeq2 scRNAseq analysis

</br>

### Load packages, load data and pre-QC

<details><summary>Expand</summary>
<p>

Automatic switch for running pipeline through Nextflow or interactively in Rstudio

```R
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
```

</br>

Set paths and pipeline parameters. Load data, packages and custom functions.

```R
{
  if (opt$runtype == "user"){

    # load custom functions
    sapply(list.files('./NF-downstream_analysis/bin/custom_functions/', full.names = T), source)

    output_path = "./output/NF-downstream_analysis/smartseq_analysis/output/"
    plot_path = "./output/NF-downstream_analysis/smartseq_analysis/output/plots/"
    merged_counts_path = './output/NF-smartseq2_alignment/merged_counts/output/'
    genome_annotations_path = './output/NF-downstream_analysis/extract_gtf_annotations/'
    gfp_counts = './output/NF-smartseq2_alignment/merged_counts/output/'
    velocyto_input = './output/NF-smartseq2_alignment/velocyto/'

    # set cores
    ncores = 8

  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')

    # load custom functions
    sapply(list.files(opt$custom_functions, full.names = T), source)
    output_path = "./output/"
    plot_path = "./output/plots/"
    merged_counts_path = './'
    genome_annotations_path = './'
    gfp_counts = './'
    velocyto_input = './'

    # set cores
    ncores = opt$cores
  }

  dir.create(output_path, recursive = T)
  dir.create(plot_path, recursive = T)

  # load required packages
  library(Antler)
  library(velocyto.R)
  library(stringr)
  library(monocle)
  library(plyr)
  library(dplyr)
  library(ggsignif)
  library(cowplot)
  library(rstatix)
  library(extrafont)
}

# set pipeline params
seed=1
perp=5
eta=200

#' Stage colors
stage_cols = setNames(c("#BBBDC1", "#6B98E9", "#05080D"), c('8', '11', '15'))
```

</br>

Load data and initialise Antler object

```R
# Load and hygienize dataset
m = Antler$new(plot_folder=plot_path, num_cores=ncores)

# load in phenoData and assayData from ../dataset -> assayData is count matrix; phenoData is metaData (i.e. replicated, conditions, samples etc)
m$loadDataset(folderpath=merged_counts_path)

pData(m$expressionSet)$cells_colors = stage_cols[as.character(pData(m$expressionSet)$timepoint)]
```

</br>

Plot pre-QC metrics

```R
m$plotReadcountStats(data_status="Raw", by="timepoint", category="timepoint", basename="preQC", reads_name="read", cat_colors=unname(stage_cols))
```

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'Cells per stage-pre')">Cells per stage</button>
  <button class="tablinks" onclick="openTab(event, 'Readcounts per cell-pre')">Readcounts per cell</button>
  <button class="tablinks" onclick="openTab(event, 'Genes per cell-pre')">Genes per cell</button>
</div>

<div id="Cells per stage-pre" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/preQC_statistics_cellNumber_timepoint_by_timepoint.png">
</div>

<div id="Readcounts per cell-pre" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/preQC_statistics_counts_timepoint_by_timepoint.png">
</div>

<div id="Genes per cell-pre" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/preQC_statistics_geneCounts_timepoint_by_timepoint.png">
</div>

</br>

</br>

Some key genes are not annotated in the GTF - manually add these gene names to the annotations file

```R
# read in annotations file
gtf_annotations = read.csv(list.files(genome_annotations_path, pattern = '*gene_annotations.csv', full.names = T), stringsAsFactors = F)

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
```

</br>

Add list of key genes to Antler object

```R
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
```

</details>

---

</br>

### Data pre-processing and QC

</br>

<details><summary>Expand</summary>
<p>

Remove gene and cell outliers

```R
#' Remove outliers genes and cells
m$removeOutliers( lowread_thres = 5e5,   # select cells with more than 500000 reads
                  genesmin = 1000,       # select cells expressing more than 1k genes
                  cellmin = 3,           # select genes expressed in more than 3 cells)
                  data_status='Raw')
```

</br>

Remove control cells

```R
annotations = list(
  "blank"=c('241112', '250184', '265102', '272111', '248185', '274173'),
  "bulk"=c('225110', '251172', '273103', '280110', '235161', '246161'),
  "human"=c('233111', '249196', '257101', '264112', '233185', '247173')
)

m$excludeCellsFromIds(which(m$getCellsNames() %in% unlist(annotations)))
```

Remove cells with more than 6% of mitochondrial read counts

```R
m$removeGenesFromRatio(
  candidate_genes=grep('^MT-', m$getGeneNames(), value=T),
  threshold = 0.06
)
```

</br>

QC plots

```R
m$plotReadcountStats(data_status="Raw", by="timepoint", category="timepoint", basename="postQC", reads_name="read", cat_colors=unname(stage_cols))

```

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'Cells per stage-post')">Cells per stage</button>
  <button class="tablinks" onclick="openTab(event, 'Readcounts per cell-post')">Readcounts per cell</button>
  <button class="tablinks" onclick="openTab(event, 'Genes per cell-post')">Genes per cell</button>
</div>

<div id="Cells per stage-post" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/postQC_statistics_cellNumber_timepoint_by_timepoint.png">
</div>

<div id="Readcounts per cell-post" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/postQC_statistics_counts_timepoint_by_timepoint.png">
</div>

<div id="Genes per cell-post" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/postQC_statistics_geneCounts_timepoint_by_timepoint.png">
</div>

</br>

</br>

Normalize readcounts to CPM and remove genes with <10 CPM

```R
m$normalize(method="Count-Per-Million")
m$excludeUnexpressedGenes(min.cells=3, data_status='Normalized')
m$removeLowlyExpressedGenes(expression_threshold=1, selection_theshold=10, data_status='Normalized')
```

</details>

---

</br>

### Transcriptomic analysis of all cells

</br>

<details><summary>Expand</summary>
<p>

Generate gene-gene correlation matrix

```R
# change plot folder
curr_plot_folder = paste0(plot_path, "all_cells/")
dir.create(curr_plot_folder)

# "fastCor" uses the tcrossprod function which by default relies on BLAS to speed up computation. This produces inconsistent result in precision depending on the version of BLAS being used.
# see "matprod" in https://stat.ethz.ch/R-manual/R-devel/library/base/html/options.html
corr.mat = fastCor(t(m$getReadcounts(data_status='Normalized')), method="spearman")
```

Identification of modules of co-correlated genes

This feature is the prime reason for using the Antler package. This enables us to heirarchically cluster genes based on gene-gene correlation providing sets of co-correlated genes. Poor quality clusters are then filtered unbiasedly through iterative filtering/clustering.

For further information on Antler gene module identification, click [here](https://juliendelile.github.io/Antler/articles/Transcriptomic-summary.html#gene-modules-identification).

```R
#' ## Gene modules identification
m$identifyGeneModules(
  method="TopCorr_DR",
  corr=corr.mat,
  corr_t = 0.3, # gene correlation threshold
  topcorr_mod_min_cell=0, # default
  topcorr_mod_consistency_thres=0.4, # default
  topcorr_mod_skewness_thres=-Inf, # default
  topcorr_min_cell_level=5,
  topcorr_num_max_final_gms=40, # heuristically set to give sufficient cluster diversity and of reasonable size in order to select modules from gene candidate list
  data_status='Normalized'
)
# corr_t = gene correlation threshold
# topcorr_mod_min_cell = minimum number of cells expressing gene module
# topcorr_mod_consistency_thres = proportion of cells expressing gene module (binarised z-scored log-transformed normalised expression levels)
# topcorr_mod_skewness_thres = test whether genes are expressed only in subset of cells
# topcorr_min_cell_level = minimum expression level in cells which express gene module
# topcorr_num_max_final_gms = maximum number of final gene modules

names(m$topCorr_DR$genemodules) <- paste0("GM ", seq(length(m$topCorr_DR$genemodules)))
```

Plot gene modules

```R
# identify cell clusters based on remaining genes
m$identifyCellClusters(method='hclust', used_genes="topCorr_DR.genemodules", data_status='Normalized')

m$plotGeneModules(
  basename='AllCells',
  curr_plot_folder = curr_plot_folder,
  displayed.gms = 'topCorr_DR.genemodules',
  displayed.geneset=NA,
  use.dendrogram='hclust',
  display.clusters=NULL,
  file_settings=list(list(type='pdf', width=20, height=20)),
  data_status='Normalized',
  gene_transformations='logscaled',
  pretty.params=list("size_factor"=3, "ngenes_per_lines" = 0, "side.height.fraction"=.15),
  display.legend = FALSE,
  extra_legend=list("text"=c('ss8-9', 'ss11-12', 'ss14-15'), "colors"=unname(stage_cols))
)

m$writeGeneModules(basename='AllCells_allGms', gms='topCorr_DR.genemodules', folder_path = curr_plot_folder)
```

<a href="{{ site.baseurl }}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/AllCells_topCorr_DR.genemodules_Normalized_logscaled.pdf">Download
PDF</a>

<a href="{{ site.baseurl }}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/AllCells_allGms_topCorr_DR.genemodules.txt">Download
gene module list</a>

![]({{ site.baseurl }}{% link /assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/AllCells_topCorr_DR.genemodules_Normalized_logscaled.png %})

</br>

Filter out cells with low summary counts

```R
m2 = m$copy()
m2$excludeCellsFromIds(m$getReadcounts('Normalized')[unlist(m$topCorr_DR$genemodules),] %>% colSums %>% {as.numeric(scale(log(.), center=TRUE, scale=T))} %>% {which(. < -1.5)})
m2$excludeUnexpressedGenes(min.cells=1, data_status="Normalized", verbose=TRUE)
```

</br>

Select gene modules based on the presence of genes known to be involved in differentiation and re-cluster cells

```R
bait_genes = c("HOXA2", "PAX6", "SOX2", "MSX1", "Pax3", "SALL1", "ETS1", "TWIST1", "HOMER2", "LMX1A", "VGLL2", "EYA2", "PRDM1", "FOXI3", "NELL1", "DLX5", "SOX8", "SOX10", "SOHO-1", "IRX4", "DLX6")

m2$dR$genemodules = Filter(function(x){any(bait_genes %in% x)}, m2$topCorr_DR$genemodules)

# cluster into 5 clusters
m2$identifyCellClusters(method='hclust', clust_name="Mansel", used_genes="dR.genemodules", data_status='Normalized', numclusters=5)

# set cluster colours
clust.colors <- c('#da70d6', '#c71585', '#b0c4de', '#afeeee', '#5f9ea0')

# read in GFP counts
gfp_counts = read.table(file=paste0(gfp_counts, 'gfpData.csv'), header=TRUE, check.names=FALSE)
```

</br>

Plot final clustering of all cells

```R
m2$plotGeneModules(
  basename='AllCellsManualGMselection',
  curr_plot_folder = curr_plot_folder,
  displayed.gms = 'dR.genemodules',
  displayed.geneset=NA,
  use.dendrogram='Mansel',
  file_settings=list(list(type='pdf', width=10, height=10)),
  data_status='Normalized',
  gene_transformations=c('log', 'logscaled'),
  extra_colors=cbind(
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
pretty.params=list("size_factor"=2, "side.height.fraction"=0.5),
display.legend = FALSE,
extra_legend=list("text"=c('ss8-9', 'ss11-12', 'ss14-15'), "colors"=unname(stage_cols))
)

m2$writeGeneModules(basename='AllCells_baitGMs', gms='dR.genemodules', folder_path = curr_plot_folder)
```

<a href="{{ site.baseurl }}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/AllCellsManualGMselection_dR.genemodules_Normalized_logscaled.pdf">Download
PDF</a>

<a href="{{ site.baseurl }}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/AllCells_baitGMs_dR.genemodules">Download
gene module list</a>

![]({{ site.baseurl }}{% link /assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/AllCellsManualGMselection_dR.genemodules_Normalized_logscaled.png %})

</br>

Plot tSNEs of cell clusters and developmental stage

```R
curr_plot_folder = paste0(plot_path, "all_cells/tsne/")
dir.create(curr_plot_folder)


png(paste0(curr_plot_folder, 'allcells_clusters_TSNE.png'), height = 15, width = 15, family = 'Arial', units = 'cm', res = 400)
tsne_plot(m2, m2$dR$genemodules, seed=seed, colour_by=m2$cellClusters$Mansel$cell_ids, colours=clust.colors, perplexity=perp, eta=eta) +
  ggtitle('Clusters') +
  theme(plot.title = element_text(hjust = 0.5))
graphics.off()


png(paste0(curr_plot_folder, 'allcells_stage_TSNE.png'), height = 15, width = 15, family = 'Arial', units = 'cm', res = 400)
tsne_plot(m2, m2$dR$genemodules, seed=seed, colour_by = pData(m2$expressionSet)$timepoint, colours = stage_cols, perplexity=perp, eta=eta) +
  ggtitle('Developmental stage') +
  theme(plot.title = element_text(hjust = 0.5))
graphics.off()
```

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'allcells_clusters_TSNE')">Cell clusters tSNE</button>
  <button class="tablinks" onclick="openTab(event, 'allcells_stage_TSNE')">Developmental stage tSNE</button>
</div>

<div id="allcells_clusters_TSNE" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/tsne/allcells_clusters_TSNE.png">
</div>

<div id="allcells_stage_TSNE" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/tsne/allcells_stage_TSNE.png">
</div>

</br>

Plot tSNEs of genes of interest

```R
gene_list = c('SOX2', 'SOX10', 'SOX8', 'PAX7', 'PAX2', 'LMX1A', 'SOX21', 'SIX1')
for(gn in gene_list){
  png(paste0(curr_plot_folder, gn, '_TSNE.png'), height = 15, width = 15, family = 'Arial', units = 'cm', res = 400)
  print(tsne_plot(m2, m2$dR$genemodules, seed=seed, colour_by = as.integer(1+100*log10(1+m2$getReadcounts(data_status='Normalized')[gn,]) / max(log10(1+m2$getReadcounts(data_status='Normalized')[gn,]))),
            colours = c("grey", "darkmagenta"), perplexity=perp, eta=eta) +
          ggtitle(gn) +
          theme(plot.title = element_text(hjust = 0.5)))
  graphics.off()
}


png(paste0(curr_plot_folder, 'GFP_TSNE.png'), height = 15, width = 15, family = 'Arial', units = 'cm', res = 400)
tsne_plot(m2, m2$dR$genemodules, seed=seed, colour_by = as.integer(1+100*log10(1+gfp_counts[m2$getCellsNames()]) / max(log10(1+gfp_counts[m2$getCellsNames()]))),
                colours = c("grey", "darkgreen"), perplexity=perp, eta=eta) +
  ggtitle('GFP') +
  theme(plot.title = element_text(hjust = 0.5))
graphics.off()
```

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'all_SOX2')">Sox2</button>
  <button class="tablinks" onclick="openTab(event, 'all_SOX10')">Sox10</button>
  <button class="tablinks" onclick="openTab(event, 'all_SOX8')">Sox8</button>
  <button class="tablinks" onclick="openTab(event, 'all_PAX7')">Pax7</button>
  <button class="tablinks" onclick="openTab(event, 'all_PAX2')">Pax2</button>
  <button class="tablinks" onclick="openTab(event, 'all_LMX1A')">Lmx1a</button>
  <button class="tablinks" onclick="openTab(event, 'all_SOX21')">Sox21</button>
  <button class="tablinks" onclick="openTab(event, 'all_SIX1')">Six1</button>
  <button class="tablinks" onclick="openTab(event, 'all_GFP')">GFP</button>

</div>

<div id="all_SOX2" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/tsne/SOX2_TSNE.png">
</div>

<div id="all_SOX10" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/tsne/SOX10_TSNE.png">
</div>

<div id="all_SOX8" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/tsne/SOX8_TSNE.png">
</div>

<div id="all_PAX7" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/tsne/PAX7_TSNE.png">
</div>

<div id="all_PAX2" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/tsne/PAX2_TSNE.png">
</div>

<div id="all_LMX1A" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/tsne/LMX1A_TSNE.png">
</div>

<div id="all_SOX21" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/tsne/SOX21_TSNE.png">
</div>

<div id="all_SIX1" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/tsne/SIX1_TSNE.png">
</div>

<div id="all_GFP" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/tsne/GFP_TSNE.png">
</div>

</br>

Plot dotplot for genes of interest

```R
curr_plot_folder = paste0(plot_path, "all_cells/")
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
  dplyr::mutate(value = scale(value)) %>%
  # calculate mean expression
  dplyr::group_by(genename, celltype) %>%
  dplyr::mutate('Scaled Average Expression'=mean(value)) %>%
  dplyr::distinct(genename, celltype, .keep_all=TRUE) %>%
  dplyr::ungroup() %>%
  # make factor levels to order genes in dotplot
  dplyr::mutate(genename = factor(genename, levels = gene_list)) %>%
  # make factor levels to order cells in dotplott
  dplyr::mutate(celltype = factor(celltype, levels = rev(c("OEP", "Late Placodal", "Neural Crest", "Neural", "Mesodermal"))))

png(paste0(curr_plot_folder, "all_cells_dotplot.png"), width=25, height=10, family = 'Arial', units = "cm", res = 400)
ggplot(dotplot_data, aes(x=genename, y=celltype, size=`Proportion of Cells Expressing`, color=`Scaled Average Expression`)) +
  geom_count() +
  scale_size_area(max_size=5) +
  scale_x_discrete(position = "top") + xlab("") + ylab("") +
  scale_color_gradient(low = "grey90", high = "blue") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size=9), axis.text.y = element_text(colour = clust.colors[c(4,3,5,2,1)], face = 'bold', size = 12),
        legend.position="bottom", legend.box = "horizontal", plot.margin=unit(c(0,1,0,0),"cm"))
graphics.off()
```

![]({{ site.baseurl }}{% link /assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/all_cells/all_cells_dotplot.png %})

</details>

---

</br>

### Transcriptomic analysis of OEP cells

<details><summary>Expand</summary>
<p>

Clusters 3-5 are composed of non-oep derived populations (they are also mostly PAX2 negative)

We exclude these cells from the subsequent analysis

```R
curr_plot_folder = paste0(plot_path, "oep_subset/")
dir.create(curr_plot_folder)

m_oep = m2$copy()
m_oep$excludeCellFromClusterIds(cluster_ids=c(3:5), used_clusters='Mansel', data_status='Normalized')

#' Some genes may not be expressed any more in the remaining cells
m_oep$excludeUnexpressedGenes(min.cells=1, data_status="Normalized", verbose=TRUE)
m_oep$removeLowlyExpressedGenes(expression_threshold=1, selection_theshold=10, data_status='Normalized')
```

Calculate and plot gene modules for OEPs

```R

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

m_oep$writeGeneModules(folder_path = curr_plot_folder, basename='OEP_allGms', gms='topCorr_DR.genemodules')

m_oep$identifyCellClusters(method='hclust', used_genes="topCorr_DR.genemodules", data_status='Normalized')

m_oep$plotGeneModules(
  curr_plot_folder = curr_plot_folder,
  basename='OEP_allGms',
  displayed.gms = 'topCorr_DR.genemodules',
  displayed.geneset=NA,
  use.dendrogram='hclust',
  display.clusters=NULL,
  file_settings=list(list(type='pdf', width=10, height=10)),
  data_status='Normalized',
  gene_transformations='logscaled',
  pretty.params=list("size_factor"=2, "ngenes_per_lines" = 8, "side.height.fraction"=0.15),
  display.legend = FALSE,
  extra_legend=list("text"=c('ss8-9', 'ss11-12', 'ss14-15'), "colors"=unname(stage_cols))
)
```

<a href="{{ site.baseurl }}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/OEP_allGms_topCorr_DR.genemodules_Normalized_logscaled.pdf">Download
PDF</a>

<a href="{{ site.baseurl }}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/OEP_allGms_topCorr_DR.genemodules.txt">Download
gene module list</a>

![]({{ site.baseurl }}{% link /assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/OEP_allGms_topCorr_DR.genemodules_Normalized_logscaled.png %})

Select gene modules based on the presence of genes known to be involved in otic and epibranchial differentiation and re-cluster cells

```R
bait_genes = c("HOMER2", "LMX1A", "SOHO-1", "SOX10", "VGLL2", "FOXI3", 'ZNF385C', 'NELL1', "CXCL14", "EYA4")

m_oep$topCorr_DR$genemodules.selected = Filter(function(x){any(bait_genes %in% x)}, m_oep$topCorr_DR$genemodules)

# cluster into 5 clusters
m_oep$identifyCellClusters(method='hclust', clust_name="Mansel", used_genes="topCorr_DR.genemodules.selected", data_status='Normalized', numclusters=5)

clust.colors <- c('#ffa07a', '#f55f20', '#dda0dd', '#48d1cc', '#b2ffe5')

m_oep$plotGeneModules(
  curr_plot_folder = curr_plot_folder,
  basename='OEP_GMselection',
  displayed.gms = 'topCorr_DR.genemodules.selected',
  displayed.geneset=NA,
  use.dendrogram='Mansel',
  file_settings=list(list(type='pdf', width=10, height=5)),
  data_status='Normalized',
  gene_transformations=c('log', 'logscaled'),
  extra_colors=cbind(
    m_oep$cellClusters$Mansel$cell_ids %>% clust.colors[.]
  ),
  display.legend = FALSE,
  pretty.params=list("size_factor"=2, "ngenes_per_lines" = 6, "side.height.fraction"=0.5),
  extra_legend=list("text"=c('ss8-9', 'ss11-12', 'ss14-15'), "colors"=unname(stage_cols))
)

# add gene modules txt
m_oep$writeGeneModules(folder_path = curr_plot_folder, basename='OEP_GMselection', gms='topCorr_DR.genemodules')
```

<a href="{{ site.baseurl }}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/OEP_GMselection_topCorr_DR.genemodules.selected_Normalized_logscaled.pdf">Download
PDF</a>

<a href="{{ site.baseurl }}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/OEP_GMselection_topCorr_DR.genemodules.txt">Download
gene module list</a>

![]({{ site.baseurl }}{% link /assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/OEP_GMselection_topCorr_DR.genemodules.selected_Normalized_logscaled.png %})

</br>

Plot tSNEs of cell clusters and developmental stage

```R
curr_plot_folder = paste0(plot_path, 'oep_subset/tsne/')
dir.create(curr_plot_folder)

png(paste0(curr_plot_folder, 'OEP_clusters_TSNE.png'), height = 15, width = 15, family = 'Arial', units = 'cm', res = 400)
tsne_plot(m_oep, m_oep$topCorr_DR$genemodules.selected, seed=seed, colour_by=m_oep$cellClusters[['Mansel']]$cell_ids, colours=clust.colors, perplexity=perp, eta=eta) +
ggtitle('Clusters') +
theme(plot.title = element_text(hjust = 0.5))
graphics.off()

png(paste0(curr_plot_folder, 'OEP_stage_TSNE.png'), height = 15, width = 15, family = 'Arial', units = 'cm', res = 400)
tsne_plot(m_oep, m_oep$topCorr_DR$genemodules.selected, seed=seed, colour_by = pData(m_oep$expressionSet)$timepoint, colours = stage_cols, perplexity=perp, eta=eta) +
ggtitle('Developmental stage') +
theme(plot.title = element_text(hjust = 0.5))
graphics.off()
```

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'oep_clusters_TSNE')">Cell clusters tSNE</button>
  <button class="tablinks" onclick="openTab(event, 'oep_stage_TSNE')">Developmental stage tSNE</button>
</div>

<div id="oep_clusters_TSNE" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/OEP_clusters_TSNE.png">
</div>

<div id="oep_stage_TSNE" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/OEP_stage_TSNE.png">
</div>

</br>

Plot expression of bait genes plus genes from bulk RNAseq and literature on tsne

```R
gene_list = c(bait_genes, 'SOX8', 'PAX2', 'TFAP2E', 'SIX1', 'ZBTB16', 'FOXG1', 'PDLIM1')
for(gn in gene_list){
png(paste0(curr_plot_folder, gn, '\_TSNE.png'), height = 15, width = 15, family = 'Arial', units = 'cm', res = 400)
print(tsne_plot(m_oep, m_oep$topCorr_DR$genemodules.selected, seed=seed, colour_by = as.integer(1+100*log10(1+m_oep$getReadcounts(data_status='Normalized')[gn,]) / max(log10(1+m_oep$getReadcounts(data_status='Normalized')[gn,]))),
colours = c("grey", "darkmagenta"), perplexity=perp, eta=eta) +
ggtitle(gn) +
theme(plot.title = element_text(hjust = 0.5)))
graphics.off()
}
```

<a href="{{ site.baseurl }}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/oep_GOI_tSNEs.zip">Download
all tSNEs</a>

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'oep_SOX8')">Sox2</button>
  <button class="tablinks" onclick="openTab(event, 'oep_PAX2')">Pax2</button>
  <button class="tablinks" onclick="openTab(event, 'oep_TFAP2E')">TFAP2E</button>
  <button class="tablinks" onclick="openTab(event, 'oep_SIX1')">Six1</button>
  <button class="tablinks" onclick="openTab(event, 'oep_ZBTB16')">ZBTB16</button>
  <button class="tablinks" onclick="openTab(event, 'oep_FOXG1')">FOXG1</button>
  <button class="tablinks" onclick="openTab(event, 'oep_PDLIM1')">PDLIM1</button>

</div>

<div id="oep_SOX8" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/SOX8_TSNE.png">
</div>

<div id="oep_PAX2" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/PAX2_TSNE.png">
</div>

<div id="oep_TFAP2E" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/TFAP2E_TSNE.png">
</div>

<div id="oep_SIX1" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/SIX1_TSNE.png">
</div>

<div id="oep_ZBTB16" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/ZBTB16_TSNE.png">
</div>

<div id="oep_FOXG1" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/FOXG1_TSNE.png">
</div>

<div id="oep_PDLIM1" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/PDLIM1_TSNE.png">
</div>

</br>

Plot tSNE co-expression plots

```R
gene_pairs <- list(c("FOXI3", "LMX1A"), c("TFAP2E", "LMX1A"), c("FOXI3", "SOX8"), c("TFAP2E", "SOX8"))
lapply(gene_pairs, function(x) {plot_tsne_coexpression(m_oep, m_oep$topCorr_DR$genemodules.selected, gene1 = x[1], gene2 = x[2], plot_folder = curr_plot_folder,
seed=seed, perplexity=perp, pca=FALSE, eta=eta, height = 15, width = 22, res = 400, units = 'cm', family = 'Arial')})

```

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'FOXI3_LMX1A')">FOXI3/LMX1A</button>
  <button class="tablinks" onclick="openTab(event, 'TFAP2E_LMX1A')">TFAP2E/LMX1A</button>
  <button class="tablinks" onclick="openTab(event, 'FOXI3_SOX8')">FOXI3/SOX8</button>
  <button class="tablinks" onclick="openTab(event, 'TFAP2E_SOX8')">TFAP2E/SOX8</button>

</div>

<div id="FOXI3_LMX1A" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/FOXI3_LMX1A_co-expression_TSNE.png">
</div>

<div id="TFAP2E_LMX1A" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/TFAP2E_LMX1A_co-expression_TSNE.png">
</div>

<div id="FOXI3_SOX8" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/FOXI3_SOX8_co-expression_TSNE.png">
</div>

<div id="TFAP2E_SOX8" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/oep_subset/tsne/TFAP2E_SOX8_co-expression_TSNE.png">
</div>

</details>

---

</br>

### Pseudotime using [Monocle2](http://cole-trapnell-lab.github.io/monocle-release/docs/)

<details><summary>Expand</summary>
<p>

In order to study the process of differentiation from OEP to Otic and Epibranchial lineages, we order the cells in pseudotime using our subset gene modules identified as important for this process.

```R
curr_plot_folder = paste0(plot_path, "monocle_plots/")
dir.create(curr_plot_folder)

# gene modules for biological process of interest
monocle.input_dims = unlist(m_oep$topCorr_DR$genemodules.selected)

# input Antler data into Monocle
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

# Dimensionality reduction and order cells along pseudotime from earliest "State" (ie DDRTree branch)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- setOrderingFilter(HSMM, monocle.input_dims)
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
HSMM <- orderCells(HSMM, root_state = which.max(table(pData(HSMM)$State, pData(HSMM)$timepoint)[, "8"]))
```

Plot cell stage over projected pseudotime coordinates

```R
png(paste0(curr_plot_folder, 'Monocle_DDRTree_samples.png'), height = 15, width = 15, family = 'Arial', units = 'cm', res = 400)
ggplot(data.frame('Component 1' = reducedDimS(HSMM)[1,], 'Component 2' = reducedDimS(HSMM)[2,],
                  'colour_by' = as.factor(pData(m_oep$expressionSet)$timepoint), check.names = FALSE), aes(x=`Component 1`, y=`Component 2`, color=colour_by)) +
  geom_point() +
  scale_color_manual(values = stage_cols) +
  theme_classic() +
  theme(legend.position = "none", axis.ticks=element_blank(), axis.text = element_blank()) +
  ggtitle('Developmental stage') +
  theme(plot.title = element_text(hjust = 0.5))
graphics.off()
```

<div class="tab">
  <button class="tablinks" onclick="openTab(event, 'Monocle_DDRTree_samples')">Pseudotime: developmental stage</button>
  <button class="tablinks" onclick="openTab(event, 'Monocle_DDRTree_Clusters')">Pseudotime: clusters</button>

</div>

<div id="Monocle_DDRTree_samples" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/monocle_plots/Monocle_DDRTree_samples.png">
</div>
<div id="Monocle_DDRTree_Clusters" class="tabcontent">
  <img src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/monocle_plots/Monocle_DDRTree_Clusters.png">
</div>

</br>

</details>
