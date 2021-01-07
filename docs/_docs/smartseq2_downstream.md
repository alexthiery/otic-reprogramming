---
layout: page
title: SmartSeq2 downstream
order: 6
---

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

</br>

### Data pre-processing and QC

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
  <button class="tablinks" onclick="openTab(event, 'Cells per stage-pre')">Cells per stage</button>
  <button class="tablinks" onclick="openTab(event, 'Readcounts per cell-pre')">Readcounts per cell</button>
  <button class="tablinks" onclick="openTab(event, 'Genes per cell-pre')">Genes per cell</button>
</div>

<div id="Cells per stage-pre" class="tabcontent">
  <embed src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/preQC_statistics_cellNumber_timepoint_by_timepoint.pdf" width="500" height="500" type="application/pdf">
</div>

<div id="Readcounts per cell-pre" class="tabcontent">
  <embed src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/preQC_statistics_counts_timepoint_by_timepoint.pdf" width="500" height="500" type="application/pdf">
</div>

<div id="Genes per cell-pre" class="tabcontent">
  <embed src="{{site.baseurl}}/assets/output/NF-downstream_analysis/smartseq_analysis/output/plots/preQC_statistics_counts_timepoint_by_timepoint.pdf" width="500" height="500" type="application/pdf">
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
