# docker run -e PASSWORD=test -p 8787:8787 --rm -v ~/dev/repos/otic-reprogramming/:/home/rstudio alexthiery/rscript:latest

custom_functions = "./bin/custom_functions/"
output_dir = "./output/antler"
input_dir = "./output/merged_counts"

library(devtools)
library(Antler)

#' load custom functions
sapply(list.files(custom_functions, full.names = T), source)

#' # Load and hygienize dataset
m <- Antler$new(output_folder = output_dir, num_cores = 4)
m$load_dataset(folder_path = input_dir)
m$convert_ensembl_id_to_name(biomart_dataset = "ggallus_gene_ensembl")

#' Add custom colours for timepoints
m$add_color_map(
  name = "timepoint",
  content = setNames(c("#BBBDC1", "#6B98E9", "#05080D"), c('8', '11', '15')))

#' Plot cell QC
m$plot_QC(
  data_status = "Raw",
  feature_1   = "timepoint",
  reads_type  = "read")

#' Remove outliers genes and cells
m$remove_outliers(
  lowread_thres = 5e5, # select cells with more than 500000 reads 
  min_genes = 1000,    # select cells expressing more than 1k genes
  min_cells = 3,       # select genes expressed in more than 3 cells)
  data_status = 'Raw')

annotations = list(
  "blank"=c('241112', '250184', '265102', '272111', '248185', '274173'),
  "bulk"=c('225110', '251172', '273103', '280110', '235161', '246161'),
  "human"=c('233111', '249196', '257101', '264112', '233185', '247173')
)

m$remove_cells(which(m$cell_names() %in% unlist(annotations)))


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

m$removeGenesFromRatio(
  candidate_genes=grep('^MT-', m$gene_names(), value=T),
  threshold = 0.06
)


#' Normalise by CPM
m$normalize(method = 'CPM')

#' Remove genes with an average level of 10 CPM or less, among positive cells
m$exclude_unexpressed_genes(min_cells = 3, min_level = 10, data_status='Normalized')


saveRDS(m, file=paste0(output_dir, 'm_clean.rds'))
#' Read m_clean RDS file if needed
#' m=readRDS(paste0(rds_path, 'm_clean.rds'))




m$gene_modules$identify(
  name                  = "initGMs",
  corr_t                = 0.3,  # the Spearman correlation treshold
  corr_min              = 0,    # min. number of genes a gene must correlate with
  mod_consistency_thres = 0.4,  # ratio of expressed genes among "positive" cells
  min_cell_level        = 5,
  num_max_final_gms     = 40,
  process_plots         = TRUE)

##### UP TO HERE #####

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


