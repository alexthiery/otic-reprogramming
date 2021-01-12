---
layout: page
label: SmartSeq2 alignment
order: 2
---

## SmartSeq2 alignment

</br>

**If you would like to re-run the alignment, follow the instructions [here]({{site.baseurl}}/general/quick_start#smartseq).**

Our SmartSeq2 single cell RNAseq data is aligned using a custom Nextflow pipeline. This pipeline is built using [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html), which allows for modularisation and re-usability of each tool.

Each process is configured to run within a separate container, meaning that package versions are hard coded into the pipeline. In order to change these package versions, the container used will need to be changed within the individual process.

Each of the processes used in the pipeline can be found by following the going to the path in the `include {process} from "<PATH>"` statements in each of the workflows below.

---

</br>

### Workflows

Our custom pipeline integrates both standard alignment and RNA velocity ([Velocyto](http://velocyto.org/)).

For simplicity, the pipeline is split into four sub-workflows:

- [Modify genome to add GFP sequence](#GFP)
- [Align reads](#align)
- [Process readcounts for downstream](#process_reads)
- [Run RNA velocity](#velocyto)

</br>

The main workflow integrates each of these sub-workflows.

<details open><summary class="simple">Collapse main workflow</summary>
<p>

```groovy
// Define DSL2
nextflow.enable.dsl=2

include {add_gfp} from "$baseDir/../workflows/add_gfp/main.nf"
include {smartseq2_align} from "$baseDir/../workflows/scRNAseq_alignment/main.nf"
include {process_counts} from "$baseDir/../workflows/process_counts/main.nf"
include {velocyto_smartseq2} from "$baseDir/../workflows/velocyto_smartseq2/main.nf"

Channel
    .value(file(params.fasta, checkIfExists: true))
    .set {ch_fasta}

Channel
    .value(file(params.gtf, checkIfExists: true))
    .set {ch_gtf}

Channel
    .value(file(params.gfp_seq, checkIfExists: true))
    .set {ch_gfp_seq}

workflow {
    // add gfp to genome and gtf
    add_gfp (ch_fasta, ch_gtf, ch_gfp_seq)

    // align using smartseq2 workflow
    smartseq2_align (add_gfp.out.genome, add_gfp.out.gtf, params.input)

    // merge and extract gfp counts
    process_counts (smartseq2_align.out.htseq_count_files)

    // run velocyto
    velocyto_smartseq2 (smartseq2_align.out.bam_files, ch_gtf)

    // view output
    process_counts.out.processed_counts | view
    velocyto_smartseq2.out.velocyto_counts | view
}
```

</details>

---

</br>

**_Modify genome to add GFP sequence_**<a name="GFP"></a>

This sub-workdflow takes a GFP sequence from the base repository and appends the genome provided in order to determine the number of GFP reads downstream.

<details open><summary class="simple">Collapse sub-workflow</summary>
<p>

```groovy
// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {add_genome_gfp; add_gtf_gfp} from "$baseDir/../modules/genome-tools/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow add_gfp {
    take:
        genome
        gtf
        gfp_seq

    main:
        add_genome_gfp (params.modules['add_genome_gfp'], genome, gfp_seq)
        add_gtf_gfp (params.modules['add_gtf_gfp'], gtf, gfp_seq)

    emit:
        genome = add_genome_gfp.out
        gtf = add_gtf_gfp.out
}
```

</details>

---

</br>

**_Align reads_**<a name="align"></a>

First we trim adaptor sequences using Cutadapt v2.10. We then align reads to GalGal6 with HISAT2 v2.2.1. Samtools v1.10 is used to convert SAM to BAM, before finally obtaining gene counts with HTSeq v0.12.4.

<details open><summary class="simple">Collapse sub-workflow</summary>
<p>

```groovy
// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {smartseq2_fastq_metadata} from "$baseDir/../luslab-nf-modules/tools/metadata/main.nf"
include {cutadapt} from "$baseDir/../luslab-nf-modules/tools/cutadapt/main.nf"
include {hisat2_build; hisat2_splice_sites; hisat2_splice_align} from "$baseDir/../luslab-nf-modules/tools/hisat2/main.nf"
include {samtools_view as samtools_view_a;samtools_view as samtools_view_b; samtools_sort} from "$baseDir/../luslab-nf-modules/tools/samtools/main.nf"
include {htseq_count} from "$baseDir/../luslab-nf-modules/tools/htseq/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow smartseq2_align {
    take:
        genome
        gtf
        sample_csv

    main:
        smartseq2_fastq_metadata (sample_csv)

        cutadapt (params.modules['cutadapt'], smartseq2_fastq_metadata.out)
        hisat2_build ( params.modules['hisat2_build'], genome )
        hisat2_splice_sites ( params.modules['hisat2_splice_sites'], gtf )
        hisat2_splice_align ( params.modules['hisat2_splice_align'], cutadapt.out.fastq, hisat2_build.out.genome_index.collect(), hisat2_splice_sites.out.splice_sites.collect() )

        samtools_view_a ( params.modules['samtools_view_a'], hisat2_splice_align.out.sam )
        samtools_sort ( params.modules['samtools_sort'], samtools_view_a.out.bam )
        samtools_view_b ( params.modules['samtools_view_b'], samtools_sort.out.bam )

        htseq_count ( params.modules['htseq_count'], samtools_view_b.out.bam, gtf )

    emit:
        bam_files = samtools_sort.out.bam
        htseq_count_files = htseq_count.out.counts
}
```

</details>

---

</br>

**_Process read counts for input into Antler_**<a name="process_reads"></a>

The merge counts process is a custom R script which formats cell names and generates the phenoData.csv, assayData.csv and gfpData.csv required by our Antler downstream pipeline.

<details open><summary class="simple">Collapse sub-workflow</summary>
<p>

```groovy
// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {r_analysis as merge_counts} from "$baseDir/../modules/r_analysis/main.nf"


/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow process_counts {
    take:
        count_files

    main:
        // group counts into a single channel for merge cell counts
        ch_all_counts = count_files
            .map { [file(it[1], checkIfExists: true)] }
            .collect()

        // merge cell counts into csv
        merge_counts (params.modules['merge_counts'], ch_all_counts)

    emit:
        processed_counts = merge_counts.out
}
```

</details>

</br>

<details><summary class="simple">Collapse merge_counts.R</summary>
<p>

```R
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
    output_path = "./output/merge_counts/output/"
    input_files <- list.files("./testData/process_counts/", pattern = "*.txt", full.names = T)

  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through nextflow\n')

    output_path = "output/"
    input_files <- list.files("./", pattern = "*.txt", full.names = T)
  }
  dir.create(output_path, recursive = T)

  library(plyr)
}

assayData = do.call(cbind, lapply(input_files, function(f){
	df = read.table(f, row.names=1)
	df[seq(nrow(df)-5),, drop=F]
}))

# ss11 ss15 WTCHG_528508_201189.txt -> 201189
# ss89 TSS_P2_E2_13-1234.txt -> P2E2
phenoData <- setNames(ldply(basename(input_files), .fun = function(x) {
  if (grepl("ss11", x)) {
    c(strsplit(tools::file_path_sans_ext(x), split = "_")[[1]][[3]], 11, 1, "sc")
  } else  if (grepl("ss15", x)) {
    c(strsplit(tools::file_path_sans_ext(x), split = "_")[[1]][[3]], 15, 1, "sc")
  } else {
    spl1 = strsplit(tools::file_path_sans_ext(x), split = "_")[[1]]
    c(paste0(spl1[2], strsplit(spl1[3], split = "-")[[1]][1]), 8, 1, "sc")
  }
}), c("cell_ID", "timepoint", "replicate_id", "treatment"))

rownames(phenoData) <- phenoData$cell_ID
phenoData$cell_ID <- NULL
colnames(assayData) <- rownames(phenoData)

# make dataframe for gfp counts
gfpData <- assayData["GFP",]

# remove gfp counts from assayData
assayData <- assayData[-which(rownames(assayData) == "GFP"),]

write.table(x=gfpData, file=paste0(output_path, 'gfpData.csv'), sep='\t', row.names=TRUE, quote=FALSE, col.names=NA)
write.table(x=assayData, file=paste0(output_path, 'assayData.csv'), sep='\t', row.names=TRUE, quote=FALSE, col.names=NA)
write.table(x=phenoData, file=paste0(output_path, 'phenoData.csv'), sep='\t', row.names=TRUE, quote=FALSE, col.names=NA)
```

</details>

---

</br>

**_Run RNA velocity_**<a name="velocyto"></a>

The RNA velocity sub-workflow takes the output BAM files from HISAT2 and groups all the cells into a single channel input. We then run Velocyto which generates a single Loom file output containing spliced, unspliced and spanning reads.

<details open><summary class="simple">Collapse sub-workflow</summary>
<p>

```groovy
// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {velocyto_run_smartseq2} from "$baseDir/../luslab-nf-modules/tools/velocyto/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow velocyto_smartseq2 {
    take:
        bam_files
        gtf

    main:
        // group bams into a single channel for velocyto
        ch_velocyto_bam = bam_files
            .map { [[sample_id:"all_cells"], file(it[1], checkIfExists: true)] }
            .groupTuple(by: 0)

        velocyto_run_smartseq2 ( params.modules['velocyto_run_smartseq2'], ch_velocyto_bam, gtf )

    emit:
        velocyto_counts = velocyto_run_smartseq2.out.velocyto
}
```

</details>

</br>
