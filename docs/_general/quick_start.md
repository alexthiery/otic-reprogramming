---
layout: page
label: Quick start
order: 1
---

## Quick start

<br/>

Reproducing our analysis is computationally intensive, therefore we recommend to run the pipeline on a HPC.

The Nextflow alignment pipelines generate numerous temporary files stored within `work` directories. You can clear these after each alignment without affecting the outputs required for downstream analysis. Without clearing the temporary files, we recommend that you have at minimum ~4TB available storage, 8CPUs and 64GB RAM.

In order to reproduce our analysis, you will need to:

- [download the GitHub project repository](#download_repo)
- [install Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- [install Docker Desktop](https://www.docker.com/get-started) (if running on local PC)
- [download Galgal6 and sequencing data](#download_data)
- [set Nextflow configuration file](#config)
- [align sequencing data](#align)
- [run downstream analysis](#downstream)

**Important!** All shell scripts are to be executed from the base directory of the project!

---

<br/>

## Download GitHub repository<a name="download_repo"></a>

<br/>

In order to reproduce our analysis, you will first need to download our otic-repregramming GitHub repository. To do this run:

```shell
git clone --recurse-submodules https://github.com/alexthiery/otic-reprogramming
```

---

<br/>

## Download data<a name="download_data"></a>

<br/>

To download the Galgal6 genome from Ensembl, run:

```shell
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-97/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz <LOCAL_PATH>
```

To download the Galgal6 gtf from Ensembl, run:

```shell
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-97/gtf/gallus_gallus/Gallus_gallus.GRCg6a.97.gtf.gz <LOCAL_PATH>
```

Once the genome files have been downloaded, they need to be unzipped before running the alignment.

**ADD SEQUENCING DATA DOWNLOAD PATHS**

---

<br/>

## Nextflow configuration file<a name="config"></a>

<br/>

Configuration properties and file paths are passed to Nextflow via configuration files.

In order to re-run the alignments and downstream analysis for this project, you will need to create a custom configuration file which contains:

- genome file path
- gtf file path
- max system resource allocation
- Docker/Singularity activation

<br/>

Below is a copy of the config file used to run the pipeline at The Francis Crick Institute HPC. Given that this alignment was run on a HPC, we configure Nextflow to run with Singularity instead of Docker.

```bash
#!/usr/bin/env nextflow

singularity {
  enabled = true
  autoMounts = true
  docker.enabled = false
}

process {
  executor = 'slurm'
}

params {
  max_memory = 224.GB
  max_cpus = 32
  max_time = 72.h

  fasta = "/camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa"
  gtf = "/camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf"
}
```

Further information on Nextflow configuration can be found [here](https://www.nextflow.io/docs/latest/config.html).

After you have saved your custom configuration file, you are ready to align the data!

---

<br/>

## Sequence alignments<a name="align"></a>

<br/>

When running the NF-core alignments, the pipeline may be interupted as the relevant docker images are downloaded. If this happens, simply re-run the shell script and the pipeline will continue from where it ended.

<br/>

### ChIPseq alignment

<br/>

Our ChIPseq data was aligned using the NF-core ChIPseq [v1.2.0](https://nf-co.re/ChIPseq/1.2.0/usage) pipeline.

Sample FASTQ paths are provided via a samplesheet csv. The template csv used to run this analysis can be found [here](https://github.com/alexthiery/otic-reprogramming/blob/master/NF-ChIP_alignment/crick_samplesheet.csv).

After saving the required files, you can run the alignment with the following shell script.

```shell
export NXF_VER=20.07.1

nextflow pull nf-core/chipseq

nextflow run nf-core/chipseq \
    -r 1.2.0 \
    -c <PATH_TO_CONFIG> \
    --input <PATH_TO_SAMPLESHEET.CSV> \
    --macs_gsize 1.05e9 \
    --skip_diff_analysis \
    --outdir output/NF-ChIP_alignment \
    --email <INSERT_EMAIL_ADDRESS> \
    -resume
```

---

<br/>

### ATACseq alignment

<br/>

Our ATACseq data was aligned using the NF-core ATACseq [v1.2.0](https://nf-co.re/ATACseq/1.2.0/usage) pipeline.

Sample FASTQ paths are provided via a samplesheet csv. The template csv used to run this analysis can be found [here](https://github.com/alexthiery/otic-reprogramming/blob/master/NF-ATAC_alignment/crick_samplesheet.csv).

After saving the required files, you can run the alignment with the following shell script.

```shell
export NXF_VER=20.07.1

nextflow pull nf-core/atacseq

nextflow run nf-core/atacseq \
    -r 1.2.0 \
    -c <PATH_TO_CONFIG> \
    --input <PATH_TO_SAMPLESHEET.CSV> \
    --macs_gsize 1.05e9 \
    --narrow_peak \
    --skip_diff_analysis \
    --outdir output/NF-ATAC_alignment \
    --email <INSERT_EMAIL_ADDRESS> \
    -resume
```

---

<br/>

### Bulk RNAseq alignment

<br/>

Our bulk RNAseq data were aligned using the NF-core RNAseq [v2.0](https://nf-co.re/RNAseq/2.0/usage) pipeline.

We need to align two separate bulk RNAseq datasets.

<br/>

**_Lmx1aE1 vs Sox3U3_**

First, we collected cells from putative otic and epibranchial fates using Lmx1aE1 and Sox3U3 enhancers respectively.

An example samplesheet csv can be found [here](https://github.com/alexthiery/otic-reprogramming/blob/master/NF-lmx1a_alignment/crick_samplesheet.csv).

Execute the following shell script to run the Lmx1a alignment.

```shell
export NXF_VER=20.07.1

nextflow pull nf-core/rnaseq

nextflow run nf-core/rnaseq \
  -r 2.0 \
  -c <PATH_TO_CONFIG> \
  --input <PATH_TO_LMX1A_SAMPLESHEET.CSV>  \
  --outdir output/NF-lmx1a_alignment \
  --email <INSERT_EMAIL_ADDRESS> \
  -resume
```

<br/>

**_Sox8 overexpression_**

Next, we collected cells following Sox8 overexpression and compared relative to control cells.

Sox8 overexpression samples are double positive Sox8-mCherry/Lmx1a-EGFP, whilst control samples are labelled with constitutive mCherry/EGFP.

An example samplesheet csv can be found [here](https://github.com/alexthiery/otic-reprogramming/blob/master/NF-sox8_alignment/crick_samplesheet.csv)

Execute the following shell script to run the Sox8 overexpression alignment.

```shell
export NXF_VER=20.07.1

nextflow pull nf-core/rnaseq

nextflow run nf-core/rnaseq \
  -r 2.0 \
  -c <PATH_TO_CONFIG> \
  --input <PATH_TO_SOX8_SAMPLESHEET.CSV>  \
  --outdir output/NF-sox8_alignment \
  --email <INSERT_EMAIL_ADDRESS> \
  -resume
```

---

<br/>

### Smartseq2 single cell RNAseq alignment<a name="smartseq"></a>

<br/>

Our Smartseq2 single cell RNAseq data was aligned using a custom DSL2 Nextflow pipeline. Details of the pipeline processes can be found [here]({{ site.baseurl }}{% link _general/smartseq2_alignment.md %}).

As with the NF-core pipelines above, our custom SmartSeq2 alignment pipeline passes the sample FASTQ files via a samplesheet csv. An example samplesheet can be found [here](https://github.com/alexthiery/otic-reprogramming/blob/master/NF-smartseq2_alignment/crick_samplesheet.csv).

As we have individual FASTQ files for each cell, it would be cumbersome to list each file individually. Instead we can pass an entire sample FASTQ directory which is then enumerated in Nextflow. For this to work, only the paths in column 2 ($data1) in the samplesheet csv should be edited.

Once you have saved your samplesheet, execute the following shell script to run the SmartSeq2 alignment.

```shell
export NXF_VER=20.07.1

nextflow run ./NF-smartseq2_alignment/main.nf \
  -c <PATH_TO_CONFIG> \
  --input <PATH_TO_SMARTSEQ2_SAMPLESHEET.CSV> \
  --outdir output/NF-smartseq2_alignment \
  -resume
```

---

<br/>

### Downstream analysis<a name="downstream"></a>

<br/>

We have built a custom Nextflow pipeline to integrate the outputs from the separate alignments.

After aligning the data, you will need to create a samplesheet csv containing the paths to all of the alignment output folders. You can find an example of this [here](https://github.com/alexthiery/otic-reprogramming/blob/master/NF-downstream_analysis/crick_samplesheet.csv).

To run the integrated downstream analysis execute the following shell script.

```shell
export NXF_VER=20.07.1

nextflow run ./NF-downstream_analysis/main.nf \
  -c <PATH_TO_CONFIG> \
  --input <PATH_TO_DOWNSTREAM_ANALYSIS_SAMPLESHEET.CSV> \
  --outdir output/NF-downstream_analysis \
  -resume
```

If you would like to interactively carry out downstream anlysis click [here]({{site.baseurl}}/downstream/downstream_intro).
