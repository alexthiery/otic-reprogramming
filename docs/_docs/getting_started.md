---
layout: page
title: Getting Started
order: 6
---

## Data availability

This repository contains the required code to run the entire alignment and downstream analysis pipeline. For those who would like to re-run the analysis, the following files should be downloaded:

- to run the analysis including the alignment, the raw fastq sequencing files can be found [here]().
- to run the downstream analysis from the UMI counts generated from 10x Genomics Cell Ranger are embedded can be found within the repository [here]("./alignmentOut/cellrangerCounts_renamed").

## Analysis pre-requisites

The pipeline is run using Nextflow and Docker to ensure reproducibility. The repository can be called directly from GitHub, so to re-run the analysis you just need to install [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) and [Docker](https://docs.docker.com/get-docker/).

If you are wanting to run the downstream analysis interactively outside of Nextflow, you still need to download Docker and you will also need to download this repository.

## Docker

Docker allows us to run our analysis through a container with all required libraries and dependencies. This ensures cross-platform reproducibility of our results.

The docker image used can be found [here](https://hub.docker.com/r/alexthiery/10x_neural_tube)

We have also included Rstudio within our docker image to allow you to run the downstream analysis interactively if you wish - for details on how to do this click [here](#interactive-downstream-analysis).

## Nextflow

You can easily re-run our entire pipeline in Nextflow using the following steps:

1. Install Nextflow (version >=20.07.1) and Docker
2. Download chick genome ([galgal5](ftp://ftp.ensembl.org/pub/release-94/fasta/gallus_gallus/dna/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa.gz))
3. For this analysis we used a manually edited GTF file as described in the paper methods. However you can also run the analysis with the GTF from ensembl ([galgal5](ftp://ftp.ensembl.org/pub/release-94/gtf/gallus_gallus/Gallus_gallus.Gallus_gallus-5.0.94.gtf.gz))
4. Download raw reads from [here]()
5. Make a sampleInfo.csv file containing the sample names and corresponding paths using this [template](sampleInfo.csv)
6. Run Nextflow using the following command

```sh
nextflow run alexthiery/10x_neural_tube \
-profile docker \
--metadata <path to sampleInfo.csv> \
--gtf <path to gtf> \
--fa <path to genome>
```

This pipeline is configured to be ran on a cluster with 224GB memory and 32CPUs by default. The `-profile` flag can be used to set either 'docker' or 'singularity', depending on the container application installed on your system. These settings can be adjusted by replacing the `-profile` flag with a custom config file as below.

```sh
nextflow run alexthiery/10x_neural_tube \
-c <path to custom.config file> \
--metadata <path to sampleInfo.csv> \
--gtf <path to gtf> \
--fa <path to genome>
```

For a template of a custom.config file, see [crick.config](conf/crick.config). Further information on Nextflow config files can be found [here](https://www.nextflow.io/docs/latest/config.html#configuration-file).

## Interactive downstream analysis

If do not want to re-run the alignment, but would like to run the downstream analysis from the count files, you can run RStudio from within the docker container.

To do this, follow these steps:

1. clone this repository to your local computer
2. start a terminal session and download the docker image - `docker pull alexthiery/10x_neural_tube:v1.0`
3. within terminal launch a docker container interactively - `docker run --rm -e PASSWORD=test -p 8787:8787 -v <PATH_TO_LOCAL_REPOSITORY>:/home/rstudio alexthiery/10x_neural_tube:v1.0`
4. go to `localhost:8787` in the browser to open RStudio
