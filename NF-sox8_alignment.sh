#!/bin/sh
#SBATCH --job-name=NF-sox8
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.07.1
ml Singularity/3.4.2
ml Graphviz

export TERM=xterm
export NXF_VER=20.07.1
export NXF_SINGULARITY_CACHEDIR=/camp/home/thierya/working/NF_singularity

## UPDATE PIPLINE
nextflow pull nf-core/rnaseq

## RUN alignment
nextflow run nf-core/rnaseq \
  -r 2.0 \
  -profile crick \
  --input './NF-sox8_alignment/crick_samplesheet.csv' \
  --fasta '/camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa' \
  --gtf '/camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf' \
  --outdir output/NF-sox8_alignment \
  --email alex.thiery@crick.ac.uk