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
  -c ./configs/crick.config \
  --aligner star \
  --input ./NF-sox8_alignment/crick_samplesheet.csv \
  --outdir output/NF-sox8_alignment \
  --email alex.thiery@crick.ac.uk \
  -resume