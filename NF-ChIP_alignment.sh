#!/bin/sh
#SBATCH --job-name=NF-chip
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
nextflow pull nf-core/chipseq

## RUN alignment
nextflow run nf-core/chipseq \
    -r 1.2.0 \
    -c ./configs/crick.config \
    --input ./NF-ChIP_alignment/crick_samplesheet.csv \
    --macs_gsize 1.05e9 \
    --skip_diff_analysis \
    --outdir output/NF-ChIP_alignment \
    --email alex.thiery@crick.ac.uk \
    -resume