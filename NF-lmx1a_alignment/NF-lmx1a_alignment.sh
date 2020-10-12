#!/bin/sh
#SBATCH --job-name=NF-lmx1a
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.01.0
ml Singularity/2.6.0-foss-2016b
ml Graphviz

## UPDATE PIPLINE
nextflow pull nf-core/rnaseq

## RUN alignment
nextflow run nf-core/rnaseq \
  -r 1.4.2 \
  -profile crick \
  -with-singularity /camp/apps/misc/stp/babs/nf-core/singularity/rnaseq/1.4.2/nfcore-rnaseq-1.4.2.img \
  --reads "/camp/home/thierya/working/raw_data/rnaoe/*_{1,2}.fastq.gz" \
  --fasta "/camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa" \
  --gtf "/camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf" \
  --outdir results-lmx1a \
  --email alex.thiery@crick.ac.uk