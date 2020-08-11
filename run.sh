#!/bin/bash
#SBATCH --job-name=otic-reprogamming
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.07.1
ml Singularity/3.4.2
ml Graphviz

nextflow pull alexthiery/otic-reprogramming -r dev

export NXF_VER=20.07.1
export NXF_SINGULARITY_CACHEDIR=/camp/home/thierya/working/NF_singularity

nextflow run alexthiery/otic-reprogramming \
-hub github \
-r dev \
-profile crick \
-resume