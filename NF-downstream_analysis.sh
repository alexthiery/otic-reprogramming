#!/bin/sh
#SBATCH --job-name=NF-downstream_analysis
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/21.10.6
ml Singularity/3.6.4
ml Graphviz

export TERM=xterm
export NXF_VER=21.10.6
export NXF_SINGULARITY_CACHEDIR=/camp/home/thierya/working/NF_singularity

nextflow run ./NF-downstream_analysis/main.nf \
  -c ./configs/crick.config \
  --input ./NF-downstream_analysis/crick_samplesheet.csv \
  --outdir output/NF-downstream_analysis \
  -resume