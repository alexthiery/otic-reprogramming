#!/bin/sh
#SBATCH --job-name=NF-ChIP
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.01.0
ml Singularity/2.6.0-foss-2016b

## UPDATE PIPLINE
nextflow pull nf-core/chipseq

# RUN alignment
nextflow run nf-core/chipseq \
-r 1.2.0 \
-profile crick \
-with-singularity /camp/apps/misc/stp/babs/nf-core/singularity/chipseq/1.2.0/nfcore-chipseq-1.2.0.img \
--input /camp/home/thierya/working/analysis/otic-reprogramming/NF-ChIP_alignment/ChIP_design.csv \
--fasta /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa \
--gtf /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf \
--macs_gsize 1.05e9 \
--skip_diff_analysis \
--outdir results-chip-no-shift \
--email alex.thiery@crick.ac.uk \
-resume