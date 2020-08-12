#!/bin/bash
#SBATCH --job-name=otic-reprogamming
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

export TERM=xterm

# # command to symlinking samples ready for test alignment
# mkdir -p /camp/home/thierya/working/raw_data/otic-reprogamming/scRNAseq/testData/ss8
# mkdir -p /camp/home/thierya/working/raw_data/otic-reprogamming/scRNAseq/testData/ss11
# mkdir -p /camp/home/thierya/working/raw_data/otic-reprogamming/scRNAseq/testData/ss15
# for file in $(find /camp/home/thierya/working/raw_data/ailin_scRNAseq/Samples/*/Files/*1234.fastq.gz | head -2); do ln -s $file /camp/home/thierya/working/raw_data/otic-reprogamming/scRNAseq/testData/ss8 ; done
# for file in $(find /camp/home/thierya/working/raw_data/ailin_scRNAseq/ss11_123fq/* | head -2); do ln -s $file /camp/home/thierya/working/raw_data/otic-reprogamming/scRNAseq/testData/ss11 ; done
# for file in $(find /camp/home/thierya/working/raw_data/ailin_scRNAseq/ss15_123fq/* | head -2); do ln -s $file /camp/home/thierya/working/raw_data/otic-reprogamming/scRNAseq/testData/ss15 ; done

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
-profile crick_test \
-resume