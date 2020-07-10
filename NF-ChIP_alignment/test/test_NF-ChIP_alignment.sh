#!/bin/bash
#SBATCH --job-name=10x
#SBATCH -t 72:00:00
#SBATCH --mail-type=ALL,ARRAY_TASKS
#SBATCH --mail-user=alex.thiery@crick.ac.uk

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.01.0
ml Singularity/3.4.2
ml Graphviz

nextflow run nf-core/chipseq \
-r 1.2.0 \
-profile crick \
--input /camp/home/thierya/scratch/otic-reprogramming/NF-ChIP_alignment/test/ChIP_design_test.csv \
--fasta /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa \
--gtf /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf \
-c /camp/home/thierya/scratch/otic-reprogramming/NF-ChIP_alignment/test/local.config \
--save_reference \
-resume

# nextflow run nf-core/chipseq \
# -r 1.2.0 \
# -profile docker \
# --input /Users/alex/dev/repos/otic-reprogramming/NF-ChIP_alignment/test/ChIP_design_test.csv \
# --fasta /Users/alex/dev/genomes/galgal5/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa \
# --gtf /Users/alex/dev/genomes/galgal5/Gallus_gallus.Gallus_gallus-5.0.94_modified.gtf \
# -c /Users/alex/dev/repos/otic-reprogramming/NF-ChIP_alignment/test/local.config \
# --save_reference \
# -resume



