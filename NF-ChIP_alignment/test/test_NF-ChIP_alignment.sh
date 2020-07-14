#!/bin/sh

export TERM=xterm

## LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.01.0
ml Singularity/3.4.2
ml Graphviz

## UPDATE PIPLINE
nextflow pull nf-core/chipseq

# RUN alignment
nextflow run nf-core/chipseq \
-r 1.2.0 \
-c /camp/home/thierya/scratch/otic-reprogramming/NF-ChIP_alignment/conf/crick.config \
--input /camp/home/thierya/scratch/otic-reprogramming/NF-ChIP_alignment/ChIP_design.csv \
--fasta /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa \
--gtf /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf \
--bwa_index /camp/home/thierya/scratch/otic-reprogramming/NF-ChIP_alignment/genome/BWAIndex/Gallus_gallus.GRCg6a.dna.toplevel.fa \
--macs_gsize 1.05e9 \
--skip_diff_analysis \
--outdir ./results \
--email alex.thiery@crick.ac.uk \
-resume

# nextflow run nf-core/chipseq \
# -r 1.2.0 \
# -c /camp/home/thierya/scratch/otic-reprogramming/NF-ChIP_alignment/conf/crick.config \
# --input /camp/home/thierya/scratch/otic-reprogramming/NF-ChIP_alignment/test/ChIP_design_test.csv \
# --fasta /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa \
# --gtf /camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf \
# --bwa_index /camp/home/thierya/scratch/otic-reprogramming/NF-ChIP_alignment/genome/BWAIndex/Gallus_gallus.GRCg6a.dna.toplevel.fa \
# -resume


nextflow run nf-core/chipseq \
-r 1.2.0 \
-profile docker \
--input /Users/alex/dev/repos/otic-reprogramming/NF-ChIP_alignment/test/ChIP_design_test.csv \
--fasta /Users/alex/dev/genomes/galgal5/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa \
--gtf /Users/alex/dev/genomes/galgal5/Gallus_gallus.Gallus_gallus-5.0.94_modified.gtf \
--bwa_index /Users/alex/dev/repos/otic-reprogramming/NF-ChIP_alignment/genome/BWAIndex/Gallus_gallus.GRCg6a.dna.toplevel.fa \
-c /Users/alex/dev/repos/otic-reprogramming/NF-ChIP_alignment/conf/local.config \
--macs_gsize 1.05e9 \
--skip_diff_analysis \
-resume




