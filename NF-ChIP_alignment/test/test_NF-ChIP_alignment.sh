#!/bin/bash

nextflow run nf-core/chipseq \
-r 1.2.0 \
-profile docker \
--input /Users/alex/dev/repos/otic-reprogramming/NF-ChIP_alignment/test/ChIP_design_test.csv \
--fasta /Users/alex/dev/genomes/galgal5/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa \
--gtf /Users/alex/dev/genomes/galgal5/Gallus_gallus.Gallus_gallus-5.0.94_modified.gtf \
-c /Users/alex/dev/repos/otic-reprogramming/NF-ChIP_alignment/test/local.config \
--save_reference \
-resume
