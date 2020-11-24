#!/bin/sh

nextflow run ./NF-downstream_analysis/main.nf \
  -c ./configs/local.config \
  --input ./NF-downstream_analysis/local_samplesheet.csv \
  --outdir output/NF-downstream_analysis \
  -resume