---
layout: page
title: ATACseq alignment
order: 3
---

<br/>

Our ATACseq data was aligned using the nf-core ATACseq [v1.2.0](https://nf-co.re/atacseq/1.2.0/usage) nextflow pipeline.

The paths to the sample fastq files are provided via a samplesheet csv. The template csv used to run this analysis can be found [here](https://github.com/alexthiery/otic-reprogramming/blob/master/NF-ATAC_alignment/crick_samplesheet.csv).

To run this alignment you will also need to generate a custom configuration file containing:

- genome file path
- gtf file path
- max system resource allocation
- activate Docker/Singularity

<br/>

Below is a copy of the config file used to run the pipeline at the Crick. Given that this alignment was run on a HPC, we configure nextflow to run with Singularity instead of Docker.

```bash
#!/usr/bin/env nextflow

singularity {
  enabled = true
  autoMounts = true
  docker.enabled = false
}

process {
  executor = 'slurm'
}

params {
  max_memory = 224.GB
  max_cpus = 32
  max_time = 72.h

  fasta = "/camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.dna.toplevel.fa"
  gtf = "/camp/home/thierya/working/genomes/galgal6/Gallus_gallus.GRCg6a.97.gtf"
}
```

<br/>

After downloading the required files, you can run the alignment with the following shell [script](https://github.com/alexthiery/otic-reprogramming/blob/master/NF-ATAC_alignment.sh).

```shell
export NXF_VER=20.07.1

nextflow pull nf-core/atacseq

nextflow run nf-core/atacseq \
 -r 1.2.0 \
 -c ./configs/crick.config \
 --input ./NF-ATAC_alignment/crick_samplesheet.csv \
 --macs_gsize 1.05e9 \
 --narrow_peak \
 --skip_diff_analysis \
 --outdir output/NF-ATAC_alignment \
 --email <INSERT_EMAIL_ADDRESS> \
 -resume
```
