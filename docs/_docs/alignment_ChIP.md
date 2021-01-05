---
layout: page
title: ChIPseq alignment
order: 3
---

<br/>

Our ChIPseq data was aligned using the nf-core ChIPseq [v1.2.0](https://nf-co.re/ChIPseq/1.2.0/usage) nextflow pipeline.

The paths to the sample fastq files are provided via a samplesheet csv. The template csv used to run this analysis can be found [here](https://github.com/alexthiery/otic-reprogramming/blob/master/NF-ChIP_alignment/crick_samplesheet.csv).

To run this alignment you will also need to generate a custom configuration file. Details can be found [here]({{site.baseurl}}{% link _docs/nextflow_configuration.md %})

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
