---
layout: page
label: Home
order: 1
---

The vertebrate inner ear is an exclusive vertebrate feature, key to the perception of sound and movement that arises from a pool of ectodermal progenitors. In our study (Buzzi et.al., 2021) we used a multi-omics approach to explore the molecular mechanisms that regulate the commitment of these progenitors and foster the acquisition of otic character.

In this website you will find the documentation of our High-Throughput-Sequencing data analysis, from sequence alignments to downstream analysis.

---

### Nextflow

</br>

Our analysis uses a combination of NF-core (ChIP, ATAC and bulk RNAseq) pipelines, as well as a custom SmartSeq2 Nextflow pipeline to align our data. We have also developed a custom downstream analysis pipeline in order to integrate the outputs from our datasets into a single reproducible workflow.

For those who would like to re-run the analysis in Nextflow, including re-aligning any of our data, click [here]({{site.baseurl}}/general/quick_start) for instructions.

---

### Docker

</br>

Our downstream analysis pipeline makes use of a custom Docker container, which also allows for interactive downstream exploration of the data in RStudio server. This allows you to explore the data within the same computing environment which our analysis was carried out in. You can find instructions on how to run the Docker container for data exploration [here]({{site.baseurl}}/downstream/downstream_intro#interactive).
