# Otic Reprogramming

<br/>

The vertebrate inner ear is an exclusive vertebrate feature, key to the perception of sound and movement that arises from a pool of ectodermal progenitors. In our study (Buzzi et.al., 2022) we used a multi-omics approach to explore the molecular mechanisms that regulate the commitment of these progenitors and foster the acquisition of otic character.

This repository contains all of the code required to re-run our analysis.

Full documentation of our High-Throughput-Sequencing data analysis, from sequence alignments to downstream analysis, can be found [here](https://alexthiery.github.io/otic-reprogramming/).

---

<br/>

### Clone repository

<br/>

In order to reproduce our analysis, you will first need to download our otic-repregramming GitHub repository. To do this run:

```shell
git clone https://github.com/alexthiery/otic-reprogramming
```

---

<br/>

### Access read counts and peaks

<br/>

RNAseq read counts (Lmx1aE1 vs Sox3U3, SOX8OE and SmartSeq2) and ChIP/ATAC peak files are embedded in the [repository](./alignment_output) - they will be downloaded automatically when the repository is cloned.

---

</br>

### Genome browser tracks

</br>

Genome browser tracks for out ATAC and ChIPseq data can be found [here](https://genome.ucsc.edu/s/alexthiery/Sox8_otic_reprogramming).

---

</br>

### qPCR

<br/>

The script used to run our qPCR analyis can be found [here](./qPCR/qpcr.R). The docker image used for the multi-omic analysis also contains the required packages to re-run the qPCR analysis `alexthiery/otic-reprogramming-r_analysis:latest`
