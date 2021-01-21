---
layout: page
label: Intro
category: Downstream analysis
order: 1
---

## Downstream analysis

</br>

In order to maintain reproducibility when carrying out the downstream analysis, we have developed a custom Nextflow pipeline which integrates the alignment outputs into a single workflow.

This workflow makes use of a [custom Docker container](https://hub.docker.com/repository/docker/alexthiery/otic-reprogramming-r_analysis) which contains all required packages.

If you would like to re-run the entire downstream analysis, you will need to first align the data by following the steps in the [_Quick start_]({{site.baseurl}}/general/quick_start) guide.

After aligning the data, you will need to create a samplesheet csv containing the paths to all of the alignment output folders. You can find an example of this [here](https://github.com/alexthiery/otic-reprogramming/blob/master/NF-downstream_analysis/crick_samplesheet.csv).

The entire downstream analysis can then be carried out by executing using the following shell script.

```shell
export NXF_VER=20.07.1

nextflow run ./NF-downstream_analysis/main.nf \
  -c <PATH_TO_CONFIG> \
  --input <PATH_TO_DOWNSTREAM_ANALYSIS_SAMPLESHEET.CSV> \
  --outdir output/NF-downstream_analysis \
  -resume
```

---

</br>

## Interactive downstream analysis<a name="interactive"></a>

</br>

If do not want to re-run the alignment, but would like to run the downstream analysis from the count files, you can run RStudio server from within the Docker container. This will ensure that you have all of the same packages and dependencies required to carry out the analysis.

RNAseq read counts (Lmx1aE1 vs Sox3U3, SOX8OE and SmartSeq2) and ChIP/ATAC peaks files are embedded in the [repository](https://github.com/alexthiery/otic-reprogramming/tree/master/alignment_output) - they will be downloaded automatically when the repository is cloned.

To interactively explore the data, follow these steps:

1. clone our GitHub repository to your local computer - `git clone --recurse-submodules https://github.com/alexthiery/otic-reprogramming`
2. start a terminal session and pull the Docker image from Dockerhub - `docker pull alexthiery/otic-reprogramming-r_analysis:latest`
3. within terminal launch a Docker container interactively - `docker run --rm -e PASSWORD=password -p 8787:8787 -v <PATH_TO_LOCAL_REPOSITORY>:/home/rstudio alexthiery/otic-reprogramming-r_analysis:latest`
4. go to `localhost:8787` in the browser to open RStudio
5. enter the username `rstudio` and the password `password`
6. access the desired R script in the `Files` tab in R studio by following the path `./NF-downstream_analysis/bin/<FILE_OF_INTEREST>`

---

</br>

## Further documentation

</br>

- [Lmx1aE1 vs Sox3U3 differential expression analysis]({{site.baseurl}}/downstream/lmx1a_downstream)
- [SmartSeq2 scRNAseq analysis]({{site.baseurl}}/downstream/smartseq2_downstream)
- [Sox8 over-expression differential expression analysis]({{site.baseurl}}/downstream/sox8_downstream)
- [Enhancer discovery pipeline]({{site.baseurl}}/downstream/enhancer_discovery)
