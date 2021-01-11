---
layout: page
label: Intro
category: Downstream analysis
order: 1
---

## Downstream analysis

</br>

In order to maintain reproducibility when carrying out the downstream analysis, we have developed a custom Nextflow pipeline which integrates the alignment outputs into a single Nextflow workflow.

This workflow makes use of a single custom Docker container which contains all required packages which can be found [here](https://hub.docker.com/repository/docker/alexthiery/otic-reprogramming-r_analysis).

If you would like to re-run the entire downstream analysis, you will need to first align the data by following the steps in the [**_Quick start_**]({{site.baseurl}}/general/quick_start) guide.

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

## Interactive downstream analysis

</br>

If do not want to re-run the alignment, but would like to run the downstream analysis from the count files, you can run RStudio from within the Docker container. This will ensure that you have all of the same packages and dependencies required to carry out the analysis.

If you would like to interactively explore the data, follow these steps:

1. download the relevant alignment output
2. [download](https://github.com/alexthiery/otic-reprogramming/archive/master.zip) our GitHub repository to your local computer
3. start a terminal session and pull the Docker image from Dockerhub - `docker pull alexthiery/otic-reprogramming-r_analysis:latest`
4. within terminal launch a Docker container interactively - `docker run --rm -e PASSWORD=password -p 8787:8787 -v <PATH_TO_LOCAL_REPOSITORY>:/home/rstudio alexthiery/otic-reprogramming-r_analysis:latest`
5. go to `localhost:8787` in the browser to open RStudio
6. enter the username `rstudio` and the password `password`
7. change the paths at the beginning of the script you wish to analyse in order to read in the data

---

</br>

## Further documentation

</br>

- [Lmx1a differential expression analysis]({{site.baseurl}}/downstream/lmx1a_downstream)
- [SmartSeq2 scRNAseq analysis]({{site.baseurl}}/downstream/smartseq2_downstream)
- [Sox8 differential expression analysis]({{site.baseurl}}/downstream/sox8_downstream)
- [Enhancer discovery pipeline]({{site.baseurl}}/downstream/enhancer_discovery)