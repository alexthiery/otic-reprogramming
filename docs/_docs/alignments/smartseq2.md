---
title: How to run Nextflow on Rosalind
category: Rosalind
order: 3
---

## Installing Nextflow on Rosalind
Nextflow is not one of the avaiable programmes already installed on Rosalind. It has to be installed in a conda environment. 

First load the anaconda module so you can use it:
`module load devtools/anaconda/2019.3-python3.7.3`
This creates a conda environment in the base directory. As you don't have permission to install anything in that directory you need to make a conda environment in your own directory. Go to your home directory and enter:
`conda create --name nextflow`
`conda environments` should return the base directory environment and your new one in your folder

You then have to restart your Rosalind session, type `exit` and then log back in as before

Run this line to activate your new conda environment:
`conda activate nextflow`

Finally, run these lines to set up a conda channel in this environment and install nextflow:
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install nextflow
```


## Check you have the latest version of Nextflow - higher than or equal to 20.07.1
You can do it by first checking which versions of Nextflow are installed in your conda environment. To do this first activate your nextflow conda environment `conda activate nextflow` and then check the nextflow version `nextflow -version`.


## Use a bash wrapper to start up a Nextflow pipeline
You need Rosalind to run Nextflow pipelines. It is convenient to use a bash wrapper like the following one to start up a Nextflow pipeline:

```
## LOAD REQUIRED MODULES
module load apps/singularity/3.5.3

## RUN PIPELINE
nextflow run path/to/pipeline_main.nf \
    --design input.csv \
    -resume
```

In this example, we first load the required Singularity module. [Singularity](https://en.wikipedia.org/wiki/Singularity_(software)) is a containerisation system used on computational clusters instead of Docker, because it is more secure than Docker. Docker images are automatically converted into the Singularity format, so you do not need to think about this.


## Run the pipeline
Next, we execute the command `nextflow run path/to/pipeline_main.nf` which starts the pipeline described in the `path/to/pipeline_main.nf` Nextflow script and provide necessary parameters to Nextflow (one dash) and to the pipeline itself (double dash).

You also need to make sure that the [Rosalind-specific config file for Nextflow](https://github.com/Streit-lab/Streit-lab.github.io/blob/master/rosalind.config) is located in the same directory as the pipeline script (`path/to/` in our example above).

To run an _nf-core pipeline_ on Rosalind, you need the same kind of a bash wrapper, but instead of a path to a Nextflow script, you need to provide the name of the pipeline (Nextflow will download the pipeline automatically):

```
## LOAD REQUIRED MODULES
module load apps/singularity/3.5.3

## RUN nf-core PIPELINE
nextflow run nf-core/chipseq \
    --design design_table.csv \
    --genome GRCm38 \
    --skipSpp \
    --skipDiffAnalysis \
    --pvalue 1e-3 \
    --email your.email@kcl.ac.uk \
    -r 1.0.0
    -c <path_to_rosalind_config> \
    -resume
```

In this example, the name of the pipeline is `nf-core/chipseq`, and we run its version `1.0.0`. In the case of the `chipseq` pipeline, `1.0.0` is not the latest version, but when you begin your project, please use the latest versions of the nf-core pipelines.

**When running an nf-core pipeline on Rosalind, please make sure that you provide the Rosalind config to the `-c` option.** 

To find out more about possible Nextflow and pipeline-specific options that you can use when running an nf-core pipeline, see the [nf-core documentation](https://nf-co.re/usage/introduction) and the pipeline-specific documentation.
