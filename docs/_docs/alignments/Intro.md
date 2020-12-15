---
title: Intro
category: Rosalind
order: 1
---

## Introduction to Rosalind

Rosalind is our high performance cluster at KCL.

They have recently re-written the [HPC wiki](https://rosalind.kcl.ac.uk/) which is very informative. This is the first place to look if you are having problems.

- I'm stuck! - Email [IT](mailto:rosalind-support@kcl.ac.uk) with a description of your problem. Make sure you include as much information as possible, including job submission script details, command line errors, paths to files involved.

### Getting access

To request access email [IT](mailto:rosalind-support@kcl.ac.uk), you will need to generate a public key to send to the HPC team. 

Here are some instructions on how to do this:

Open a terminal window. At the shell prompt, type the following command:

`ssh-keygen -t rsa`

The ssh-keygen program will prompt you for the location of the key file. Press Return to accept the defaults. You can optionally specify a passphrase to protect your key material. Press Return to omit the passphrase. The output of the program will look similar to this:

```
Enter file in which to save the key (/Users/alex/.ssh/id_rsa):
Created directory '/Users/alex/.ssh'.
Enter passphrase (empty for no passphrase):
Enter same passphrase again:
Your identification has been saved in /Users/alex/.ssh/id_rsa.
Your public key has been saved in /Users/alex/.ssh/id_rsa.pub.
```

### Once you have an account, how to log on

You can log on to Rosalind via ssh:

`ssh <username>@login.rosalind.kcl.ac.uk`

Once you logon you'll find yourself in your home directory `/users/<username>`.

You will have been given a personal scratch directory, which is where you should run most of your analyses for storage purposes `/scratch/users/<username>`.

You should also have access to streit lab scratch storage `/scratch/groups/streit`. This is where we store lab raw sequencing data, genomes, as well as downloaded nextflow pipelines and singularity images. We use this lab scratch in order to reduce the amount of file duplication across our accounts.

It is recommmended to create symbolic links from your personal scratch `ln -s /scratch/users/<username> /users/<username>/scratch` and streit lab scratch `ln -s /scratch/groups/streit/ /users/<username>/streit_scratch` to your home directory. This makes it far easier to move between the different directories on Rosalind (`cd ~/scratch` for personal scratch, `cd ~/streit_scratch` for streit scratch, and `cd ~/` for home).

After doing this I would recommend also creating the symbolic links for placing nextflow and singularity files in the joint lab scratch - this will save major headaches later on!

First let's make the singularity symlink: `ln -s /scratch/groups/streit/singularity_cache /users/<username>/.singularity`
Next, symlink for nextflow pipelines: `ln -s /scratch/groups/streit/nextflow /users/<username>/.nextflow`

When you pull nextflow pipelines from github, they will now automatically be placed in the lab scratch!

### Directory structure on Rosalind

On Rosalind you will have access to three main directories: scratch, lab scratch and users. The lab scratch folder contains data that we all have access to and so we like to keep this tidy. We ask that whenever someone finds themselves downloading a resource file, such as a genome, an annotation .etc that they try and add this to the shared folder so that we can save duplicating large files in our personal directories.

The structure of the shared folder is:

```
/scratch/groups/streit
├── nextflow
|   └── this will be where nextflow pipelines are pulled to from github
├── singularity_cache
|   └── here we store all lab singularity images
├── raw_sequencing_data
|   └── here we store all raw sequencing data
├── ref
    └── genomes
        └── chick
        └── mouse
        └── etc
```

### Doing analysis and running jobs

The "quickstart" version is that Rosalind is divided into "partitions" with certain nodes dedicated for certain things.
KCL members have access to the grc "cpu" and "gpu" partitions.
Mostly you will submit batch jobs to the "cpu" partition.
Never run anything substantial (ie. more than simple head, cd, mkdir .etc commands) on the login node (this is where you are when you first ssh into the cluster).

Note that "threads" are equivalent to --cpus-per-task and the max seems to be 32 on the cpu partition.
We can only have *one* interactive session running at a time. If you can't get one, it might be that you have one running that you forgot about.

To check what you've got running:
`squeue -u <username>`

To cancel a job:
`scancel <JOBID>`

To cancel all user jobs:
`scancel -u <username>`