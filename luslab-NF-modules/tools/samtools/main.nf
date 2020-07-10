#!/usr/bin/env nextflow

// Specify DSL2
nextflow.preview.dsl = 2

// Samtools index
process samtools_index {
    publishDir "${params.outdir}/samtools/index",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-samtools:latest'

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("*.bam.bai"), emit: baiFiles
 
    script:

    // Check main args string exists and strip whitespace
    args = ""
    if(params.samtools_index_args && params.samtools_index_args != '') {
        ext_args = params.samtools_index_args
        args += " " + ext_args.trim()
    }

    // Construct CL line
    index_command = "samtools index ${args} -@ ${task.cpus} ${bam[0]}"

    // Log
    if (params.verbose){
        println ("[MODULE] samtools/index command: " + index_command)
    }
    
    """
    ${index_command}
    """
}

// Samtools view - only works for bam files - requires bam and bai
process samtools_view {
    publishDir "${params.outdir}/samtools/view",
        mode: "copy", overwrite: true

    container 'luslab/nf-modules-samtools:latest'

    input:
        tuple val(sample_id), path(bam_bai)

    output:
        tuple val(sample_id), path("*.b*"), emit: bamFiles
 
    script:

    // Check main args string exists and strip whitespace
    args = "-b -h"
    if(params.samtools_view_args && params.samtools_view_args != '') {
        ext_args = params.samtools_view_args
        args += " " + ext_args.trim()
    }

    // Construct CL line
    view_command = "samtools view ${args} -@ ${task.cpus} -o ${bam_bai[0].simpleName}.filt.bam ${bam_bai[0]} ${params.samtools_view_region}"

    // Log
    if (params.verbose){
        println ("[MODULE] samtools/view command: " + view_command)
    }
    
    """
    ${view_command}
    samtools index -@ ${task.cpus} ${bam_bai[0].simpleName}.filt.bam
    """
}
