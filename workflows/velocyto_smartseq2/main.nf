#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {velocyto_run_smartseq2} from "$baseDir/../luslab-nf-modules/tools/velocyto/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow velocyto_smartseq2 {
    take:
        bam_files
        gtf

    main:
        // group bams into a single channel for velocyto
        ch_velocyto_bam = bam_files
            .map { [[sample_id:"all_cells"], file(it[1], checkIfExists: true)] }
            .groupTuple(by: 0)

        velocyto_run_smartseq2 ( params.modules['velocyto_run_smartseq2'], ch_velocyto_bam, gtf )

    emit:
        velocyto_counts = velocyto_run_smartseq2.out.velocyto
}