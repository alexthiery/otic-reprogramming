#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {add_genome_gfp; add_gtf_gfp} from "$baseDir/../modules/genome-tools/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow add_gfp {
    take:
        genome
        gtf
        gfp_seq

    main:
        add_genome_gfp (params.modules['add_genome_gfp'], genome, gfp_seq)
        add_gtf_gfp (params.modules['add_gtf_gfp'], gtf, gfp_seq)

    emit:
        genome = add_genome_gfp.out
        gtf = add_gtf_gfp.out
}