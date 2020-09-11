#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include {add_gfp} from "$baseDir/workflows/add_gfp/main.nf"
include {smartseq2_align} from "$baseDir/workflows/scRNAseq_alignment/main.nf"

Channel
    .value(file(params.genome, checkIfExists: true))
    .set {ch_genome}

Channel
    .value(file(params.gtf, checkIfExists: true))
    .set {ch_gtf}

Channel
    .value(file(params.gfp_seq, checkIfExists: true))
    .set {ch_gfp_seq}

workflow {
    // add gfp to genome and gtf
    add_gfp (ch_genome, ch_gtf, ch_gfp_seq)
    // align using smartseq2 workflow
    // smartseq2_align (add_gfp.out.genome, add_gfp.out.gtf, params.smartseq2_sample_csv)
    // smartseq2_align.out.velocyto_counts | view
    // smartseq2_align.out.merged_counts | view
}




/*------------------------------------------------------------------------------------*/
/* Workflow to run peaks intersect
--------------------------------------------------------------------------------------*/

// include {peak_intersect} from "$baseDir/workflows/peak_intersect/main.nf"


// workflow {
//     peak_intersect (params.peak_intersect_sample_csv, ch_genome, ch_gtf)
// }
