#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// include {smartseq2_align} from "$baseDir/workflows/scRNAseq_alignment/main.nf"

// Channel
//     .value(file(params.genome, checkIfExists: true))
//     .set {ch_genome}

// Channel
//     .value(file(params.gtf, checkIfExists: true))
//     .set {ch_gtf}

// workflow {
//     smartseq2_align (ch_genome, ch_gtf, params.smartseq2_sample_csv)
//     smartseq2_align.out.velocyto_counts | view
//     smartseq2_align.out.merged_counts | view
// }

// include {merge_counts} from "$baseDir/custom-nf-modules/rscript/main.nf"





// include {peak_intersect} from "$baseDir/workflows/peak_intersect/main.nf"

// workflow {
//     peak_intersect (params.peak_intersect_sample_csv, ch_genome, ch_gtf)
// }



// /*------------------------------------------------------------------------------------*/
// /* Run tests
// --------------------------------------------------------------------------------------*/

// // Run workflow
// workflow {
//     // intersect peak files
//     bedtools_intersect(ch_K27Ac, ch_atacData)

//     bedtools_subtract(bedtools_intersect.out, ch_K27me3)
    
//     // View outputs
//     bedtools_intersect.out | view
// }





include {merge_counts} from "/Users/alex/dev/repos/otic-reprogramming/rscript/main.nf"
// workflow for just merging the output of htseq count for testing

Channel
    .fromPath("/Users/alex/dev/repos/otic-reprogramming/output/htseq_count/*txt")
    .map { [[sample_id:"all_cells"], file(it, checkIfExists: true)] }
    .groupTuple(by: 0)
    .set {ch_all_counts}


workflow {
    merge_counts (params.modules['merge_counts'], ch_all_counts)
}