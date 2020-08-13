#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// include {smartseq2_align} from "$baseDir/workflows/scRNAseq_alignment/main.nf"

// Channel
//     .value(file(params.genome))
//     .set {ch_genome}

// Channel
//     .value(file(params.gtf))
//     .set {ch_gtf}

// workflow {
//     smartseq2_align (ch_genome, ch_gtf, params.smartseq2_sample_csv)
//     smartseq2_align.out.velocyto_counts | view
//     smartseq2_align.out.merged_counts | view
// }





include {peak_intersect} from "$baseDir/workflows/peak_intersect/main.nf"

// (ch_testData_ChIP.filter { it[0] == 'K27Ac' }.subscribe { println it }, ch_testData_ATAC.filter { it[0] == 'Sample1' }.subscribe { println it })

workflow {
    peak_intersect (params.peak_intersect_sample_csv)
}




// /*------------------------------------------------------------------------------------*/
// /* Define params
// --------------------------------------------------------------------------------------*/

// params.bedtools_intersect_args = '-wa -wb'
// params.verbose = true

// /*------------------------------------------------------------------------------------*/
// /* Module inclusions 
// --------------------------------------------------------------------------------------*/


// // Define test data input channel
// Channel
//     .from(testData)
//     .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
//     .set {ch_testData}

// ch_testData
//     .filter { it[0] == 'K27Ac' }
//     .set {ch_K27Ac}

// ch_testData
//     .filter { it[0] == 'K27me3' }
//     .set {ch_K27me3}

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