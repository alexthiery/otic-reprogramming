#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

// include {add_gfp} from "$baseDir/workflows/add_gfp/main.nf"
// include {smartseq2_align} from "$baseDir/workflows/scRNAseq_alignment/main.nf"
// include {process_counts} from "$baseDir/workflows/process_counts/main.nf"
// include {velocyto_smartseq2} from "$baseDir/workflows/velocyto_smartseq2/main.nf"

// Channel
//     .value(file(params.genome, checkIfExists: true))
//     .set {ch_genome}

// Channel
//     .value(file(params.gtf, checkIfExists: true))
//     .set {ch_gtf}

// Channel
//     .value(file(params.gfp_seq, checkIfExists: true))
//     .set {ch_gfp_seq}

// workflow {
//     // add gfp to genome and gtf
//     add_gfp (ch_genome, ch_gtf, ch_gfp_seq)

//     // align using smartseq2 workflow
//     smartseq2_align (add_gfp.out.genome, add_gfp.out.gtf, params.smartseq2_sample_csv)

//     // merge and extract gfp counts
//     process_counts (smartseq2_align.out.htseq_count_files)

//     // run velocyto
//     velocyto_smartseq2 (smartseq2_align.out.bam_files, ch_gtf)

//     // view output
//     process_counts.out.processed_counts | view
//     velocyto_smartseq2.out.velocyto_counts | view
// }



// // /*------------------------------------------------------------------------------------*/
// // /* Workflow to run peaks intersect
// // --------------------------------------------------------------------------------------*/
// Channel
//     .value(file(params.genome, checkIfExists: true))
//     .set {ch_genome}

// Channel
//     .value(file(params.gtf, checkIfExists: true))
//     .set {ch_gtf}

// include {peak_intersect} from "$baseDir/workflows/peak_intersect/main.nf"
// include {fastq_metadata as parse_metadata} from "$baseDir/luslab-nf-modules/tools/metadata/main.nf"
// include {r_analysis as enhancer_profile} from "$baseDir/modules/r_analysis/main.nf"

// workflow {

//     // identify putative enhancers (overlap ATAC + ChIP) and plot peak profiles across enhancers
//     parse_metadata (params.peak_intersect_sample_csv)
//     peak_intersect (parse_metadata.out, ch_genome, ch_gtf)
//     enhancer_profile( params.modules['enhancer_profile'], parse_metadata.out.map{ [it[1]]}.flatten().collect().combine(peak_intersect.out.putative_enhancers))
// }
        


/*------------------------------------------------------------------------------------*/
/* Workflow to run sox8 DEA
--------------------------------------------------------------------------------------*/

// include {r_analysis as lmx1a_dea} from "$baseDir/modules/r_analysis/main.nf"
// include {r_analysis as sox8_dea} from "$baseDir/modules/r_analysis/main.nf"
// include {r_analysis as enhancer_profile} from "$baseDir/modules/r_analysis/main.nf"

// Channel
//     .fromPath(params.lmx1a_counts)
//     .set { ch_lmx1a_counts }

// Channel
//     .fromPath(params.sox8_counts)
//     .set { ch_sox8_counts }


// Channel
//     .fromPath(params.chip_bigwig)
//     .set { ch_chip_bigwig }

// workflow {
//     // lmx1a_dea( params.modules['lmx1a_dea'], ch_lmx1a_counts )
//     // sox8_dea( params.modules['sox8_dea'], ch_sox8_counts )
//     enhancer_profile( params.modules['enhancer_profile'], ch_chip_bigwig.combine(peak_intersect.out.putative_enhancers) )

// }




// // Run merge scRNAseq test
// /*------------------------------------------------------------------------------------*/
// /* Module inclusions
// --------------------------------------------------------------------------------------*/
// include {process_counts} from "$baseDir/workflows/process_counts/main.nf"

// /*------------------------------------------------------------------------------------*/
// /* Set input channel
// --------------------------------------------------------------------------------------*/

// testData = [
//     [[sample_id:"1"], "$baseDir/testData/process_counts/ss8-TSS_P1_A1.txt"],
//     [[sample_id:"2"], "$baseDir/testData/process_counts/ss8-TSS_P1_A2.txt"],
//     [[sample_id:"3"], "$baseDir/testData/process_counts/ss8-TSS_P1_A3.txt"],
//     [[sample_id:"4"], "$baseDir/testData/process_counts/ss8-TSS_P1_A4.txt"],
//     [[sample_id:"5"], "$baseDir/testData/process_counts/ss8-TSS_P1_A5.txt"],
//     [[sample_id:"6"], "$baseDir/testData/process_counts/ss8-TSS_P1_A6.txt"]
// ]
// Channel
//     .from(testData)
//     .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
//     .set { ch_testData }

// /*------------------------------------------------------------------------------------*/
// /* Workflow to run process_counts DEA
// --------------------------------------------------------------------------------------*/

// workflow {
//     process_counts( ch_testData )
// }







// Run merge smartseq bulk test
/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {r_analysis as merge_smartseq_bulk} from "$baseDir/modules/r_analysis/main.nf"
include {process_counts} from "$baseDir/workflows/process_counts/main.nf"

/*------------------------------------------------------------------------------------*/
/* Set input channel
--------------------------------------------------------------------------------------*/

Channel
    .fromPath("$baseDir/output/htseq_count/*.txt")
    .map { [[sample_id:"test"], file(it, checkIfExists: true) ]}
    .set { ch_testData }

Channel
    .fromPath(params.sox8_counts)
    .set { ch_sox8_counts }

/*------------------------------------------------------------------------------------*/
/* Workflow to run process_counts DEA
--------------------------------------------------------------------------------------*/

workflow {
    process_counts( ch_testData )
    // merge smartseq single cell counts with bulk data and output in Antler format
    merge_smartseq_bulk( params.modules['merge_smartseq_bulk'], process_counts.out.processed_counts.listFiles().combine(ch_sox8_counts).flatten().collect() )
}