#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Workflow to full downstream analysis
--------------------------------------------------------------------------------------*/

Channel
    .value(file(params.genome, checkIfExists: true))
    .set {ch_genome}

Channel
    .fromPath(params.gtf)
    .map { [[:], [file(it, checkIfExists: true)]]}
    .set {ch_gtf}

Channel
    .fromPath(params.lmx1a_counts)
    .set { ch_lmx1a_counts }

Channel
    .fromPath(params.sox8_counts)
    .set { ch_sox8_counts }


// need to add all inputs for antler together here
// Channel
//     .fromPath(params.smartseq2_counts)
//     .set { ch_smartseq2_counts }


include {enhancer_analysis} from "$baseDir/../workflows/enhancer_analysis/main.nf"
include {fastq_metadata as parse_metadata} from "$baseDir/../luslab-nf-modules/tools/metadata/main.nf"
include {r_analysis as lmx1a_dea; r_analysis as sox8_dea; r_analysis as smartseq_analysis} from "$baseDir/../modules/r_analysis/main.nf"

workflow {
    //  Run differential expression analysis for lmx1a vs sox3U3
    lmx1a_dea( params.modules['lmx1a_dea'], ch_lmx1a_counts )

    //  Run differential expression analysis for sox8 over expression vs control
    sox8_dea( params.modules['sox8_dea'], ch_sox8_counts )

    // Identify putative enhancers (overlap ATAC + ChIP) and run peak profiles, motif enrichment and functional enrichment analysis
    parse_metadata( params.enhancer_analysis_sample_csv )
    enhancer_analysis( parse_metadata.out, ch_genome, ch_gtf )

    //  Run smartseq2 Antler analysis
    // smartseq_analysis( params.modules['smartseq_analysis'], ch_smartseq2_counts )
}