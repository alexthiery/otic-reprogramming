#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include {add_gfp} from "$baseDir/../workflows/add_gfp/main.nf"
include {smartseq2_align} from "$baseDir/../workflows/scRNAseq_alignment/main.nf"
include {process_counts} from "$baseDir/../workflows/process_counts/main.nf"
include {velocyto_smartseq2} from "$baseDir/../workflows/velocyto_smartseq2/main.nf"

Channel
    .value(file(params.fasta, checkIfExists: true))
    .set {ch_fasta}

Channel
    .value(file(params.gtf, checkIfExists: true))
    .set {ch_gtf}

Channel
    .value(file(params.gfp_seq, checkIfExists: true))
    .set {ch_gfp_seq}

workflow {
    // add gfp to genome and gtf
    add_gfp (ch_fasta, ch_gtf, ch_gfp_seq)

    // align using smartseq2 workflow
    smartseq2_align (add_gfp.out.genome, add_gfp.out.gtf, params.input)

    // merge and extract gfp counts
    process_counts (smartseq2_align.out.htseq_count_files)

    // run velocyto
    velocyto_smartseq2 (smartseq2_align.out.bam_files, ch_gtf)

    // view output
    process_counts.out.processed_counts | view
    velocyto_smartseq2.out.velocyto_counts | view
}