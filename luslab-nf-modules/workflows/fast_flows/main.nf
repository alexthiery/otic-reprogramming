#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

include {region2bed} from '../../tools/luslab_genome_tools/main.nf'
include {seqtk_subseq} from '../../tools/seqtk/main.nf'
include {samtools_faidx} from '../../tools/samtools/main.nf'

// Define workflow to subset and index a genome region fasta file
workflow subset_genome {
    take: fasta
    take: region
    main:
        // Create bed from region
        region2bed( region )

        // Subset fasta file
        seqtk_subseq( params.modules['seqtk_subseq'], fasta, region2bed.out.bed )

        // Index the fasta file
        samtools_faidx( params.modules['samtools_faidx'], seqtk_subseq.out.subset )

    emit:
        fasta = samtools_faidx.out.fasta
}