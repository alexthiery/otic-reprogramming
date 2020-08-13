#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include {fastq_metadata} from "$baseDir/luslab-nf-modules/tools/metadata/main.nf"
include {bedtools_intersect} from "$baseDir/luslab-nf-modules/tools/bedtools/main.nf"
include {homer_annotatePeaks} from "$baseDir/luslab-nf-modules/tools/homer/main.nf"



/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow peak_intersect {
    take:
        sample_csv

    main:
        fastq_metadata (sample_csv)
        fastq_metadata.out | view

        bedtools_intersect(fastq_metadata.out.filter { it[0] == 'K27Ac' }, fastq_metadata.out.filter { it[0] == 'ATAC' })
        bedtools_intersect.out | view
        // cutadapt (params.modules['cutadapt'], smartseq2_fastq_metadata.out)

    // emit:
    //     velocyto_counts = velocyto_run_smartseq2.out.velocyto
    //     merged_counts = merge_counts.out.counts
}
