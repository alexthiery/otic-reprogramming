#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include {fastq_metadata} from "$baseDir/luslab-nf-modules/tools/metadata/main.nf"
include {bedtools_intersect2} from "$baseDir/luslab-nf-modules/tools/bedtools/main.nf"
include {homer_annotatePeaks} from "$baseDir/luslab-nf-modules/tools/homer/main.nf"



/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow peak_intersect {
    take:
        sample_csv

    main:
        fastq_metadata (sample_csv)

        bedtools_intersect2(params.modules['bedtools_intersect2'], fastq_metadata.out.filter{ it[0].sample_id == 'K27Ac' }, fastq_metadata.out.filter{ it[0].sample_id == 'ATAC' }.map{ it[1] } )
        bedtools_intersect2.out | view
        // cutadapt (params.modules['cutadapt'], smartseq2_fastq_metadata.out)

    // emit:
    //     velocyto_counts = velocyto_run_smartseq2.out.velocyto
    //     merged_counts = merge_counts.out.counts
}
