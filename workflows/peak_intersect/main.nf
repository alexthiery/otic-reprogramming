#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include {fastq_metadata} from "$baseDir/luslab-nf-modules/tools/metadata/main.nf"
include {bedtools_intersect} from "$baseDir/luslab-nf-modules/tools/bedtools/main.nf"
include {bedtools_subtract} from "$baseDir/luslab-nf-modules/tools/bedtools/main.nf"
include {homer_annotate_peaks; homer_find_motifs} from "$baseDir/luslab-nf-modules/tools/homer/main.nf"



/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow peak_intersect {
    take:
        sample_csv
        genome
        gtf

    main:
        fastq_metadata (sample_csv)

        bedtools_intersect(params.modules['bedtools_intersect'], fastq_metadata.out.filter{ it[0].sample_id == 'K27Ac' }, fastq_metadata.out.filter{ it[0].sample_id == 'ATAC' }.map{ it[1] } )

        bedtools_subtract(params.modules['bedtools_subtract'], bedtools_intersect.out, fastq_metadata.out.filter{ it[0].sample_id == 'K4me3' }.map{ it[1] } )

        homer_annotate_peaks(params.modules['homer_annotate_peaks'], bedtools_subtract.out, genome, gtf)

        homer_annotate_peaks.out | view

        homer_find_motifs(params.modules['homer_find_motifs'], bedtools_subtract.out, genome)

        homer_find_motifs.out | view
        // cutadapt (params.modules['cutadapt'], smartseq2_fastq_metadata.out)

    // emit:
    //     velocyto_counts = velocyto_run_smartseq2.out.velocyto
    //     merged_counts = merge_counts.out.counts
}
