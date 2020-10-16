#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include {awk} from "$baseDir/luslab-nf-modules/tools/luslab_linux_tools/main.nf"
include {bedtools_intersect} from "$baseDir/luslab-nf-modules/tools/bedtools/main.nf"
include {bedtools_subtract} from "$baseDir/luslab-nf-modules/tools/bedtools/main.nf"
include {homer_annotate_peaks; homer_find_motifs} from "$baseDir/luslab-nf-modules/tools/homer/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow peak_intersect {
    take:
        peaks
        genome
        gtf

    main:


        // Intersect ATAC peaks with H3K27Ac peaks
        bedtools_intersect(params.modules['bedtools_intersect'], peaks.filter{ it[0].sample_id == 'ATAC' }, peaks.filter{ it[0].sample_id == 'H3K27Ac' }.map{ it[1] } )

        // Remove any peaks which also have hits for H3K27me3
        bedtools_subtract(params.modules['bedtools_subtract'], bedtools_intersect.out, peaks.filter{ it[0].sample_id == 'H3K27me3' }.map{ it[1] } )

        // Annotate remaining peaks
        homer_annotate_peaks(params.modules['homer_annotate_peaks'], bedtools_subtract.out, genome, gtf)

        // Remove peaks in promoter regions (±2kb TSS) or exons
        awk(params.modules['awk'], homer_annotate_peaks.out)

        // Run motif enrichment analysis on remaining peaks
        homer_find_motifs(params.modules['homer_find_motifs'], awk.out.file, genome)

    emit:
        putative_enhancers = awk.out.file_no_meta
}
