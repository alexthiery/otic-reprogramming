#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include {awk as awk_enhancer_filter; awk as awk_gtf_filter; cut} from "$baseDir/luslab-nf-modules/tools/luslab_linux_tools/main.nf"
include {bedtools_intersect} from "$baseDir/luslab-nf-modules/tools/bedtools/main.nf"
include {bedtools_subtract} from "$baseDir/luslab-nf-modules/tools/bedtools/main.nf"
include {homer_annotate_peaks; homer_find_motifs} from "$baseDir/luslab-nf-modules/tools/homer/main.nf"
include {r_analysis as enhancer_profile; r_analysis as plot_motifs; r_analysis as functional_enrichment_analysis; r_analysis as peak_annotations_frequency} from "$baseDir/modules/r_analysis/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow enhancer_analysis {
    take:
        peaks
        genome
        gtf

    main:
        // Keep only protein coding genes in gtf
        awk_gtf_filter(params.modules['awk_gtf_filter'], gtf)

        // Intersect ATAC peaks with H3K27Ac peaks
        bedtools_intersect(params.modules['bedtools_intersect'], peaks.filter{ it[0].sample_id == 'ATAC' }, peaks.filter{ it[0].sample_id == 'H3K27Ac' }.map{ it[1] } )

        // Remove any peaks which also have hits for H3K27me3
        bedtools_subtract(params.modules['bedtools_subtract'], bedtools_intersect.out, peaks.filter{ it[0].sample_id == 'H3K27me3' }.map{ it[1] } )

        // Annotate remaining peaks
        homer_annotate_peaks(params.modules['homer_annotate_peaks'], bedtools_subtract.out, genome, awk_gtf_filter.out.file_no_meta)

        // Plot distribution of intersected peak annotations
        peak_annotations_frequency(params.modules['peak_annotations_frequency'], homer_annotate_peaks.out.map{it[1]})

        // Remove peaks in protein coding promoters and  protein coding exons
        awk_enhancer_filter(params.modules['awk_enhancer_filter'], homer_annotate_peaks.out)

        // Run functional enrichment analysis on annotated putative enhancers
        functional_enrichment_analysis(params.modules['functional_enrichment_analysis'], awk_enhancer_filter.out.file_no_meta)

        // Plot ChIP and ATAC profile across enhancers
        enhancer_profile( params.modules['enhancer_profile'], peaks.map{ [it[1]]}.flatten().collect().combine(awk_enhancer_filter.out.file_no_meta))

        // Convert awk output to bed file
        cut(params.modules['cut'], awk_enhancer_filter.out.file)

        // Run motif enrichment analysis on remaining peaks
        homer_find_motifs(params.modules['homer_find_motifs'], awk_enhancer_filter.out.file, genome)

        // Generate motif plot
        plot_motifs( params.modules['plot_motifs'], homer_find_motifs.out.enrichedMotifs.map{it[1]} )
}
