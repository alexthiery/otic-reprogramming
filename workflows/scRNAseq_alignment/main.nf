#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {smartseq2_fastq_metadata} from "$baseDir/luslab-nf-modules/tools/metadata/main.nf"
include {cutadapt} from "$baseDir/luslab-nf-modules/tools/cutadapt/main.nf"
include {hisat2_build; hisat2_splice_sites; hisat2_splice_align} from "$baseDir/luslab-nf-modules/tools/hisat2/main.nf"
include {samtools_view as samtools_view_a;samtools_view as samtools_view_b; samtools_sort} from "$baseDir/luslab-nf-modules/tools/samtools/main.nf"
include {htseq_count} from "$baseDir/luslab-nf-modules/tools/htseq/main.nf"
include {velocyto_run_smartseq2} from "$baseDir/luslab-nf-modules/tools/velocyto/main.nf"
include {merge_counts} from "$baseDir/custom-nf-modules/runR/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow smartseq2_align {
    take:
        genome
        gtf
        sample_csv

    main:
        smartseq2_fastq_metadata (sample_csv)
        
        cutadapt (params.modules['cutadapt'], smartseq2_fastq_metadata.out)
        hisat2_build ( params.modules['hisat2_build'], genome )
        hisat2_splice_sites ( params.modules['hisat2_splice_sites'], gtf )
        hisat2_splice_align ( params.modules['hisat2_splice_align'], cutadapt.out.fastq, hisat2_build.out.genome_index.collect(), hisat2_splice_sites.out.splice_sites.collect() )

        samtools_view_a ( params.modules['samtools_view_a'], hisat2_splice_align.out.sam )
        samtools_sort ( params.modules['samtools_sort'], samtools_view_a.out.bam )
        samtools_view_b ( params.modules['samtools_view_b'], samtools_sort.out.bam )
        
        // group bams into a single channel for velocyto
        ch_velocyto_bam = samtools_sort.out.bam
            .map { [[sample_id:"all_cells"], file(it[1], checkIfExists: true)] }
            .groupTuple(by: 0)

        velocyto_run_smartseq2 ( params.modules['velocyto_run_smartseq2'], ch_velocyto_bam, gtf )

        htseq_count ( params.modules['htseq_count'], samtools_view_b.out.bam, gtf )

        // merge cell counts into csv
        merge_counts (params.modules['merge_counts'], htseq_count.out.counts.collect())

    emit:
        velocyto_counts = velocyto_run_smartseq2.out.velocyto
        merged_counts = merge_counts.out.counts
}

