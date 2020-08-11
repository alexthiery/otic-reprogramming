#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/
params.modules['cutadapt'].args = '-a CTGTCTCTTATA'
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {smartseq2_fastq_metadata} from "$baseDir/../../luslab-nf-modules/tools/metadata/main.nf"
include {cutadapt} from "$baseDir/../../luslab-nf-modules/tools/cutadapt/main.nf"
include {hisat2_build; hisat2_splice_sites; hisat2_splice_align} from "$baseDir/../../luslab-nf-modules/tools/hisat2/main.nf"
include {samtools_view as samtools_view_a;samtools_view as samtools_view_b; samtools_sort} from "$baseDir/../../luslab-nf-modules/tools/samtools/main.nf"
include {htseq_count} from "$baseDir/../../luslab-nf-modules/tools/htseq/main.nf"
include {velocyto_run_smartseq2} from "$baseDir/../../luslab-nf-modules/tools/velocyto/main.nf"

// include {assert_channel_count} from "$baseDir/../../luslab-nf-modules/workflows/test_flows/main.nf"


/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

Channel
    .value(file(params.genome))
    .set {ch_genome}

Channel
    .value(file(params.gtf))
    .set {ch_gtf}

workflow {
    smartseq2_fastq_metadata (params.sample_csv)
    cutadapt (params.modules['cutadapt'], smartseq2_fastq_metadata.out)
    hisat2_build ( params.modules['hisat2_build'], ch_genome )
    hisat2_splice_sites ( params.modules['hisat2_splice_sites'], ch_gtf )
    hisat2_splice_align ( params.modules['hisat2_splice_align'], cutadapt.out.fastq, hisat2_build.out.genome_index.collect(), hisat2_splice_sites.out.splice_sites.collect() )

    samtools_view_a ( params.modules['samtools_view_a'], hisat2_splice_align.out.sam )
    samtools_sort ( params.modules['samtools_sort'], samtools_view_a.out.bam )
    samtools_view_b ( params.modules['samtools_view_b'], samtools_sort.out.bam )
    
    velocyto_run_smartseq2 ( params.modules['velocyto_run_smartseq2'], samtools_sort.out.bam, ch_gtf )

    htseq_count ( params.modules['htseq_count'], samtools_view_b.out.bam, ch_gtf )

    // // Collect file names and view output
    htseq_count.out.counts | view

    // //Check count
    // assert_channel_count ( htseq_count.out.counts, "sam", 2)
}




// command to move subset of files for testing
// for file in $(find /Volumes/lab-luscomben/home/users/thierya/raw_data/ailin_scRNAseq/Samples/*/Files/*1234.fastq.gz | head -2); do rsync -azP $file /Users/alex/dev/repos/otic-reprogramming/data/ss8_9 ; done
// for file in $(find /Volumes/lab-luscomben/home/users/thierya/raw_data/ailin_scRNAseq/ss11_123fq/* | head -4); do rsync -azP $file /Users/alex/dev/repos/otic-reprogramming/data/ss11_123fq ; done
// // for file in $(find /Volumes/lab-luscomben/home/users/thierya/raw_data/ailin_scRNAseq/ss15_123fq/* | head -4); do rsync -azP $file /Users/alex/dev/repos/otic-reprogramming/data/ss15_123fq ; done
