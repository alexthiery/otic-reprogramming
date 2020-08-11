#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

include {smartseq2_align} from "$baseDir/workflows/scRNAseq_alignment/main.nf"

Channel
    .value(file(params.genome))
    .set {ch_genome}

Channel
    .value(file(params.gtf))
    .set {ch_gtf}

workflow {
    smartseq2_align (ch_genome, ch_gtf, params.sample_csv)
}




// command to move subset of files for testing
// for file in $(find /Volumes/lab-luscomben/home/users/thierya/raw_data/ailin_scRNAseq/Samples/*/Files/*1234.fastq.gz | head -2); do rsync -azP $file /Users/alex/dev/repos/otic-reprogramming/data/ss8_9 ; done
// for file in $(find /Volumes/lab-luscomben/home/users/thierya/raw_data/ailin_scRNAseq/ss11_123fq/* | head -4); do rsync -azP $file /Users/alex/dev/repos/otic-reprogramming/data/ss11_123fq ; done
// // for file in $(find /Volumes/lab-luscomben/home/users/thierya/raw_data/ailin_scRNAseq/ss15_123fq/* | head -4); do rsync -azP $file /Users/alex/dev/repos/otic-reprogramming/data/ss15_123fq ; done

// mkdir -p /camp/home/thierya/working/raw_data/test/ss8_9
// mkdir -p /camp/home/thierya/working/raw_data/test/ss11_123fq
// mkdir -p /camp/home/thierya/working/raw_data/test/ss15_123fq
// for file in $(find /camp/home/thierya/working/raw_data/ailin_scRNAseq/Samples/*/Files/*1234.fastq.gz | head -2); do rsync -azP $file /camp/home/thierya/working/raw_data/test/ss8_9 ; done
// for file in $(find /camp/home/thierya/working/raw_data/ailin_scRNAseq/ss11_123fq/* | head -2); do rsync -azP $file /camp/home/thierya/working/raw_data/test/ss11_123fq ; done
// for file in $(find /camp/home/thierya/working/raw_data/ailin_scRNAseq/ss15_123fq/* | head -2); do rsync -azP $file /camp/home/thierya/working/raw_data/test/ss15_123fq ; done

// workflow {
//     smartseq2_align (ch_genome, ch_gtf, params.sample_csv, params.modules['cutadapt'], params.modules['hisat2_build'], params.modules['hisat2_splice_sites'], params.modules['hisat2_splice_align'],
//     params.modules['samtools_view_a'], params.modules['samtools_sort'], params.modules['samtools_view_b'], params.modules['velocyto_run_smartseq2'], params.modules['htseq_count'])
// }

// /*------------------------------------------------------------------------------------*/
// /* Define params
// --------------------------------------------------------------------------------------*/

// params.bedtools_intersect_args = '-wa -wb'
// params.verbose = true

// /*------------------------------------------------------------------------------------*/
// /* Module inclusions 
// --------------------------------------------------------------------------------------*/

// include {bedtools_intersect} from 'luslab-nf-modules/bedtools/homer/main.nf'
// include {homer_annotatePeaks} from 'luslab-nf-modules/tools/homer/main.nf'

// /*------------------------------------------------------------------------------------*/
// /* Define input channels
// /*------------------------------------------------------------------------------------*/

// // Define test data
// chipData = [
//     ['K4me1', "$baseDir/testData/ss8-K4me1_R1_peaks.broadPeak"],
//     ['K4me3', "$baseDir/testData/ss8-K4me3_R1_peaks.broadPeak"],
//     ['K27Ac', "$baseDir/testData/ss8-K27Ac_R1_peaks.broadPeak"],
//     ['K27me3', "$baseDir/testData/ss8-K27me3_R1_peaks.broadPeak"]
// ]

// atacData = [
//     ['ATAC', "$baseDir/testData/ss8_R1.mLb.clN_peaks.broadPeak"]
// ]

// // Define test data input channel
// Channel
//     .from(testData)
//     .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
//     .set {ch_testData}

// ch_testData
//     .filter { it[0] == 'K27Ac' }
//     .set {ch_K27Ac}

// ch_testData
//     .filter { it[0] == 'K27me3' }
//     .set {ch_K27me3}

// /*------------------------------------------------------------------------------------*/
// /* Run tests
// --------------------------------------------------------------------------------------*/

// // Run workflow
// workflow {
//     // intersect peak files
//     bedtools_intersect(ch_K27Ac, ch_atacData)

//     bedtools_subtract(bedtools_intersect.out, ch_K27me3)
    
//     // View outputs
//     bedtools_intersect.out | view
// }