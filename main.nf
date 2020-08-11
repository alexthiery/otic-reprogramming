#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.bedtools_intersect_args = '-wa -wb'
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include {bedtools_intersect} from 'luslab-nf-modules/bedtools/homer/main.nf'
include {homer_annotatePeaks} from 'luslab-nf-modules/tools/homer/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
chipData = [
    ['K4me1', "$baseDir/testData/ss8-K4me1_R1_peaks.broadPeak"],
    ['K4me3', "$baseDir/testData/ss8-K4me3_R1_peaks.broadPeak"],
    ['K27Ac', "$baseDir/testData/ss8-K27Ac_R1_peaks.broadPeak"],
    ['K27me3', "$baseDir/testData/ss8-K27me3_R1_peaks.broadPeak"]
]

atacData = [
    ['ATAC', "$baseDir/testData/ss8_R1.mLb.clN_peaks.broadPeak"]
]

// Define test data input channel
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_testData}

ch_testData
    .filter { it[0] == 'K27Ac' }
    .set {ch_K27Ac}

ch_testData
    .filter { it[0] == 'K27me3' }
    .set {ch_K27me3}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

// Run workflow
workflow {
    // intersect peak files
    bedtools_intersect(ch_K27Ac, ch_atacData)

    bedtools_subtract(bedtools_intersect.out, ch_K27me3)
    
    // View outputs
    bedtools_intersect.out | view
}