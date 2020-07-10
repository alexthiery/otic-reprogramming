#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for get_crosslinks...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.fai = "$baseDir/input/GRCh38.primary_assembly.genome_chr6_34000000_35000000.fa.fai"

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include getcrosslinks from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

//Define test data 
testDataBam = [
    ['Sample1', "$baseDir/input/sample1.bam"],
    ['Sample2', "$baseDir/input/sample2.bam"]
]

testDataBamBai = [
    ['Sample1', "$baseDir/input/sample1.bam", "$baseDir/input/sample1.bam.bai"],
    ['Sample2', "$baseDir/input/sample2.bam", "$baseDir/input/sample2.bam.bai"]
]

//Define test data input channels

// Fai input channel
Channel
    .fromPath(params.fai, checkIfExists: true)
    .set {ch_fai}

// Bam input channel
Channel
    .from(testDataBam)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .combine( ch_fai )
    .set {ch_bam}

// Bam/bai input channel
Channel
    .from(testDataBamBai)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists:true)] ] }
    .combine( ch_fai )
    .set {ch_bam_bai}


/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run getcrosslinks
    getcrosslinks ( ch_bam )
    //getcrosslinks ( ch_bam_bai )

    // Collect file names and view output
    getcrosslinks.out.crosslinkBed | view
}