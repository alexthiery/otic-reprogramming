#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for bedtools...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.regions_file = [["$baseDir/input/regions_GENCODE_v30.gtf.gz"]]
params.bedtools_intersect_args = '-wa -wb -s'
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions 
--------------------------------------------------------------------------------------*/

include bedtools_intersect from '../main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
/*------------------------------------------------------------------------------------*/

// Define test data
testData = [
    ['Sample1', "$baseDir/input/sample1.xl.bed.gz"],
    ['Sample2', "$baseDir/input/sample2.xl.bed.gz"],
    ['Sample3', "$baseDir/input/sample3.xl.bed.gz"],
    ['Sample4', "$baseDir/input/sample4.xl.bed.gz"],
    ['Sample5', "$baseDir/input/sample5.xl.bed.gz"],
    ['Sample6', "$baseDir/input/sample6.xl.bed.gz"]
]

// Define regions file input channel
Channel
    .from(params.regions_file)
    .map { row -> file(row[0], checkIfExists: true)}
    .set {ch_test_regions_file}

// Define test data input channel
Channel
    .from(testData)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .combine(ch_test_regions_file)
    .set {ch_test_crosslinks}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run bedtools_intersect
    bedtools_intersect ( ch_test_crosslinks )

    // Collect file names and view output
    bedtools_intersect.out.annotatedBed | view
}