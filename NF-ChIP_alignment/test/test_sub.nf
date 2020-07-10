#!/usr/bin/env nextflow

// Define DSL2
nextflow.preview.dsl=2

// Log
log.info ("Starting tests for seqtk...")

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

params.seqtk_subsample_args = ''
params.seqtk_subsample_seed = 999
params.seqtk_subsample_number = 100
params.verbose = true

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

include seqtk_subsample from '/Users/alex/dev/repos/otic-reprogramming/luslab-nf-modules/tools/seqtk/main.nf'

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// Define test data
testDataPairedEnd = [
    ['ss11-input', "/Volumes/Extreme SSD/ChIP/SS_11-12/input/WTCHG_45718_279_1_sequence.txt.gz", "/Volumes/Extreme SSD/ChIP/SS_11-12/input/WTCHG_45718_279_2_sequence.txt.gz"],
    ['ss11-K27Ac', "/Volumes/Extreme SSD/ChIP/SS_11-12/K27Ac/WTCHG_45718_214_1_sequence.txt.gz", "/Volumes/Extreme SSD/ChIP/SS_11-12/K27Ac/WTCHG_45718_214_2_sequence.txt.gz"],
    ['ss11-K4me3', "/Volumes/Extreme SSD/ChIP/SS_11-12/K4me3/WTCHG_45718_276_1_sequence.txt.gz", "/Volumes/Extreme SSD/ChIP/SS_11-12/K4me3/WTCHG_45718_276_2_sequence.txt.gz"]
]

//Define test data input channel

Channel
    .from(testDataPairedEnd)
    .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
    .set {ch_fastq_paired_end}

/*------------------------------------------------------------------------------------*/
/* Run tests
--------------------------------------------------------------------------------------*/

workflow {
    // Run seqtk paired end
    seqtk_subsample ( ch_fastq_paired_end )
    seqtk_subsample.out.sampledReads | view
}

//nextflow run ./NF-ChIP_alignment/test/test_sub.nf --outdir /Users/alex/dev/repos/otic-reprogramming/NF-ChIP_alignment/test/test-data