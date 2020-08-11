#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

import {bedtools_intersect} from 

testData_ChIP = [
    ['K4me1', "/Users/alex/dev/repos/otic-reprogramming/testData/ChIP/ss8-K4me1_R1_peaks.broadPeak"],
    ['K4me3', "/Users/alex/dev/repos/otic-reprogramming/testData/ChIP/ss8-K4me3_R1_peaks.broadPeak"],
    ['K27Ac', "/Users/alex/dev/repos/otic-reprogramming/testData/ChIP/ss8-K27Ac_R1_peaks.broadPeak"],
    ['K27me3', "/Users/alex/dev/repos/otic-reprogramming/testData/ChIP/ss8-K27me3_R1_peaks.broadPeak"]
]

testData_ATAC = [
    ['Sample1', "/Users/alex/dev/repos/otic-reprogramming/testData/ATAC/ss8_R1.mLb.clN_peaks.narrowPeak"]
]

// Define test data input channel
Channel
    .from(testData_ChIP)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_testData_ChIP}

Channel
    .from(testData_ATAC)
    .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
    .set {ch_testData_ATAC}

workflow {

}

(ch_testData_ChIP.filter { it[0] == 'K27Ac' }.subscribe { println it }, ch_testData_ATAC.filter { it[0] == 'Sample1' }.subscribe { println it })