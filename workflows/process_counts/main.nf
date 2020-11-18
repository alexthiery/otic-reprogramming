#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {r_analysis as merge_counts} from "$baseDir/../modules/r_analysis/main.nf"


/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow process_counts {
    take:
        count_files

    main:
        // group counts into a single channel for merge cell counts
        ch_all_counts = count_files
            .map { [file(it[1], checkIfExists: true)] }
            .collect()

        // merge cell counts into csv
        merge_counts (params.modules['merge_counts'], ch_all_counts)

    emit:
        processed_counts = merge_counts.out
}