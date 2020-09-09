#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Generic awk process

process awk {
    input:
      tuple val(meta), path(input_file)
      val args

    output:
        tuple val(meta), path("*.*"), emit: file
        path "*.*", emit: fileNoMeta

    script: 
    """
    FILE=$input_file
    awk ${args} $input_file > "\${FILE%.*}"
    """
}