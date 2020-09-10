#!/usr/bin/env nextflow

process merge_counts {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'alexthiery/modules-rscript:latest'

    input:
        val opts
        tuple val(meta), path(input)

    output:
        tuple val(meta), path("*.csv"), emit: counts

    script:
        //SHELL

        merge_counts_command = "Rscript ${baseDir}/bin/mergeCounts.R"

        if (params.verbose){
            println ("[MODULE] merge counts command: " + merge_counts_command)
        }

        //SHELL
        """
        ${merge_counts_command}
        """
}

