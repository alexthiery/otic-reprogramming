#!/usr/bin/env nextflow

process merge_counts {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'alexthiery/10x_neural_tube:v1.0'

    input:
        val opts
        path(input)

    output:
        path("${opts.output}")

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
