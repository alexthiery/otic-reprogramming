#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process r_analysis {
    publishDir "${params.outdir}/${opts.publish_dir}",
    mode: "copy",
    overwrite: true,
    saveAs: { filename ->
                    if (opts.publish_results == "none") null
                    else filename }

    container "alexthiery/modules-rscript:latest"

    input:
        val opts
        path input

    output:
        path("output")

    script:

    args = ""
        if(opts.args && opts.args != '') {
            ext_args = opts.args
            args += ext_args.trim()
        }

    """
    Rscript ${opts.script} --cores ${task.cpus} --runtype nextflow ${args}
    """
}