#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2


process add_genome_gfp {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    input:
        val(opts)
        path(fa)

    output:
        path("${fa.baseName}_GFP.fa")

    script:

        //SHELL
        """
        wget -O GFP_orig.fa ${opts.gfp_url}
        cat GFP_orig.fa | sed 's/>.*/>GFP/' > GFP.fa
        cp ${fa} ${fa.baseName}_GFP.fa
        cat GFP.fa >> ${fa.baseName}_GFP.fa
        """
}




process add_gtf_gfp {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    input:
        val(opts)
        path(gtf)

    output:
        path("${gtf.baseName}_GFP.gtf")

    script:

        //SHELL
        """
        wget -O GFP_orig.fa ${opts.gfp_url}
        cat GFP_orig.fa | sed 's/>.*/>GFP/' > GFP.fa
        echo -e "GFP\tunknown\texon\t1\t\$(cat GFP.fa | grep -v "^>" | tr -d "\n" | wc -c)\t.\t+\t.\tgene_id 'GFP'; transcript_id 'GFP'; gene_name 'GFP'; gene_biotype 'protein_coding';" > GFP.gtf
        cp ${gtf} ${gtf.baseName}_GFP.gtf
        cat GFP.gtf >> ${gtf.baseName}_GFP.gtf
        """
}