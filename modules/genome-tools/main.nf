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
    
    container 'nfcore/base:1.10.2'

    input:
        val(opts)
        path(fa)
        path(gfp_seq)

    output:
        path("${fa.baseName}_GFP.fa")

    script:
        """
        cat ${gfp_seq} | sed 's/>.*/>GFP/' > GFP.fa
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

    container 'nfcore/base:1.10.2'

    input:
        val(opts)
        path(gtf)
        path(gfp_seq)

    output:
        path("${gtf.baseName}_GFP.gtf")

    script:

        //SHELL
        """
        cat ${gfp_seq} | sed 's/>.*/>GFP/' > GFP.fa
        echo -e 'GFP\tunknown\texon\t1\t'"\$(cat GFP.fa | grep -v "^>" | tr -d "\n" | wc -c)"'\t.\t'"${opts.strand}"'\t.\tgene_id "GFP"; transcript_id "GFP"; gene_name "GFP"; gene_biotype "protein_coding";' > GFP.gtf
        cp ${gtf} ${gtf.baseName}_GFP.gtf
        cat GFP.gtf >> ${gtf.baseName}_GFP.gtf
        """
}



process extract_gtf_annotations {
    publishDir "${params.outdir}/${opts.publish_dir}",
        mode: "copy", 
        overwrite: true,
        saveAs: { filename ->
                      if (opts.publish_results == "none") null
                      else filename }

    container 'nfcore/base:1.10.2'

    input:
        val(opts)
        path(gtf)

    output:
        path("${gtf.baseName}_gene_annotations.csv")

    script:
    """
    #!/opt/conda/bin/python
    
    import re

    # create output file to write to
    outfile = open("${gtf.baseName}_gene_annotations.csv", 'a')

    outfile.write("Gene stable ID,Gene name"+"\\n")

    out_genes = []

    with open("${gtf}", 'rt') as gtf:
        for line in gtf:
            # only search lines with gene_id in order to skip header lines
            if 'gene_id' in line:
                gene_id = re.sub('.*gene_id "', '', line)
                gene_id = re.sub('".*', '', gene_id).rstrip()

                if gene_id not in out_genes:
                    out_genes.append(gene_id)
                    
                    if 'gene_name' in line:
                        gene_name = re.sub('.*gene_name "', '', line)
                        gene_name = re.sub('".*', '', gene_name).rstrip()

                        outfile.write(",".join([gene_id, gene_name])+"\\n")

                    else:
                        outfile.write(",".join([gene_id, gene_id])+"\\n")
    """
}