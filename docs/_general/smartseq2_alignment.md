---
layout: page
label: SmartSeq2 alignment
order: 2
---

## SmartSeq2 alignment

</br>

**If you would like to re-run the alignment, follow the instructions [here]({{site.baseurl}}/general/quick_start#smartseq).**

Our SmartSeq2 single cell RNAseq data is aligned using a custom Nextflow pipeline. This pipeline is built using [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html), which allows for modularisation and re-usability of each tool.

Each of the processes used in the pipeline can be found by following the going to the path in the `include {process} from "<PATH>"` statements in each of the workflows below.

</br>

---

### Workflows

Our custom pipeline integrates both standard alignment and RNA velocity ([Velocyto](http://velocyto.org/)).

For simplicity, the pipeline is split into four sub-workflows:

- [Modify genome to add GFP sequence](#GFP)
- [Align reads](#align)
- [Process readcounts for downstream](#process_reads)
- [Run RNA velocity](#velocyto)

</br>

The main workflow integrates each of these sub-workflows.

<details open><summary class="simple">Expand main workflow</summary>
<p>

```groovy
// Define DSL2
nextflow.enable.dsl=2

include {add_gfp} from "$baseDir/../workflows/add_gfp/main.nf"
include {smartseq2_align} from "$baseDir/../workflows/scRNAseq_alignment/main.nf"
include {process_counts} from "$baseDir/../workflows/process_counts/main.nf"
include {velocyto_smartseq2} from "$baseDir/../workflows/velocyto_smartseq2/main.nf"

Channel
    .value(file(params.fasta, checkIfExists: true))
    .set {ch_fasta}

Channel
    .value(file(params.gtf, checkIfExists: true))
    .set {ch_gtf}

Channel
    .value(file(params.gfp_seq, checkIfExists: true))
    .set {ch_gfp_seq}

workflow {
    // add gfp to genome and gtf
    add_gfp (ch_fasta, ch_gtf, ch_gfp_seq)

    // align using smartseq2 workflow
    smartseq2_align (add_gfp.out.genome, add_gfp.out.gtf, params.input)

    // merge and extract gfp counts
    process_counts (smartseq2_align.out.htseq_count_files)

    // run velocyto
    velocyto_smartseq2 (smartseq2_align.out.bam_files, ch_gtf)

    // view output
    process_counts.out.processed_counts | view
    velocyto_smartseq2.out.velocyto_counts | view
}
```

</details>

---

</br>

**_Modify genome to add GFP sequence_**<a name="GFP"></a>

This sub-workdflow takes a GFP sequence from the base repository and appends the genome provided in order to determine the number of GFP reads downstream.

<details open><summary class="simple">Expand sub-workflow</summary>
<p>

```groovy
// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {add_genome_gfp; add_gtf_gfp} from "$baseDir/../modules/genome-tools/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow add_gfp {
    take:
        genome
        gtf
        gfp_seq

    main:
        add_genome_gfp (params.modules['add_genome_gfp'], genome, gfp_seq)
        add_gtf_gfp (params.modules['add_gtf_gfp'], gtf, gfp_seq)

    emit:
        genome = add_genome_gfp.out
        gtf = add_gtf_gfp.out
}
```

</details>

---

</br>

**_Align reads_**<a name="align"></a>

<details open><summary class="simple">Expand sub-workflow</summary>
<p>

```groovy
// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {smartseq2_fastq_metadata} from "$baseDir/../luslab-nf-modules/tools/metadata/main.nf"
include {cutadapt} from "$baseDir/../luslab-nf-modules/tools/cutadapt/main.nf"
include {hisat2_build; hisat2_splice_sites; hisat2_splice_align} from "$baseDir/../luslab-nf-modules/tools/hisat2/main.nf"
include {samtools_view as samtools_view_a;samtools_view as samtools_view_b; samtools_sort} from "$baseDir/../luslab-nf-modules/tools/samtools/main.nf"
include {htseq_count} from "$baseDir/../luslab-nf-modules/tools/htseq/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow smartseq2_align {
    take:
        genome
        gtf
        sample_csv

    main:
        smartseq2_fastq_metadata (sample_csv)

        cutadapt (params.modules['cutadapt'], smartseq2_fastq_metadata.out)
        hisat2_build ( params.modules['hisat2_build'], genome )
        hisat2_splice_sites ( params.modules['hisat2_splice_sites'], gtf )
        hisat2_splice_align ( params.modules['hisat2_splice_align'], cutadapt.out.fastq, hisat2_build.out.genome_index.collect(), hisat2_splice_sites.out.splice_sites.collect() )

        samtools_view_a ( params.modules['samtools_view_a'], hisat2_splice_align.out.sam )
        samtools_sort ( params.modules['samtools_sort'], samtools_view_a.out.bam )
        samtools_view_b ( params.modules['samtools_view_b'], samtools_sort.out.bam )

        htseq_count ( params.modules['htseq_count'], samtools_view_b.out.bam, gtf )

    emit:
        bam_files = samtools_sort.out.bam
        htseq_count_files = htseq_count.out.counts
}
```

</details>

---

</br>

**_Process read counts for input into Antler_**<a name="process_reads"></a>

<details open><summary class="simple">Expand sub-workflow</summary>
<p>

```groovy
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
```

</details>

---

</br>

**_Run RNA velocity_**<a name="velocyto"></a>

<details open><summary class="simple">Expand sub-workflow</summary>
<p>

```groovy
// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {velocyto_run_smartseq2} from "$baseDir/../luslab-nf-modules/tools/velocyto/main.nf"

/*------------------------------------------------------------------------------------*/
/* Define sub workflow
--------------------------------------------------------------------------------------*/

workflow velocyto_smartseq2 {
    take:
        bam_files
        gtf

    main:
        // group bams into a single channel for velocyto
        ch_velocyto_bam = bam_files
            .map { [[sample_id:"all_cells"], file(it[1], checkIfExists: true)] }
            .groupTuple(by: 0)

        velocyto_run_smartseq2 ( params.modules['velocyto_run_smartseq2'], ch_velocyto_bam, gtf )

    emit:
        velocyto_counts = velocyto_run_smartseq2.out.velocyto
}
```

</details>

</br>
