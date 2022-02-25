#!/usr/bin/env nextflow

// Import groovy libs
import groovy.transform.Synchronized

// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/
include {projectHeader} from "$baseDir/../modules/projectHeader/main.nf"
include {extract_gtf_annotations} from "$baseDir/../modules/genome-tools/main.nf"
include {enhancer_analysis} from "$baseDir/../workflows/enhancer_analysis/main.nf"
include {r_analysis as lmx1a_dea; r_analysis as sox8_dea; r_analysis as smartseq_analysis} from "$baseDir/../modules/r_analysis/main.nf"

/*-----------------------------------------------------------------------------------------------------------------------------
Init
-------------------------------------------------------------------------------------------------------------------------------*/
// Show banner
log.info projectHeader()

// Header log info
def summary = [:]
summary['Run Name']               = workflow.runName
summary['Input File']          = params.input
summary['Fasta File']             = params.fasta
summary['GTF File']               = params.gtf
summary['Max Resources']          = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine)     summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output Dir']             = params.outdir
summary['Launch Dir']             = workflow.launchDir
summary['Working Dir']            = workflow.workDir
summary['Script Dir']             = workflow.projectDir
summary['User']                   = workflow.userName
log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join('\n')
log.info "-\033[2m--------------------------------------------------\033[0m-"


/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

Channel
    .value(file(params.fasta, checkIfExists: true))
    .set {ch_fasta}

Channel
    .value(file(params.gtf, checkIfExists:true))
    .set {ch_gtf}

Channel
    .fromPath( params.input )
    .splitCsv(header:true)
    .map { row -> processRow(row) }
    .set { metadata }

metadata
    .filter{ it[0].sample_id == 'sox8_alignment_out' }
    .map { row -> row[1].collect{ file(it+"/star/featurecounts.merged.counts.tsv", checkIfExists: true) } }
    .set { ch_sox8_readcounts }

metadata
    .filter{ it[0].sample_id == 'lmx1a_alignment_out' }
    .map { row -> row[1].collect{ file(it+"/star/featurecounts.merged.counts.tsv", checkIfExists: true) } }
    .set { ch_lmx1a_readcounts }

metadata
    .filter{ it[0].sample_id == 'ChIP_alignment_out' }
    .map { row -> [row[0], row[1].collect{ file(it+"/bwa/mergedLibrary/bigwig", checkIfExists: true) }] }
    .map { row -> listFiles(row, '.*.bigWig') }
    .flatMap { row -> enumerateDir(row) }
    .set { ch_chip_bigwig }

metadata
    .filter{ it[0].sample_id == 'ATAC_alignment_out' }
    .map { row -> [row[0], row[1].collect{ file(it+"/bwa/mergedLibrary/bigwig", checkIfExists: true) }] }
    .map { row -> listFiles(row, '.*.bigWig') }
    .flatMap { row -> enumerateDir(row) }
    .set { ch_atac_bigwig }

metadata
    .filter{ it[0].sample_id == 'ChIP_alignment_out' }
    .map { row -> [row[0], row[1].collect{ file(it+"/bwa/mergedLibrary/macs/broadPeak", checkIfExists: true) }] }
    .map { row -> listFiles(row, '.*.broadPeak') }
    .flatMap { row -> enumerateDir(row) }
    .set { ch_chip_peaks }

metadata
    .filter{ it[0].sample_id == 'ATAC_alignment_out' }
    .map { row -> [row[0], row[1].collect{ file(it+"/bwa/mergedLibrary/macs/narrowPeak", checkIfExists: true) }] }
    .map { row -> listFiles(row, '.*.narrowPeak') }
    .flatMap { row -> enumerateDir(row) }
    .set { ch_atac_peaks }

metadata
    .filter{ it[0].sample_id == 'smartseq2_alignment_out' }
    .map { row -> [row[0], row[1].collect{ file(it+"/merged_counts/output", checkIfExists: true) }] }
    .map { row -> listFiles(row, '.*.csv') }
    .flatMap { row -> row[1] }
    .set { ch_smartseq2_counts }

metadata
    .filter{ it[0].sample_id == 'smartseq2_alignment_out' }
    .map { row -> row[1].collect{ file(it+"/velocyto/onefilepercell_ss8-TSS_P2_C10_75_and_others_TVIUH.loom", checkIfExists: true) } }
    .set { ch_smartseq2_velocyto }


// /*------------------------------------------------------------------------------------*/
// /* Workflow to full downstream analysis
// --------------------------------------------------------------------------------------*/

workflow {
    //  Run differential expression analysis for lmx1a vs sox3U3
    lmx1a_dea( params.modules['lmx1a_dea'], ch_lmx1a_readcounts )

    //  Run differential expression analysis for sox8 over expression vs control
    sox8_dea( params.modules['sox8_dea'], ch_sox8_readcounts )

    // Identify putative enhancers (overlap ATAC + ChIP) and run peak profiles, motif enrichment and functional enrichment analysis
    enhancer_analysis( ch_chip_bigwig, ch_atac_bigwig, ch_chip_peaks, ch_atac_peaks, ch_fasta, ch_gtf.map { [[:], [it]]} )

    // Extract gene annotations from gtf
    extract_gtf_annotations( params.modules['extract_gtf_annotations'], ch_gtf )

    //  Run smartseq2 Antler analysis
    smartseq_analysis( params.modules['smartseq_analysis'], ch_smartseq2_counts
                                                                .combine(ch_smartseq2_velocyto)
                                                                .combine(extract_gtf_annotations.out) )
}

// /*------------------------------------------------------------------------------------*/
// /*  Define custom functions for parsing samplesheet.csv
// --------------------------------------------------------------------------------------*/

def processRow(LinkedHashMap row) {
    def meta = [:]
    meta.sample_id = row.sample_id

    for (Map.Entry<String, ArrayList<String>> entry : row.entrySet()) {
        String key = entry.getKey();
        String value = entry.getValue();
    }

    def array = [ meta, [ row.data ] ]
    return array
}

@Synchronized
def listFiles(row, glob){
    file_array = []
    files = row[1].get(0).listFiles()
    for(def file:files){
        if(file.toString().matches(glob)){
            file_array.add(file)
        }
    }
    array = [row[0], [file_array]]
    return array
}

@Synchronized
def enumerateDir(metadata){
    def array = []
    for (def files : metadata[1].flatten()){
        String s1 = files.getName().split("_")[0]
        temp_meta = metadata[0].getClass().newInstance(metadata[0])
        temp_meta.remove("sample_id")
        temp_meta.put("sample_id", s1)
        array.add([ temp_meta, [file(files, checkIfExists: true)]])
    }
    return array
}