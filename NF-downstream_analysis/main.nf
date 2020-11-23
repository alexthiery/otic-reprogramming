#!/usr/bin/env nextflow

// Import groovy libs
import groovy.transform.Synchronized

// Define DSL2
nextflow.enable.dsl=2

// /*------------------------------------------------------------------------------------*/
// /* Workflow to full downstream analysis
// --------------------------------------------------------------------------------------*/

Channel
    .value(file(params.genome, checkIfExists: true))
    .set {ch_genome}

Channel
    .fromPath(params.gtf)
    .map { [[:], [file(it, checkIfExists: true)]]}
    .set {ch_gtf}

// Channel
//     .fromPath(params.lmx1a_counts)
//     .set { ch_lmx1a_counts }

// Channel
//     .fromPath(params.sox8_counts)
//     .set { ch_sox8_counts }


// // need to add all inputs for antler together here
// // Channel
// //     .fromPath(params.smartseq2_counts)
// //     .set { ch_smartseq2_counts }


// include {enhancer_analysis} from "$baseDir/../workflows/enhancer_analysis/main.nf"
// include {fastq_metadata as parse_metadata} from "$baseDir/../luslab-nf-modules/tools/metadata/main.nf"
// include {r_analysis as lmx1a_dea; r_analysis as sox8_dea; r_analysis as smartseq_analysis} from "$baseDir/../modules/r_analysis/main.nf"

// workflow {
//     //  Run differential expression analysis for lmx1a vs sox3U3
//     lmx1a_dea( params.modules['lmx1a_dea'], ch_lmx1a_counts )

//     //  Run differential expression analysis for sox8 over expression vs control
//     sox8_dea( params.modules['sox8_dea'], ch_sox8_counts )

//     // Identify putative enhancers (overlap ATAC + ChIP) and run peak profiles, motif enrichment and functional enrichment analysis
//     parse_metadata( params.enhancer_analysis_sample_csv )
//     enhancer_analysis( parse_metadata.out, ch_genome, ch_gtf )

//     //  Run smartseq2 Antler analysis
//     // smartseq_analysis( params.modules['smartseq_analysis'], ch_smartseq2_counts )
// }




Channel
    .fromPath( params.input )
    .splitCsv(header:true)
    .map { row -> processRow(row) }
    .set { metadata }

// metadata
//     .filter{ it[0].sample_id == 'sox8_alignment_out' }
//     .map { row -> row[1].collect{ file(it+"/star/featurecounts/featurecounts.merged.counts.tsv", checkIfExists: true) } }
//     .set { ch_sox8_readcounts }

// metadata
//     .filter{ it[0].sample_id == 'lmx1a_alignment_out' }
//     .map { row -> row[1].collect{ file(it+"/star/featurecounts/featurecounts.merged.counts.tsv", checkIfExists: true) } }
//     .set { ch_lmx1a_readcounts }

metadata
    .filter{ it[0].sample_id == 'ChIP_alignment_out' }
    .map { row -> [row[0], row[1].collect{ file(it+"/bwa/mergedLibrary/bigwig") }] }
    .map { row -> listFiles(row, '.*.bigWig') }
    .flatMap { row -> enumerateDir(row) }
    .set { chip_bigwig }

metadata
    .filter{ it[0].sample_id == 'ATAC_alignment_out' }
    .map { row -> [row[0], row[1].collect{ file(it+"/bwa/mergedLibrary/bigwig") }] }
    .map { row -> listFiles(row, '.*.bigWig') }
    .flatMap { row -> enumerateDir(row) }
    .set { atac_bigwig }

metadata
    .filter{ it[0].sample_id == 'ChIP_alignment_out' }
    .map { row -> [row[0], row[1].collect{ file(it+"/bwa/mergedLibrary/macs/broadPeak") }] }
    .map { row -> listFiles(row, '.*.broadPeak') }
    .flatMap { row -> enumerateDir(row) }
    .set { chip_peaks }

metadata
    .filter{ it[0].sample_id == 'ATAC_alignment_out' }
    .map { row -> [row[0], row[1].collect{ file(it+"/bwa/mergedLibrary/macs/narrowPeak") }] }
    .map { row -> listFiles(row, '.*.narrowPeak') }
    .flatMap { row -> enumerateDir(row) }
    .set { atac_peaks }

include {enhancer_analysis} from "$baseDir/../workflows/enhancer_analysis/main.nf"
// // // include {fastq_metadata as parse_metadata} from "$baseDir/../luslab-nf-modules/tools/metadata/main.nf"
// // include {r_analysis as lmx1a_dea; r_analysis as sox8_dea; r_analysis as smartseq_analysis} from "$baseDir/../modules/r_analysis/main.nf"

workflow {
    // //  Run differential expression analysis for lmx1a vs sox3U3
    // lmx1a_dea( params.modules['lmx1a_dea'], ch_lmx1a_readcounts )

    // //  Run differential expression analysis for sox8 over expression vs control
    // sox8_dea( params.modules['sox8_dea'], ch_sox8_readcounts )

    // Identify putative enhancers (overlap ATAC + ChIP) and run peak profiles, motif enrichment and functional enrichment analysis
    enhancer_analysis( chip_bigwig, atac_bigwig, chip_peaks, atac_peaks, ch_genome, ch_gtf )

    //  Run smartseq2 Antler analysis
    // smartseq_analysis( params.modules['smartseq_analysis'], ch_smartseq2_counts )
}







//  Custom functions for parsing samplesheet.csv
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
