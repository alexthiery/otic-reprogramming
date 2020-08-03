#!/usr/bin/env nextflow

// Define DSL2
nextflow.enable.dsl=2

/*------------------------------------------------------------------------------------*/
/* Define params
--------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------*/
/* Module inclusions
--------------------------------------------------------------------------------------*/

// include {cutadapt} from '../main.nf' 

/*------------------------------------------------------------------------------------*/
/* Define input channels
--------------------------------------------------------------------------------------*/

// command to move subset of files for testing
// for file in $(find /Volumes/lab-luscomben/home/users/thierya/raw_data/ailin_scRNAseq/Samples/*/Files/*1234.fastq.gz | head -2); do rsync -azP $file /Users/alex/dev/repos/otic-reprogramming/data/ss8_9 ; done
// for file in $(find /Volumes/lab-luscomben/home/users/thierya/raw_data/ailin_scRNAseq/ss11_123fq/* | head -4); do rsync -azP $file /Users/alex/dev/repos/otic-reprogramming/data/ss11_123fq ; done
// // for file in $(find /Volumes/lab-luscomben/home/users/thierya/raw_data/ailin_scRNAseq/ss15_123fq/* | head -4); do rsync -azP $file /Users/alex/dev/repos/otic-reprogramming/data/ss15_123fq ; done


include {fastq_metadata} from "$baseDir/../../luslab-nf-modules/tools/metadata/main.nf"


workflow {
    fastq_metadata ("$baseDir/sample.csv")

    fastq_metadata.out.view()
}



// def unlistFastqDir (dir, strip1 = false, strip2 = false) {

//     filesDir = file(dir)
//     filesDir.sort()

//     def list = []
//     if (strip2){
//         println "sequencing mode: paired end"
//         for( def file : filesDir ) {
//             if (file.getName().contains(strip1)){
//                 s1 = file.getName().replaceAll(strip1, "")
//             }
//             else if (file.getName().contains(strip2)){
//                 s1 = file.getName().replaceAll(strip2, "")
//             }
//             if (list.isEmpty()){
//                 list.add([s1, file])
//             }
//             else {
//                 x = 0
//                 for (def it : list) {
//                     if(it.contains(s1)){
//                         it.add(file)
//                         x = 1
//                         }
//                     };
//                 if(x == 1){
//                     continue
//                 }
//                 else{
//                     list.add([s1, file])
//                 }
//             }
//         }
//     }
//     else {
        // println "sequencing mode: single end"
        // for( def file : filesDir ) {
        //     String s1 = file.getName().replaceAll(strip1, "")
        //     list.add([s1, file])
        // }
//     }
//     return list
// }











// Channel

// println samples[1][1]




// //Define test data input channels

// //Single end
// Channel
//     .from(testDataSingleEnd)
//     .map { row -> [ row[0], file(row[1], checkIfExists: true) ] }
//     .set {ch_fastq_single_end}

// //Paired-end
// Channel
//     .from(testDataPairedEnd)
//     .map { row -> [ row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]]}
//     .set {ch_fastq_paired_end}
// /*------------------------------------------------------------------------------------*/
// /* Run tests
// --------------------------------------------------------------------------------------*/
  
// workflow {
//     // Run cutadapt
//     //cutadapt ( ch_fastq_single_end )
//     cutadapt ( ch_fastq_paired_end )

//     // Collect file names and view output
//     cutadapt.out.trimmedReads | view
// }