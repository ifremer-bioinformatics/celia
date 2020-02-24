#!/usr/bin/env nextflow
println "Workflow for project : $params.projectName"
println "Workflow description : $workflow.manifest.description"
println "Workflow gitLab URL : $workflow.manifest.homePage"
println "Workflow authors : $workflow.manifest.author"
println "Workflow source code : $workflow.projectDir"
println "Cmd line: $workflow.commandLine"
println "Workflow working/temp directory : $workflow.workDir"
println "Workflow output/publish directory : $params.outdir"
println "Workflow configuration file : $workflow.configFiles"
println "Directory containing raw data : $params.rawdata_dir"

Channel.fromFilePairs(params.rawdata_dir, checkIfExists:true, flat:true)
  .ifEmpty { exit 1, error "${params.rawdata_dir} is empty - no read files supplied" }
  .into { fastqc_reads ; unicyler_reads }

/* Check quality of raw reads */

process fastqc_quality_process {
  beforeScript "${params.fastqc_env}"

  // publishDir "${params.outdir}/${params.quality_check_dirname}/fastqc_output" , mode: 'copy', pattern: '*_fastqc.{zip,html}'
  // publishDir "${params.outdir}/${params.quality_check_dirname}/fastqc_output" , mode: 'copy', pattern: 'fastqc_process_*.log'
  // publishDir "${params.outdir}/${params.publish_dirname}/00_quality_check/fastqc_output", mode: 'copy', pattern: '*_fastqc.{zip,html}'
  // publishDir "${params.outdir}/${params.publish_dirname}", mode: 'copy', pattern : 'completecmd_fastqc_*', saveAs : { completecmd_fastqc -> "cmd/00_a_${task.process}_complete.sh" }

  input :
    set id, file(R1), file(R2) from fastqc_reads

    output :
      file "*_fastqc.{zip,html}" into fastqc_results
      // file 'completecmd_fastqc_*' into completecmd_fastqc
      file 'fastqc_process_*.log' into fastqc_process_log

    //Run only if process is activated in params.config file
    // when :
    //    params.quality_check_enable

    script :
    """
    fastqc ${R1} ${R2} -t ${task.cpus} >& fastqc_process_!{id}.log 2>&1
    """
}

process multiqc_quality_process {
    beforeScript "${params.multiqc_env}"

    // publishDir "${params.outdir}/${params.quality_merge_dirname}/multiqc_output" , mode: 'copy', pattern: 'multiqc_report.html'
    // publishDir "${params.outdir}/${params.publish_dirname}/00_quality_check/multiqc_output", mode: 'copy', pattern: 'multiqc_report.html'
    // publishDir "${params.outdir}/${params.publish_dirname}", mode: 'copy', pattern : 'completecmd_multiqc', saveAs : { completecmd_multiqc -> "cmd/00_b_${task.process}_complete.sh" }

    input :
      file ('fastqc/*') from fastqc_results.collect().ifEmpty([])

    output :
      file "multiqc_report.html" into multiqc_report
      file "multiqc_data"
      // file 'completecmd_multiqc' into completecmd_multiqc

    //Run only if process is activated in params.config file
    // when :
        // params.quality_check_enable

     script :
     """
     multiqc . >& multiqc_process.log 2>&1
     """
}

/* Genome assembly using Unicycler */

// process assembly_process {
    // tag "$id"
    // beforeScript "${params.unicycler_env}"

    // publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern: 'assembly.gfa'
    // publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : 'assembly.fasta'
    // publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : 'unicycler.log'
    // publishDir "${params.outdir}/${params.publish_dirname}", mode: 'copy', pattern : 'completecmd_filtering_*', saveAs : { completecmd_filtering -> "cmd/01_${task.process}_complete.sh" }
    // publishDir "${params.outdir}/${params.publish_dirname}/01_filtering_process", mode: 'copy', pattern : '*_trimming_report.txt'

    // input :
    //     set genome_name, file(read1) , file(read2) from unicyler_reads
    //
    // output :
    //   file "assembly.gfa" into assembly_gfa
    //   file "assembly.fasta" into assembly_fasta
    //   file "unicycler.log" into assembly_log

    //Run only if process is activated in params.config file
    // when :
    //    params.filtering_process_enable

    // shell :
    // """
    // unicycler ${read1} ${read2} --min_fasta_length ${params.assembly_process.min_fasta_length} -t ${task.cpus} --keep 0 --mode normal >& filtering_process_!{id}.log 2>&1
    // """
// }

/* Detection and removal of potential vectors, adpatators or contamination */
