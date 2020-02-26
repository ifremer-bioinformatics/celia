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
  .into { fastqc_reads ; unicycler_reads ; bowtie2_reads }

/* Check quality of raw reads */

process fastqc {
  beforeScript "${params.fastqc_env}"

  publishDir "${params.outdir}/${params.quality_check_dirname}" , mode: 'copy', pattern : '*_fastqc.{zip,html}'
  publishDir "${params.outdir}/${params.quality_check_dirname}" , mode: 'copy', pattern : 'fastqc_process_*.log'

  input :
    set id, file(R1), file(R2) from fastqc_reads

  output :
    file "*_fastqc.{zip,html}" into fastqc_results
    file 'fastqc_process_*.log' into fastqc_process_log

  shell :
  """
  fastqc ${R1} ${R2} -t ${task.cpus} >& fastqc_process_!{id}.log 2>&1
  """
}

process multiqc {
  beforeScript "${params.multiqc_env}"

  publishDir "${params.outdir}/${params.quality_merge_dirname}" , mode: 'copy', pattern : 'multiqc_report.html'

  input :
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])

  output :
    file "multiqc_report.html" into multiqc_report
    file "multiqc_data"

   shell :
   """
   multiqc . >& multiqc_process.log 2>&1
   """
}

/* Genome assembly using Unicycler */

process unicycler {
  beforeScript "${params.unicycler_env}"

  publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : 'assembly/*.gfa' , saveAs : { unicycler_gfa -> "${assembly_name}/assembly.gfa" }
  publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : 'assembly/*.fasta' , saveAs : { unicycler_fasta -> "${assembly_name}/assembly.fasta" }
  publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : 'assembly/*.log' , saveAs : { unicycler_log -> "${assembly_name}/unicycler.log" }
  // publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : "${assembly_name}/*.log" , saveAs : { unicycler_log -> "${assembly_name}/unicycler.log" }


  input :
    set assembly_name, file(read1) , file(read2) from unicycler_reads

  output :
    set assembly_name, file("assembly/assembly.fasta") into unicycler_fasta
    file "assembly/assembly.gfa" into unicycler_gfa
    file "assembly/unicycler.log" into unicycler_log
    // file "${assembly_name}/unicycler.log" into unicycler_log

  shell :
  """
  unicycler -1 ${read1} -2 ${read2} -o assembly/ --min_fasta_length ${params.unicycler.min_fasta_length} -t ${task.cpus} --keep 0 --mode normal >& assembly_process_!{assembly_name}.log 2>&1
  #unicycler -1 ${read1} -2 ${read2} -o ${assembly_name} --min_fasta_length ${params.unicycler.min_fasta_length} -t ${task.cpus} --keep 0 --mode normal >& assembly_process_!{assembly_name}.log 2>&1
  """
}

/* Detection and removal of potential vectors, adpatators or contamination */

// process blast {
//   beforeScript "${params.blast_env}"
//
//   publishDir "${params.outdir}/${params.contamination_check_dirname}" , mode: 'copy', pattern : '*.m6'
//
//   input :
//     set assembly_name, file(fasta) from unicycler_fasta
//
//   output :
//     file "*.m6" into genome_univec
//
//   shell :
//   """
//   blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -db ${params.blast.univec_db} -query ${fasta} -out ${assembly_name}.m6 -outfmt 6 >& blast_univec_process_!{genome_name}.log 2>&1
//   """
// }

// process remoVecSec {
//   beforeScript "${params.blast_env}"
//
//   input :
//     set val(assembly_id), file(fasta) from unicycler_fasta
//
//   output :
//     // file "${genome_name}.m6" into genome_univec
//   shell
//   """
//   remower.py [-h] --genomefile GENOMEFILE [--dbvec DBVEC] [--dbmito DBMITO] [--dbcont DBCONT] [--dist DIST]
//   """
//
// }

/* Quality and metrics of assembly */

// process bowtie2 {
//   beforeScript "${params.bowtie2_env}"
//
//   publishDir "${params.outdir}/${params.assembly_mapping_dirname}", mode: 'copy', pattern : '*.bam'
//   publishDir "${params.outdir}/${params.assembly_mapping_dirname}", mode: 'copy', pattern : '*.bai'
//
//   input :
//     set assembly_name, file(fasta) from unicycler_fasta
//     set reads_id, file(read1) , file(read2) from bowtie2_reads
//     // attention au merge des ids
//   output :
//     file "*.bam" into bowtie2_bam
//     file "*.bai" into bowtie2_index
//
//   shell :
//   """
//   bowtie2-build ${fasta} ${assembly_name} -p 4 > bowtie2-build.log 2>&1
//
//   bowtie2 --phred33 --non-deterministic -t -p ${task.cpus} --fr \
//   -1 ${read1} -2 ${read2} --sensitive -x ${assembly_name} | \
//   samtools sort -@ 1 -m 4G -O bam -l 0 -T tmp - | \
//   samtools view -b -o ${assembly_name}.bam - >& bowtie2-mapping.log 2>&1
//
//   samtools index ${assembly_name}.bam >& samtools-index.log 2>&1
//   """
// }

// process mosDepth {
//   beforeScript "${params.mosdepth_env}"
//
//   publishDir "${params.outdir}/${params.assembly_coverage_dirname}", mode: 'copy', pattern : '*.mosdepth.global.dist.txt'
//   publishDir "${params.outdir}/${params.assembly_coverage_dirname}", mode: 'copy', pattern : '*.mosdepth.summary.txt'
//   publishDir "${params.outdir}/${params.assembly_coverage_dirname}", mode: 'copy', pattern : '*.per-base.bed.gz'
//
//   input :
//     set assembly_id, file(bam) from bowtie2_bam
//
//   output :
//     file "*.mosdepth.global.dist.txt" into mosdepth_global
//     file "*.mosdepth.summary.txt" into mosdepth_summary
//     file "*.per-base.bed.gz" into mosdepth_bed
//
//   shell :
//   """
//   mosdepth ${assembly_id} ${bam} >& mosdepth.log 2>&1
//   """
// }
//
// process busco {
//   beforeScript "${params.busco_env}"
//
//   publishDir "${params.outdir}/${params.assembly_completness_dirname}", mode: 'copy', pattern : '${assembly_id}/short_summary*'
//   publishDir "${params.outdir}/${params.assembly_completness_dirname}", mode: 'copy', pattern : '${assembly_id}/run_*/full_table.csv'
//   publishDir "${params.outdir}/${params.assembly_completness_dirname}", mode: 'copy', pattern : '${assembly_id}/run_*/missing_busco_list.tsv'
//
//   input :
//     set assembly_id, file(fasta) from unicycler_fasta
//
//   output :
//     file "${assembly_id}/short_summary*" into busco_short_summary
//     file "${assembly_id}/run_*/full_table.csv" into busco_full_summary
//     file "${assembly_id}/run_*/missing_busco_list.tsv" into busco_missing_list
//
//   shell :
//   """
//   busco -c ${task.cpus} --force --offline -m genome -i ${fasta} -o ${assembly_id} -l ${params.busco.db_path}/${params.busco.db_name} >& busco.log 2>&1
//   """
// }
//
// process fastANI {
//   beforeScript "${params.fastani_env}"
//
//   input :
//     set assembly_id, file(fasta) from unicycler_fasta
//
//   output :
//     file "${assembly_id}.ani" into fastANI_summary
//
//   shell :
//   """
//   fastANI -q ${fasta} --rl ${refGenome} -o ${assembly_id}.ani >& fastANI.log 2>&1
//   """
// }
