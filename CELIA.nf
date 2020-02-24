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
  .into { fastqc_reads ; unicycler_reads }

/* Check quality of raw reads */

process fastqc {
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

process multiqc {
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

process unicycler {
  beforeScript "${params.unicycler_env}"

  // publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern: 'assembly.gfa'
  // publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : 'assembly.fasta'
  // publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : 'unicycler.log'
  // publishDir "${params.outdir}/${params.publish_dirname}", mode: 'copy', pattern : 'completecmd_filtering_*', saveAs : { completecmd_filtering -> "cmd/01_${task.process}_complete.sh" }
  // publishDir "${params.outdir}/${params.publish_dirname}/01_filtering_process", mode: 'copy', pattern : '*_trimming_report.txt'

  input :
    set genome_name, file(read1) , file(read2) from unicycler_reads

  output :
    file "${genome_name}/assembly.gfa" into unicycler_gfa
    file "${genome_name}/assembly.fasta" into unicycler_fasta
    file "unicycler.log" into unicycler_log

  // Run only if process is activated in params.config file
  // when :
  //    params.filtering_process_enable

  shell :
  """
  unicycler -1 ${read1} -2 ${read2} -o ${genome_name} --min_fasta_length ${params.unicycler.min_fasta_length} -t ${task.cpus} --keep 0 --mode normal >& assembly_process_!{genome_name}.log 2>&1
  """
}

/* Detection and removal of potential vectors, adpatators or contamination */

process blast {
  beforeScript "${params.blast_env}"

  input :
    set assembly_id, file(fasta) from unicycler_fasta

  // output :
  //   file "${genome_name}.m6" into genome_univec

  // Run only if process is activated in params.config file
  // when :
  //    params.filtering_process_enable

  shell :
  """
  #blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -db ${params.blast.univec_db} -query ${fasta} -out ${assembly_id}.m6 -outfmt 6 >& blast_univec_process_!{genome_name}.log 2>&1
  echo ${fasta}
  echo ${assembly_id}
  """
}

process remoVecSec {
  beforeScript "${params.blast_env}"

  input :
    set val(assembly_id), file(fasta) from unicycler_fasta

  output :
    // file "${genome_name}.m6" into genome_univec
  shell
  """
  remower.py [-h] --genomefile GENOMEFILE [--dbvec DBVEC] [--dbmito DBMITO] [--dbcont DBCONT] [--dist DIST]
  """

}

/* Quality and metrics of assembly */

process bowtie2 {
  beforeScript "${params.bowtie2_env}"

  input :
    set assembly_id, file(fasta) from assembly_fasta
    set reads_id, file(read1) , file(read2) from unicyler_reads

  output :
    file "*.bam" into bowtie2_bam
    // file "*.bai" into index_bowtie2

  script :
  """
  bowtie2-build ${fasta} ${assembly_id} -p 4 > bowtie2-build.log 2>&1

  #echo ${fasta%.*}
  echo ${assembly_id}

  bowtie2 --phred33 --non-deterministic -t -p ${task.cpus} --fr \
  -1 ${read1} -2 ${read2} --sensitive -x ${assembly_id} | \
  samtools sort -@ ${task.cpus} -m 4G -O bam -l 0 -T tmp - | \
  samtools view -b -o ${assembly_id}.bam - >& bowtie2-mapping.log 2>&1

  samtools index ${assembly_id}.bam
  """
}

process mosDepth {
  beforeScript "${params.mosdepth_env}"

  input :
    set name, file(bam) from bowtie2_bam

  output :
    file "*.mosdepth.global.dist.txt" into mosdepth_global
    file "*.mosdepth.summary.txt" into mosdepth_summary
    file "*.per-base.bed.gz" into mosdepth_bed

  script :
  """
  mosdepth ${name} ${bam} >& mosdepth.log 2>&1
  """
}

process busco {
  beforeScript "${params.busco_env}"

  input :
    set assembly_id, file(fasta) from unicycler_fasta

  output :
    file "${assembly_id}/short_summary*" into busco_short_summary
    file "${assembly_id}/run_*/full_table.csv" into busco_full_summary
    file "${assembly_id}/run_*/missing_busco_list.tsv" into busco_missing_list

  script :
  """
  busco -c ${task.cpus} --force --offline -m genome -i ${fasta} -o ${assembly_id} -l ${dbPath}/${dbName} >& busco.log 2>&1
  """
}

process fastANI {
  beforeScript "${params.fastani_env}"

  input :
    set assembly_id, file(fasta) from unicycler_fasta

  output :
    file "${assembly_id}.ani" into fastANI_summary

  script :
  """
  fastANI -q ${fasta} --rl ${refGenome} -o ${assembly_id}.ani >& fastANI.log 2>&1
  """
}
