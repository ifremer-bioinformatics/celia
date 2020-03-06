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

// process unicycler {
//   beforeScript "${params.unicycler_env}"
//
//   publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : "${assembly_name}/*.log" , saveAs : { unicycler_log -> "${assembly_name}/unicycler.log" }
//   publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : "${assembly_name}/*.gfa" , saveAs : { unicycler_gfa -> "${assembly_name}/assembly.gfa" }
//   publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : "${assembly_name}/*.fasta" , saveAs : { unicycler_fasta -> "${assembly_name}/assembly.fasta" }
//
//   input :
//     set assembly_name, file(read1) , file(read2) from unicycler_reads
//
//   output :
//     set assembly_name, file("${assembly_name}/*.fasta") into unicycler_fasta
//     file "${assembly_name}/*.gfa" into unicycler_gfa
//     file "${assembly_name}/*.log" into unicycler_log
//
//   shell :
//   """
//   unicycler -1 ${read1} -2 ${read2} -o ${assembly_name} --min_fasta_length ${params.unicycler.min_fasta_length} -t ${task.cpus} --keep 0 --mode normal >& assembly_process_!{assembly_name}.log 2>&1
//   """
// }

process shovill {
  beforeScript "${params.shovill_env}"

  publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : "${assembly_name}/*.log" , saveAs : { shovill_log -> "${assembly_name}/shovill.log" }
  publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : "${assembly_name}/*.gfa" , saveAs : { shovill_gfa -> "${assembly_name}/contigs.gfa" }
  publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : "${assembly_name}/*.fasta" , saveAs : { shovill_fasta -> "${assembly_name}/contigs.fasta" }

  input :
    set assembly_name, file(read1) , file(read2) from unicycler_reads

  output :
    set assembly_name, file("${assembly_name}/*.fasta") into shovill_fasta
    file "${assembly_name}/*.gfa" into shovill_gfa
    file "${assembly_name}/*.log" into shovill_log

  shell :
  """
  shovill --R1 ${read1} --R2 ${read2} --depth 0 --outdir ${assembly_name}/ --minlen ${params.shovill.min_fasta_length} --tmpdir ${TMPDIR} --keepfiles --cpus ${task.cpus} --ram ${task.memory.toGiga()} >& assembly_process_!{assembly_name}.log 2>&1
  """
}


/* Detection and removal of potential vectors, adpatators or contaminants */

process blast {
  beforeScript "${params.blast_env}"

  publishDir "${params.outdir}/${params.contamination_check_dirname}" , mode: 'copy', pattern : '*.m6'

  input :
    // set assembly_name, file(fasta) from unicycler_fasta
    set assembly_name, file(fasta) from shovill_fasta

  output :
    set assembly_name, file("*.m6"), file(fasta) into univec_blast_fasta

  shell :
  """
  blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -db ${params.blast.univec_db} -query ${fasta} -out ${assembly_name}.m6 -outfmt ${params.blast.outfmt} >& blast_univec_process_!{assembly_name}.log 2>&1
  """
}

process remoVecSec {
  beforeScript "${params.biopython_env}"

  publishDir "${params.outdir}/${params.contamination_rm_dirname}" , mode: 'copy', pattern : '*.clean.fasta'

  input :
    set assembly_name, file(blast), file(fasta) from univec_blast_fasta

  output :
    set assembly_name, file("*.clean.fasta") into vecscreen_fasta

  shell
  """
  ${baseDir}/lib/vecscreen.py -i ${blast} -f ${fasta} -o ${assembly_name}.clean.fasta
  """

}

/* Create channels from decontaminated fasta to genome metrics/quality steps */
vecscreen_fasta.into {
  fasta_busco
  fasta_bowtie2
  fasta_ani
  fasta_prokka
}

/* Quality and metrics of assembly */

process bowtie2 {
  beforeScript "${params.bowtie2_env}"

  publishDir "${params.outdir}/${params.assembly_mapping_dirname}", mode: 'copy', pattern : '*.bam'
  publishDir "${params.outdir}/${params.assembly_mapping_dirname}", mode: 'copy', pattern : '*.bai'
  publishDir "${params.outdir}/${params.assembly_mapping_dirname}", mode: 'copy', pattern : '*.bowtie2-mapping.log'

  input :
    set assembly_name, file(fasta) from fasta_bowtie2
    set reads_id, file(read1) , file(read2) from bowtie2_reads
  output :
    set assembly_name, file("*.bam"), file("*.bai") into bowtie2_bam
    file "*.bowtie2-mapping.log" into bowtie2_logs

  shell :
  """
  bowtie2-build ${fasta} ${assembly_name} -p 4 > ${assembly_name}.bowtie2-build.log 2>&1

  bowtie2 --phred33 --non-deterministic -t -p ${task.cpus} --fr -1 ${read1} -2 ${read2} --sensitive -x ${assembly_name} 2> ${assembly_name}.bowtie2-mapping.log | \
  samtools sort -@ 1 -m 4G -O bam -l 0 -T tmp - 2> ${assembly_name}.samtools-sort.log | \
  samtools view -b -o ${assembly_name}.bam - 2> ${assembly_name}.samtools-view.log

  samtools index ${assembly_name}.bam >& ${assembly_name}.samtools-index.log 2>&1
  """
}

process mosDepth {
  beforeScript "${params.mosdepth_env}"

  publishDir "${params.outdir}/${params.assembly_coverage_dirname}", mode: 'copy', pattern : '*.mosdepth.global.dist.txt'
  publishDir "${params.outdir}/${params.assembly_coverage_dirname}", mode: 'copy', pattern : '*.mosdepth.summary.txt'
  publishDir "${params.outdir}/${params.assembly_coverage_dirname}", mode: 'copy', pattern : '*.per-base.bed.gz'

  input :
    set assembly_name, file(bam), file(bai) from bowtie2_bam

  output :

    file "*.mosdepth.global.dist.txt" into mosdepth_global
    file "*.mosdepth.summary.txt" into mosdepth_summary
    file "*.per-base.bed.gz" into mosdepth_bed

  shell :
  """
  mosdepth ${assembly_name} ${bam} >& mosdepth.log 2>&1
  """
}

process busco {
  beforeScript "${params.busco_env}"

  publishDir "${params.outdir}/${params.assembly_completness_dirname}", mode: 'copy', pattern : "${assembly_name}/short_summary*"
  publishDir "${params.outdir}/${params.assembly_completness_dirname}", mode: 'copy', pattern : "${assembly_name}/run_*/full_table.tsv"
  publishDir "${params.outdir}/${params.assembly_completness_dirname}", mode: 'copy', pattern : "${assembly_name}/run_*/missing_busco_list.tsv"

  input :
    set assembly_name, file(fasta) from fasta_busco

  output :
    file "${assembly_name}/short_summary*" into busco_short_summary
    file "${assembly_name}/run_*/full_table.tsv" into busco_full_summary
    file "${assembly_name}/run_*/missing_busco_list.tsv" into busco_missing_list

  shell :
  """
  busco -c ${task.cpus} --force --offline -m genome -i ${fasta} -o ${assembly_name} -l ${params.busco.db_path}/${params.busco.db_name} >& busco.log 2>&1
  """
}

process fastANI {
  beforeScript "${params.fastani_env}"

  publishDir "${params.outdir}/${params.wgs_similarity_ANI_dirname}", mode: 'copy', pattern : "*.ani"

  input :
    set assembly_name, file(fasta) from fasta_ani

  output :
    file "*.ani" into fastANI_summary

  shell :
  """
  fastANI -q ${fasta} --rl ${params.fastani.db} -o ${assembly_name}.ani >& fastANI.log 2>&1
  """
}

// process prokka {
//   beforeScript "${params.prokka_env}"
//
//   // publishDir "${params.outdir}/${params.gene_prediction_dirname}", mode: 'copy', pattern : "prokka/${assembly_name}*"
//   publishDir "${params.outdir}/${params.gene_prediction_dirname}", mode: 'copy', pattern : "prokka/${assembly_name}*" , saveAs : { prokka_annotation -> "${assembly_name}*" }
//
//   input :
//     set assembly_name, file(fasta) from fasta_prokka
//
//   output :
//     file "prokka/${assembly_name}*" into prokka_annotation
//
//   shell :
//   """
//   prokka --outdir prokka --prefix ${assembly_name} --centre Ifremer --compliant --addgenes --gffver 3 --kingdom ${params.prokka.kingdom} --gcode ${params.prokka.gcode} --mincontiglen ${params.prokka.min_contig_length} --cpus ${task.cpus} ${fasta} >& prokka.log 2>&1
//   """
// }
