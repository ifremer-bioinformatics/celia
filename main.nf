#!/usr/bin/env nextflow

/*
========================================================================================
                          CELIA: automatiC gEnome assembLy marIne prokAryotes
========================================================================================
 CELIA Analysis Pipeline.
 #### Homepage / Documentation
 https://gitlab.ifremer.fr/bioinfo/CELIA
----------------------------------------------------------------------------------------
*/

def helpMessage() {
  // Add to this help message with new command line parameters
  log.info SeBiMERHeader()
  log.info"""
  Usage:

  The typical command for running the pipeline after filling the conf/base.config file is as follows :
    nextflow run main.nf

    Mandatory arguments:
    --rawdata_dir [path]                    Path to input directory with raaw data files

    Other options:
    --outdir [path]                         The output directory where the results will be saved
    -w/--work-dir                           The temporary directory where intermediate data will be saved
    -name [str]                             Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    --projectName [str]                     Name of the project being analyzed

    Quality checking:
    --quality_check_enable [bool]           Quality checking step | FastQC, MultiQC (default: true)
    --quality_check_post_enable [bool]      Quality checking step of genome assembly | BUSCO, Bowtie2, MosDepth (default: true)
    --odb_path [path]                       Path to all BUSCO ODB databases
    --odb_name [str]                        Specify the name of the BUSCO lineage to be used

    Genome assembly:
    --min_ctg_length [int]                  Exclude contigs from the FASTA file which are shorter than this length (default: 500)
    --assembly_mode [str]                   Unicycler graph mode: normal, conservative, bold (default: normal)

    Vectors detection:
    --univec_db [path]                      UniVec database location
    --outfmt [int/str]                      Blast output format (default: 6)

    ANI calculation:
    --ani_enable [bool]                     Compute ANI scores (default: true)
    --ani_db [str]                          A file containing list of reference genome files, one genome per line

    Plasmid detection:
    --plasmid_enable [bool]                 Identification and characterization of bacterial plasmid contigs (default: true)
    --platon_db [path]                      Platon plasmid database
    --platon_mode [str]                     Applied filter mode. Sensitivity, specificity or accuracy (default: accuracy)

    Annotation:
    --annotation_enable [bool]              Make structural and functional annotation (default: true)
    --evalue                                Similarity e-value cut-off (default: 1e-09)
    --center                                Sequencing centre ID. (default: Ifremer)
    --kingdom                               Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
    --gcode                                 Genetic code / Translation table (set if --kingdom is set) (default '11')

  """.stripIndent()
}

// Show help message
if (params.help) {
  helpMessage()
  exit 0
}

/*
* SET UP CONFIGURATION VARIABLES
*/

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
  custom_runName = workflow.runName
}

Channel.fromFilePairs(params.rawdata_dir, checkIfExists:true, flat:true)
  .ifEmpty { exit 1, error "${params.rawdata_dir} is empty - no read files supplied" }
  .into { fastqc_reads ; unicycler_reads ; bowtie2_reads }

/*
* PIPELINE INFO
*/
// Header log info
log.info SeBiMERHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Project Name']     = params.projectName
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName

// if (params.email || params.email_on_fail) {
//   summary['E-mail Address']    = params.email
//   summary['E-mail on failure'] = params.email_on_fail
// }

log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

/*
* STEP 0 - Check quality of raw reads
*/

process fastqc {
  label 'fastqc'
  beforeScript "${params.fastqc_env}"

  publishDir "${params.outdir}/${params.quality_check_dirname}" , mode: 'copy', pattern : '*_fastqc.{zip,html}'
  publishDir "${params.outdir}/${params.quality_check_dirname}" , mode: 'copy', pattern : 'fastqc_process_*.log'

  input:
    set id, file(R1), file(R2) from fastqc_reads

  output:
    file "*_fastqc.{zip,html}" into fastqc_results
    file 'fastqc_process_*.log' into fastqc_process_log

  when:
    params.quality_check_enable

  shell:
  """
  fastqc ${R1} ${R2} -t ${task.cpus} >& fastqc_process_!{id}.log 2>&1
  """
}

process multiqc {
  beforeScript "${params.multiqc_env}"

  publishDir "${params.outdir}/${params.quality_merge_dirname}" , mode: 'copy', pattern : 'multiqc_report.html'

  input:
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])

  output:
    file "multiqc_report.html" into multiqc_report
    file "multiqc_data"

  when:
    params.quality_check_enable

  shell:
  """
  multiqc . >& multiqc_process.log 2>&1
  """
}

/*
* STEP 1 - Genome assembly using Unicycler
*/

process unicycler {
  label 'assembly'
  beforeScript "${params.unicycler_env}"

  publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : "${assembly_name}/*.log" , saveAs : { unicycler_log -> "${assembly_name}/unicycler.log" }
  publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : "${assembly_name}/*.gfa" , saveAs : { unicycler_gfa -> "${assembly_name}/assembly.gfa" }
  publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : "${assembly_name}/*.fasta" , saveAs : { unicycler_fasta -> "${assembly_name}/assembly.fasta" }

  input:
    set assembly_name, file(read1) , file(read2) from unicycler_reads

  output:
    set assembly_name, file("${assembly_name}/*.fasta") into unicycler_fasta
    set assembly_name, file("${assembly_name}/*.gfa") into unicycler_gfa
    file "${assembly_name}/*.log" into unicycler_log

  shell:
  """
  unicycler -1 ${read1} -2 ${read2} -o ${assembly_name} --min_fasta_length ${params.min_ctg_length} -t ${task.cpus} --keep 0 --mode ${params.assembly_mode} >& assembly_process_!{assembly_name}.log 2>&1
  """
}

// process shovill {
//   beforeScript "${params.shovill_env}"
//
//   publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : "${assembly_name}/*.log" , saveAs : { shovill_log -> "${assembly_name}/shovill.log" }
//   publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : "${assembly_name}/*.gfa" , saveAs : { shovill_gfa -> "${assembly_name}/contigs.gfa" }
//   publishDir "${params.outdir}/${params.assembly_dirname}", mode: 'copy', pattern : "${assembly_name}/*.fa" , saveAs : { shovill_fasta -> "${assembly_name}/contigs.fasta" }
//
//   input :
//     set assembly_name, file(read1) , file(read2) from unicycler_reads
//
//   output :
//     set assembly_name, file("${assembly_name}/*.fa") into shovill_fasta
//     file "${assembly_name}/*.gfa" into shovill_gfa
//     file "${assembly_name}/*.log" into shovill_log
//
//   shell :
//   """
//   shovill --R1 ${read1} --R2 ${read2} --depth 0 --outdir ${assembly_name}/ --minlen ${params.shovill.min_fasta_length} --tmpdir ${TMPDIR} --keepfiles --cpus ${task.cpus} --ram ${task.memory.toGiga()} >& assembly_process_!{assembly_name}.log 2>&1
//   """
// }

process bandage {
  beforeScript "${params.bandage_env}"

  publishDir "${params.outdir}/${params.assembly_dirname}" , mode: 'copy', saveAs : { bandage_plot -> "${assembly_name}/${assembly_name}.svg" }

  input:
    set assembly_name, file(gfa) from unicycler_gfa

  output:
    file("*.svg")

  shell:
  """
  Bandage image ${gfa} ${assembly_name}.svg >& bandage.log 2>&1
  """
}

/*
* STEP 2 - Detection and removal of potential vectors, adpatators or contaminants
*/

process blast {
  beforeScript "${params.blast_env}"

  publishDir "${params.outdir}/${params.contamination_check_dirname}" , mode: 'copy', pattern : '*.m6'

  input:
    set assembly_name, file(fasta) from unicycler_fasta

  output:
    set assembly_name, file("*.m6"), file(fasta) into univec_blast_fasta

  shell:
  """
  blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -db ${params.univec_db} -query ${fasta} -out ${assembly_name}.m6 -outfmt ${params.outfmt} >& blast_univec_process_!{assembly_name}.log 2>&1
  """
}

process remoVecSec {
  beforeScript "${params.biopython_env}"

  publishDir "${params.outdir}/${params.contamination_rm_dirname}" , mode: 'copy', pattern : '*.clean.fasta'

  input:
    set assembly_name, file(blast), file(fasta) from univec_blast_fasta

  output:
    set assembly_name, file("*.clean.fasta") into vecscreen_fasta

  shell:
  """
  vecscreen.py -i ${blast} -f ${fasta} -o ${assembly_name}.clean.fasta >& vecscreen.log 2>&1
  """

}

/* Create channels from decontaminated fasta to genome metrics/quality steps */
vecscreen_fasta.into {
  fasta_busco
  fasta_bowtie2
  fasta_ani
  fasta_prokka
  fasta_platon
}

/*
* STEP 3 - Assembly quality and metrics
*/

process bowtie2 {
  label 'bowtie2'
  beforeScript "${params.bowtie2_env}"

  publishDir "${params.outdir}/${params.assembly_mapping_dirname}", mode: 'copy', pattern : '*.bam'
  publishDir "${params.outdir}/${params.assembly_mapping_dirname}", mode: 'copy', pattern : '*.bai'
  publishDir "${params.outdir}/${params.assembly_mapping_dirname}", mode: 'copy', pattern : '*.bowtie2-mapping.log'

  input:
    set id, file(fasta), file(read1), file(read2) from fasta_bowtie2.join(bowtie2_reads)

  output:
    set id, file("*.bam"), file("*.bai") into bowtie2_bam
    file "*.bowtie2-mapping.log" into bowtie2_logs

  when:
    params.quality_check_post_enable

  shell:
  """
  bowtie2-build ${fasta} ${id} -p 4 > ${id}.bowtie2-build.log 2>&1

  bowtie2 -t -p ${task.cpus} --no-unal -1 ${read1} -2 ${read2} --sensitive -x ${id} 2> ${id}.bowtie2-mapping.log | \
  samtools sort -@ 1 -m 4G -o ${id}.bam -O bam -T tmp - 2> ${id}.samtools-sort.log 

  samtools index ${id}.bam >& ${id}.samtools-index.log 2>&1
  """
}

process mosDepth {
  beforeScript "${params.mosdepth_env}"

  publishDir "${params.outdir}/${params.assembly_coverage_dirname}", mode: 'copy', pattern : '*.mosdepth.global.dist.txt'
  publishDir "${params.outdir}/${params.assembly_coverage_dirname}", mode: 'copy', pattern : '*.mosdepth.summary.txt'
  publishDir "${params.outdir}/${params.assembly_coverage_dirname}", mode: 'copy', pattern : '*.per-base.bed.gz'

  input:
    set assembly_name, file(bam), file(bai) from bowtie2_bam

  output:

    file "*.mosdepth.global.dist.txt" into mosdepth_global
    file "*.mosdepth.summary.txt" into mosdepth_summary
    file "*.per-base.bed.gz" into mosdepth_bed

  when:
    params.quality_check_post_enable

  shell:
  """
  mosdepth ${assembly_name} ${bam} >& mosdepth.log 2>&1
  """
}

process busco {
  label 'busco'
  beforeScript "${params.busco_env}"

  publishDir "${params.outdir}/${params.assembly_completness_dirname}", mode: 'copy', pattern : "${assembly_name}/short_summary*"
  publishDir "${params.outdir}/${params.assembly_completness_dirname}", mode: 'copy', pattern : "${assembly_name}/run_*/full_table.tsv"
  publishDir "${params.outdir}/${params.assembly_completness_dirname}", mode: 'copy', pattern : "${assembly_name}/run_*/missing_busco_list.tsv"

  input:
    set assembly_name, file(fasta) from fasta_busco

  output:
    file "${assembly_name}/short_summary*" into busco_short_summary
    file "${assembly_name}/run_*/full_table.tsv" into busco_full_summary
    file "${assembly_name}/run_*/missing_busco_list.tsv" into busco_missing_list

  when:
    params.quality_check_post_enable

  shell:
  """
  busco -c ${task.cpus} --force --offline -m genome -i ${fasta} -o ${assembly_name} -l ${params.odb_path}/${params.odb_name} >& busco.log 2>&1
  """
}

/*
* STEP 4 - ANI
*/

process fastANI {
  beforeScript "${params.fastani_env}"

  publishDir "${params.outdir}/${params.wgs_similarity_ANI_dirname}", mode: 'copy', pattern : "*.ani"

  input:
    set assembly_name, file(fasta) from fasta_ani

  output:
    file "*.ani" into fastANI_summary

  when:
    params.ani_enable

  shell:
  """
  fastANI -q ${fasta} --rl ${params.ani_db} -o ${assembly_name}.ani >& fastANI.log 2>&1
  """
}

/*
* STEP 5 - Plasmid detection
*/

process platon {
  label 'platon'
  beforeScript "${params.platon_env}"

  publishDir "${params.outdir}/${params.plasmid_detection_dirname}", mode: 'copy', pattern : "*.fasta"
  publishDir "${params.outdir}/${params.plasmid_detection_dirname}", mode: 'copy', pattern : "*.tsv"
  publishDir "${params.outdir}/${params.plasmid_detection_dirname}", mode: 'copy', pattern : "*.json"

  input:
    set assembly_name, file(fasta) from fasta_platon

  output:
    file "*.fasta" //into fastANI_summary
    file "*.tsv" //into fastANI_summary
    file "*.json" //into fastANI_summary

  when:
    params.plasmid_enable

  shell:
  """
  platon --db ${params.platon_db} --mode ${params.platon_mode} --threads ${task.cpus} ${fasta} >& platon.log 2>&1
  """
}

/*
* STEP 6 - Structural and functional annotation
*/

process prokka {
  label 'prokka'
  beforeScript "${params.prokka_env}"

  publishDir "${params.outdir}/${params.gene_prediction_dirname}", mode: 'copy'

  input:
    set assembly_name, file(fasta) from fasta_prokka

  output:
    file "prokka/*" into prokka_annotation

  when:
    params.annotation_enable

  shell:
  """
  prokka --outdir prokka --prefix ${assembly_name} --centre ${params.center} --compliant --addgenes --gffver 3 --kingdom ${params.kingdom} --gcode ${params.gcode} --mincontiglen ${params.min_ctg_length} --cpus ${task.cpus} ${fasta} >& prokka.log 2>&1
  """
}

/* Other functions */
def SeBiMERHeader() {
    // Log colors ANSI codes
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_cyan}--------------------------------------------------${c_reset}-
    ${c_blue}    __  __  __  .       __  __  ${c_reset}
    ${c_blue}   \\   |_  |__) | |\\/| |_  |__)  ${c_reset}
    ${c_blue}  __\\  |__ |__) | |  | |__ |  \\  ${c_reset}
                                            ${c_reset}
    ${c_yellow}  CELIA: automatiC gEnome assembLy marIne prokAryotes${c_reset}
    -${c_cyan}--------------------------------------------------${c_reset}-
    """.stripIndent()
}
