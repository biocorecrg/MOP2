#!/usr/bin/env nextflow

nextflow.preview.dsl=2

/* 
 * Define the pipeline parameters
 *
 */

// Pipeline version
version = '2.0'

params.help            = false
params.resume          = false

log.info """

╔╦╗┌─┐┌─┐┌┬┐┌─┐┬─┐  ┌─┐┌─┐  ╔═╗╔═╗╦═╗╔═╗╔═╗
║║║├─┤└─┐ │ ├┤ ├┬┘  │ │├┤   ╠═╝║ ║╠╦╝║╣ ╚═╗
╩ ╩┴ ┴└─┘ ┴ └─┘┴└─  └─┘└    ╩  ╚═╝╩╚═╚═╝╚═╝
                                                                                       
====================================================
BIOCORE@CRG Preprocessing of Nanopore direct RNA - N F  ~  version ${version}
====================================================

kit                       : ${params.kit}
flowcell                  : ${params.flowcell}
fast5                     : ${params.fast5}
reference                 : ${params.reference}
annotation                : ${params.annotation}

ref_type                  : ${params.ref_type}
seq_type                  : ${params.seq_type}

output                    : ${params.output}
qualityqc                 : ${params.qualityqc}
granularity               : ${params.granularity}

basecaller                : ${params.basecaller}
basecaller_opt            : ${params.basecaller_opt}
GPU                       : ${params.GPU}
demultiplexing            : ${params.demultiplexing} 
demultiplexing_opt        : ${params.demultiplexing_opt} 
demulti_fast5		      : ${params.demulti_fast5}

filter                    : ${params.filter}
filter_opt                : ${params.filter_opt}
mapper                    : ${params.mapper}
mapper_opt                : ${params.mapper_opt}
map_type                  : ${params.map_type}

counter                   : ${params.counter}
counter_opt               : ${params.counter_opt}

downsampling			  : ${params.downsampling}

variant_caller            : ${params.variant_caller}
variant_opt               : ${params.variant_opt}

email                     : ${params.email}
"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"
granularity = 1000000000
if (params.granularity == "") granularity = 1000000000

// check multi5 and GPU usage. GPU maybe can be removed as param if there is a way to detect it
if (params.GPU != "ON" && params.GPU != "OFF") exit 1, "Please specify ON or OFF in GPU processors are available"

if (params.map_type != "unspliced" && params.map_type != "spliced") exit 1, "Mapping type NOT supported! Please choose either 'spliced' or 'unspliced'"

// check input files
reference = file(params.reference)
if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"
config_report = file("$baseDir/config.yaml")
if( !config_report.exists() ) exit 1, "Missing config.yaml file!"
logo = file("$baseDir/../docs/logo_small.png")
deeplexicon_folder = file("$baseDir/deeplexicon/")


basecaller   		= params.basecaller
basecaller_opt  	= params.basecaller_opt
demultiplexer 		= params.demultiplexing
demultiplexer_opt   = params.demultiplexing_opt
mapper      		= params.mapper
mapper_opt   		= params.mapper_opt
counter_opt   		= params.counter_opt 
granularity			= params.granularity 
gpu				    = params.GPU

if (params.basecaller != "guppy" && gpu != "OFF") {
	log.info "GPU param will be set to OFF since the basecaller is not guppy"
	gpu = "OFF"
}


// if you are using GPU analyse the whole dataset, otherwise make batch of 4,000 sequences if they are single fast5
// or single batches of multi fast5 sequences
//multi5_type_for_granularity.merge(fast5_4_granularity.collect()).map{
//    (params.GPU == "YES" ? params.granularity  : (it[0] == 0 ? it[1..-1].collate(4000) : it[1..-1].collate(1)) )
//}
// Output folders
outputFastq    = "${params.output}/fastq_files"
outputFast5    = "${params.output}/fast5_files"
outputQual     = "${params.output}/QC_files"
outputMultiQC  = "${params.output}/report"
outputMapping  = "${params.output}/alignment"
outputCRAM     = "${params.output}/cram_files"
outputCounts   = "${params.output}/counts"
outputVars     = "${params.output}/variants"
outputAssigned = "${params.output}/assigned"
outputReport   = file("${outputMultiQC}/multiqc_report.html")

/*
* move old multiQCreport
*/
if( outputReport.exists() ) {
  log.info "Moving old report to multiqc_report.html multiqc_report.html.old"
  outputReport.moveTo("${outputMultiQC}/multiqc_report.html.old")
}

/*
 * Creates the channels that emits fast5 files
 */
Channel
    .fromPath( params.fast5)                                             
    .ifEmpty { error "Cannot find any file matching: ${params.fast5}" }
    .set {fast5_files}


/*
* This is default value in case guppy will be used for RNA demultiplexing
*/
params.barcodekit = ""
if (demultiplexer == "") demultiplexer = "OFF"

if (params.granularity == "") params.granularity = 1000000000

if (params.ref_type == "genome") {
	if (params.annotation != "") {
		annotation = file(params.annotation)
		if( !annotation.exists() ) exit 1, "Missing annotation file: ${params.annotation}!"
	}
}
 
include TEST_FAST5 from "${baseDir}/modules/test_fast5"
include BASECALL_ONT from "${baseDir}/modules/modules.nf" addParams(GPU_OPTION: gpu)


fast5_files.map { 
    def filepath = file(it)
    def file_parts = "${filepath}".tokenize("/")
    def folder_name  = filepath[-2]
 	[folder_name, it]
}.groupTuple().set{ fast5_per_folder}



workflow flow1 {	
    fast5_type = TEST_FAST5(fast5_files)
	if (fast5_type == 0) granularity = 4000
	if (gpu == "YES") granularity = 1000000000
	
	fast5_per_folder.map{
		def folder_name = it[0]
		def buffer_files = it[1].flatten().collate(granularity)
		[folder_name, buffer_files]
	}.transpose().set{ fast5_4_analysis }

	if (demultiplexer == "OFF") BASECALL_ONT(fast5_4_analysis, fast5_type, basecaller, basecaller_opt, params.seq_type)
 	
}

workflow {
	flow1()
	 
}

workflow.onComplete {
    println "Pipeline BIOCORE@CRG Master of Pore completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}