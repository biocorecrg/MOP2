#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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

granularity				  : ${params.annotation}

ref_type                  : ${params.ref_type}
seq_type                  : ${params.seq_type}

output                    : ${params.output}
qualityqc                 : ${params.qualityqc}

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
gpu				    = params.GPU


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

if (params.ref_type == "genome") {
	if (params.annotation != "") {
		annotation = file(params.annotation)
		if( !annotation.exists() ) exit 1, "Missing annotation file: ${params.annotation}!"
	}
}
 
def subworkflowsDir = "${baseDir}/BioNextflow/subworkflows"

def guppy_basecall_label = (params.GPU == 'ON' ? 'basecall_gpus' : 'basecall_cpus')

include { GET_WORKFLOWS; BASECALL as BASECALL_GUPPY; BASECALL_DEMULTI as GUPPY_BASECALL_DEMULTI } from "${subworkflowsDir}/basecalling/guppy.nf" addParams(EXTRAPARS_BC: params.basecaller_opt, EXTRAPARS_DEM: params.demultiplexing_opt, LABEL: guppy_basecall_label, GPU_OPTION: gpu, OUTPUT: outputFast5)




fast5_files.map { 
    def filepath = file(it)
    def file_parts = "${filepath}".tokenize("/")
    def folder_name  = filepath[-2]
 	[folder_name, it]
}.groupTuple().set{ fast5_per_folder}



workflow flow1 {		
	fast5_per_folder.map{
		def folder_name = it[0]
		def buffer_files = it[1].flatten().collate(params.granularity)
		[folder_name, buffer_files]
	}.transpose().set{ fast5_4_analysis }

	GET_WORKFLOWS(params.flowcell, params.kit).view()

	if (demultiplexer == "OFF") BASECALL_GUPPY(fast5_4_analysis, params.flowcell, params.kit)
 	if (demultiplexer == "GUPPY") GUPPY_BASECALL_DEMULTI (fast5_4_analysis, params.flowcell, params.kit)
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