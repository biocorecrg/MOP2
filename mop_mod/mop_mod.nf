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

╔╦╗╔═╗╔═╗  ╔╦╗┌─┐┌┬┐
║║║║ ║╠═╝  ║║║│ │ ││
╩ ╩╚═╝╩    ╩ ╩└─┘─┴┘
                                                                                       
====================================================
BIOCORE@CRG Master of Pores 2. Detection of RNA modification - N F  ~  version ${version}
====================================================

*****************   Input files    *******************
input_path                              : ${params.input_path}
comparison                              : ${params.comparison}

********** reference has to be the genome *************
reference                               : ${params.reference}
output                                  : ${params.output}

coverage                                : ${params.coverage}
pars_tools								: ${params.pars_tools}

************************* Flows *******************************
epinano                             	: ${params.epinano}
nanocompore                             : ${params.nanocompore}
tombo_lsc                               : ${params.tombo_lsc}
tombo_msc                               : ${params.tombo_msc}


email                                   : ${params.email}
"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// check input files
reference = file(params.reference)
if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"

def local_modules = file("$baseDir/local_modules.nf")
def subworkflowsDir = "${baseDir}/../BioNextflow/subworkflows"

model_folder = file("$baseDir/models/")
if( !model_folder.exists() ) exit 1, "Missing folders with EpiNano's models!"
joinScript = file("$baseDir/bin/join.r")
mergeTomboScript = file("$baseDir/bin/Merge_Tombo_wigs_MoP_v2.R")

def flows = [:]
flows["nanomod"] = params.epinano
flows["nanocompore"] = params.nanocompore
flows["tombo_lsc"] = params.tombo_msc
flows["tombo_msc"] = params.tombo_lsc

// Output folders
outputEpinanoFlow    = "${params.output}/epinano_flow"


include { indexReference; callVariants; checkRef } from "${local_modules}"
include { CALC_VAR_FREQUENCIES as EPINANO_CALC_VAR_FREQUENCIES } from "${subworkflowsDir}/chem_modification/epinano" addParams(LABEL: 'big_cpus', OUTPUT: outputEpinanoFlow)
include { EVENTALIGN as NANOPOLISH_EVENTALIGN } from "${subworkflowsDir}/chem_modification/nanopolish" addParams(LABEL: 'big_cpus')

include { GET_VERSION } from "${subworkflowsDir}/chem_modification/epinano" addParams(LABEL: 'big_cpus', OUTPUT: outputEpinanoFlow)
include { makeEpinanoPlots as makeEpinanoPlots_mis; makeEpinanoPlots as makeEpinanoPlots_ins; makeEpinanoPlots as makeEpinanoPlots_del } from "${local_modules}"


// Create a channel for tool options
pars_tools = file(params.pars_tools)
if( !pars_tools.exists() ) exit 1, "Missing tools options config: '$pars_tools'"

def progPars = [:]
def tooList = [:]
allLines  = pars_tools.readLines()

for( line : allLines ) {
    list = line.split("\t")
    if (list.length <3) {
		 error "ERROR!!! Tool option file has to be tab separated\n" 
	}
    if (!(list[0] =~ /#/ )) {
		progPars["${list[0]}--${list[1]}"] = list[2].replace("\"", "").replace('$baseDir', "${baseDir}").replace('${baseDir}', "${baseDir}")
    }  
}

/*
 * Creates the channels with comparisons
 */

compfile = file(params.comparison)
if( !compfile.exists() ) exit 1, "Missing comparison file: ${compfile}. Specify path with --comparisons"

 Channel
    .from(compfile.readLines())
    .map { line ->
        list = line.split("\t")
        if (list[0]!= "") {
            def sampleID = list[0]
            def ctrlID = list[1]
            [ sampleID, ctrlID ]
        }
    }.set {comparisons}

// Output folders
output_epinano_flow = "${params.output}/Epinano_flow"
output_tombo_lsc_flow = "${params.output}/Tombo_lsc_flow"
output_tombo_msc_flow = "${params.output}/Tombo_msc_flow"
output_nanocompore_flow = "${params.output}/Nanocompore_flow"


workflow {	
	comparisons.flatten().unique().set{unique_samples}
	
	unique_samples.map {
 	  	 [it, file("${params.input_path}/alignment/${it}_s.bam")]
	}.transpose().set{bams}
	unique_samples.map {
 	  	 [it, file("${params.input_path}/alignment/${it}_s.bam.bai")]
	}.transpose().set{bais}
	unique_samples.map {
 	  	 [it, file("${params.input_path}/fastq_files/${it}.fq.gz")]
	}.transpose().set{fastqs}
	unique_samples.map {
 	  	 [it, file("${params.input_path}/QC_files/${it}_final_summary.stats")]
	}.transpose().set{summaries}
	unique_samples.map {
 	  	 [it, file("${params.input_path}/fast5_files/${it}/")]
	}.transpose().set{fast5_folders}

	unique_samples.map {
    [it, file("${params.input_path}/fast5_files/${it}/*.fast5")]
	}.transpose().set{fast5_files}
//fast5_files.view()

	ref_file = checkRef(reference)

	if (params.epinano == "YES") {
		epinano_flow(bams, ref_file, comparisons)
	}
	if (params.nanocompore == "YES") {
		event_align = NANOPOLISH_EVENTALIGN(fast5_folders, bams, bais, fastqs, summaries, ref_file)
		//event_align.groupTuple() 
		//nanocompore_flow(bams, reference)
	}
	
	
	
}

workflow epinano_flow {
    take:
		bams
		reference
		comparisons
		
	main:	
	indexes = indexReference(reference)
	variants = callVariants(bams, reference.combine(indexes))
	per_site_vars = EPINANO_CALC_VAR_FREQUENCIES(variants).per_site_vars
	per_site_vars.combine(per_site_vars).map {
		[ it[0], it[2], it[1], it[3] ]
	}.join(comparisons, by:[0,1]).set{per_site_for_plots}

	makeEpinanoPlots_ins(per_site_for_plots, "ins")
	makeEpinanoPlots_mis(per_site_for_plots, "mis")
	makeEpinanoPlots_del(per_site_for_plots, "del")
}




/*
*  Finish message
*/
workflow.onComplete {
    println "Pipeline BIOCORE@CRG Master of Pore completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

/*
* Mail notification
*/

if (params.email == "yourmail@yourdomain" || params.email == "") { 
    log.info 'Skipping the email\n'
}
else {
    log.info "Sending the email to ${params.email}\n"

    workflow.onComplete {

    def msg = """\
        Pipeline BIOCORE@CRG Master of Pore 2 modification module's execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        Error report: ${workflow.errorReport ?: '-'}
        """
        .stripIndent()

        sendMail(to: params.email, subject: "Master of Pore 2 execution", body: msg)
    }
}
