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

// include functions, outdirs from other files
evaluate(new File("../outdirs.nf"))
def local_modules = file("$baseDir/../local_modules.nf")

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

include { getParameters; mapIDPairs } from "${local_modules}" 

// Create a channel for tool options
progPars = getParameters(params.pars_tools)

include { CALC_VAR_FREQUENCIES as EPINANO_CALC_VAR_FREQUENCIES } from "${subworkflowsDir}/chem_modification/epinano" addParams(LABEL: 'big_cpus', OUTPUT: outputEpinanoFlow, EXTRAPARS: progPars["epinano--epinano"])
include { EVENTALIGN as NANOPOLISH_EVENTALIGN } from "${subworkflowsDir}/chem_modification/nanopolish" addParams(LABEL: 'big_cpus',  OUTPUT: outputNanoPolComFlow, EXTRAPARS: progPars["nanocompore--nanopolish"])
include { SAMPLE_COMPARE as NANOCOMPORE_SAMPLE_COMPARE } from "${subworkflowsDir}/chem_modification/nanocompore" addParams(LABEL: 'big_cpus',  OUTPUT: outputNanoPolComFlow, EXTRAPARS: progPars["nanocompore--nanocompore"])
include { RESQUIGGLE_RNA as TOMBO_RESQUIGGLE_RNA } from "${subworkflowsDir}/chem_modification/tombo.nf" addParams(LABEL: 'big_cpus', EXTRAPARS: progPars["tombo_resquiggling--tombo"])
include { GET_MODIFICATION_MSC as TOMBO_GET_MODIFICATION_MSC } from "${subworkflowsDir}/chem_modification/tombo.nf" addParams(LABEL: 'big_cpus', EXTRAPARS: progPars["tombo_msc--tombo"])
include { GET_MODIFICATION_LSC as TOMBO_GET_MODIFICATION_LSC } from "${subworkflowsDir}/chem_modification/tombo.nf" addParams(LABEL: 'big_cpus', EXTRAPARS: progPars["tombo_lsc--tombo"])

include { GET_VERSION as EPINANO_VER } from "${subworkflowsDir}/chem_modification/epinano" 
include { GET_VERSION as NANOPOLISH_VER } from "${subworkflowsDir}/chem_modification/nanopolish" 
include { GET_VERSION as NANOCOMPORE_VER } from "${subworkflowsDir}/chem_modification/nanocompore" 
include { GET_VERSION as TOMBO_VER } from "${subworkflowsDir}/chem_modification/tombo.nf"

include { indexReference; callVariants; checkRef; bedGraphToWig as bedGraphToWig_msc; bedGraphToWig as bedGraphToWig_lsc } from "${local_modules}"
include {  mergeTomboWigs as mergeTomboWigsPlus; mergeTomboWigs as mergeTomboWigsMinus} addParams(OUTPUT: outputTomboFlow) from "${local_modules}"
include { makeEpinanoPlots as makeEpinanoPlots_mis; makeEpinanoPlots as makeEpinanoPlots_ins; makeEpinanoPlots as makeEpinanoPlots_del } addParams(OUTPUT: outputEpinanoFlow) from "${local_modules}"

include { multiToSingleFast5 } addParams(LABEL: 'big_cpus') from "${local_modules}"



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

	ref_file = checkRef(reference)

	if (params.epinano == "YES") {
		epinano_flow(bams, ref_file, comparisons)
	}
	if (params.nanocompore == "YES") {
		compore_polish_flow(comparisons, fast5_folders, bams, bais, fastqs, summaries, ref_file) 
	}

	if (params.tombo_lsc == "YES" || params.tombo_msc == "YES") {
		tombo_data = tombo_common_flow(fast5_files, ref_file, comparisons)
	    
		if (params.tombo_msc == "YES") {
			tombo_msc_flow(tombo_data, ref_file)
			wiggle_msc = bedGraphToWig_msc(tombo_msc_flow.out.bed_graphs.transpose()).map{
				["${it[0]}_msc", it[1] ]
			}
			stat_msc = tombo_msc_flow.out.dampened_wiggles.transpose().map{
				["${it[0]}_msc", it[1] ]
			}				
		}
		if (params.tombo_lsc == "YES") {
			tombo_lsc_flow(tombo_data, ref_file)
			wiggle_lsc = bedGraphToWig_lsc(tombo_lsc_flow.out.bed_graphs.transpose()).map{
				["${it[0]}_lsc", it[1] ]
			}
			stat_lsc = tombo_lsc_flow.out.dampened_wiggles.transpose().map{
				["${it[0]}_lsc", it[1] ]
			}
		}

			wiggle_msc.mix(wiggle_lsc).branch {
        		sampleplus: it[1] =~ /\.sample\.plus\./
        		sampleminus: it[1] =~ /\.sample\.minus\./
        		controlplus: it[1] =~ /\.control\.plus\./
        		controlminus: it[1] =~ /\.control\.minus\./
    		}.set{combo_tombo}
			stat_lsc.mix(stat_msc).branch {
        		plus: it[1] =~ /\.plus\./
        		minus: it[1] =~ /\.minus\./
    		}.set{combo_stats}

    		mergeTomboWigsPlus("plus", mergeTomboScript, combo_tombo.sampleplus.join(combo_tombo.controlplus).join(combo_stats.plus))
    		mergeTomboWigsMinus("minus", mergeTomboScript, combo_tombo.sampleminus.join(combo_tombo.controlminus).join(combo_stats.minus))
	}
	
	all_ver = EPINANO_VER().mix(NANOPOLISH_VER())
	.mix(NANOCOMPORE_VER()).mix(TOMBO_VER())
	.collectFile(name: 'tool_version.txt', newLine: false, storeDir:params.output)

}

workflow tombo_common_flow {
    take:
	fast5_files
	ref_file
	comparisons
	
	main:
	fast5_files.map{
		["${it[0]}---${it[1].simpleName}", it[1]]
	}.set{fast5_reshaped}
	
	single_fast5_folders = multiToSingleFast5(fast5_reshaped)
	resquiggle = TOMBO_RESQUIGGLE_RNA(single_fast5_folders, ref_file)
	
	resquiggle.join(single_fast5_folders).map{
		def ids = it[0].split("---")
		["${ids[0]}", it[1], it[2]]
	}.groupTuple().map{
		[it[0], [it[1], it[2]]]
	}.set{reshape_resquiggle}
	
	data_for_tombo = mapIDPairs(comparisons, reshape_resquiggle).map{
		[it[0], it[1], it[2][0], it[2][1], it[3][0], it[3][1]]
	}
	
	emit:
		data_for_tombo
}

workflow tombo_msc_flow {
    take:
	data_for_tombo
	reference
	
	main:
	TOMBO_GET_MODIFICATION_MSC(data_for_tombo, reference)
	bed_graphs = TOMBO_GET_MODIFICATION_MSC.out.bedgraphs
	dampened_wiggles = TOMBO_GET_MODIFICATION_MSC.out.dampened_wiggles
	
	emit:
		bed_graphs
		dampened_wiggles
	
}

workflow tombo_lsc_flow {
    take:
	data_for_tombo
	reference
	
	main:
	TOMBO_GET_MODIFICATION_LSC(data_for_tombo, reference)
	bed_graphs = TOMBO_GET_MODIFICATION_LSC.out.bedgraphs
	dampened_wiggles = TOMBO_GET_MODIFICATION_LSC.out.dampened_wiggles
	
	emit:
		bed_graphs
		dampened_wiggles

}

workflow compore_polish_flow {
    take:
		comparisons
		fast5_folders
		bams
		bais
		fastqs
		summaries
		ref_file		
	
	main:	
		concat_events = NANOPOLISH_EVENTALIGN(fast5_folders, bams, bais, fastqs, summaries, ref_file)
		combs_events = mapIDPairs(comparisons, concat_events)
		NANOCOMPORE_SAMPLE_COMPARE(combs_events, ref_file)
	
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
