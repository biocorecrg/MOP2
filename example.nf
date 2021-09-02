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

GPU                       : ${params.GPU}

basecalling               : ${params.basecalling} 
demultiplexing            : ${params.demultiplexing} 
demulti_fast5		      : ${params.demulti_fast5}

filtering                 : ${params.filtering}
mapping                   : ${params.mapping}

map_type                  : ${params.map_type}

counting                  : ${params.counting}

downsampling			  : ${params.downsampling}

variantcall          	  : ${params.variantcall}

email                     : ${params.email}
"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// check multi5 and GPU usage. GPU maybe can be removed as param if there is a way to detect it
if (params.GPU != "ON" && params.GPU != "OFF") exit 1, "Please specify ON or OFF in GPU processors are available"

// check input files
reference = file(params.reference)
if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"
config_report = file("$baseDir/config.yaml")
if( !config_report.exists() ) exit 1, "Missing config.yaml file!"
logo = file("$baseDir/../docs/logo_small.png")

def gpu				    = params.GPU

def tools = [:]
tools["basecalling"] = params.basecalling
tools["demultiplexing"] = params.demultiplexing
tools["mapping"] = params.mapping
tools["filtering"] = params.filtering
tools["counting"] = params.counting
tools["variantcall"] = params.variantcall

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

if (params.ref_type == "genome") {
	if (params.annotation != "") {
		annotation = file(params.annotation)
		if( !annotation.exists() ) exit 1, "Missing annotation file: ${params.annotation}!"
	}
}
 
def subworkflowsDir = "${baseDir}/BioNextflow/subworkflows"
def guppy_basecall_label = (params.GPU == 'ON' ? 'basecall_gpus' : 'basecall_cpus')
def deeplexi_basecall_label = (params.GPU == 'ON' ? 'demulti_gpus' : 'demulti_cpus')

// Create a channel for tool options
pars_tools = file("${baseDir}/tool_opt.tsv")
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

// Check tools
checkTools(tools, progPars)

include { GET_WORKFLOWS; BASECALL as GUPPY_BASECALL; BASECALL_DEMULTI as GUPPY_BASECALL_DEMULTI } from "${subworkflowsDir}/basecalling/guppy" addParams(EXTRAPARS_BC: progPars["basecalling--guppy"], EXTRAPARS_DEM: progPars["demultiplexing--guppy"], LABEL: guppy_basecall_label, GPU_OPTION: gpu)
include { DEMULTIPLEX as DEMULTIPLEX_DEEPLEXICON } from "${subworkflowsDir}/demultiplexing/deeplexicon" addParams(EXTRAPARS: progPars["demultiplexing--deeplexicon"], LABEL:deeplexi_basecall_label, GPU_OPTION: gpu)
include { extracting_demultiplexed_fastq; extracting_demultiplexed_fast5} from "${baseDir}/local_modules"
include { FILTER as NANOFILT_FILTER} from "${subworkflowsDir}/trimming/nanofilt" addParams(EXTRAPARS: progPars["filtering--nanofilt"])
include { MAP as GRAPHMAP} from "${subworkflowsDir}/alignment/graphmap" addParams(EXTRAPARS: progPars["mapping--graphmap"])
include { MAP as GRAPHMAP2} from "${subworkflowsDir}/alignment/graphmap2" addParams(EXTRAPARS: progPars["mapping--graphmap2"])
include { MAP as MINIMAP2} from "${subworkflowsDir}/alignment/minimap2" addParams(EXTRAPARS: progPars["mapping--minimap2"])
include { FASTQCP as FASTQC} from "${subworkflowsDir}/qc/fastqc" addParams(LABEL: 'big_cpus')
include { SORT as SAMTOOLS_SORT; INDEX as SAMTOOLS_INDEX } from "${subworkflowsDir}/misc/samtools" addParams(LABEL: 'big_cpus')
include { MOP_QC as NANOPLOT_QC } from "${subworkflowsDir}/qc/nanoplot" 
include { COUNT as NANOCOUNT } from "${subworkflowsDir}/read_count/nanocount" addParams(EXTRAPARS: progPars["counting--nanocount"])
include { COUNT_AND_ANNO as HTSEQ_COUNT } from "${subworkflowsDir}/read_count/htseq" addParams(EXTRAPARS: progPars["counting--htseq"])
include { concatenateFastQFiles } from "${baseDir}/local_modules"
include { MinIONQC } from "${baseDir}/local_modules"
include { bam2stats } from "${baseDir}/local_modules"
include { AssignReads } from "${baseDir}/local_modules"
include { countStats } from "${baseDir}/local_modules"
include { joinCountsStats } from "${baseDir}/local_modules"

fast5_files.map { 
    def filepath = file(it)
    def file_parts = "${filepath}".tokenize("/")
    def folder_name  = filepath[-2]
 	[folder_name, it]
}.groupTuple().set{ fast5_per_folder}
	

/*
* Simple flow of basecalling
*/
workflow flow1 {
    take: 
    	fast5_4_analysis
    main:
		outbc = GUPPY_BASECALL (fast5_4_analysis, params.flowcell, params.kit)
	    basecalled_fastq = outbc.basecalled_fastq

		// Optional fastq filtering	
		if (params.filtering == "nanofilt") {
			nanofilt = NANOFILT_FILTER(outbc.basecalled_fastq)
			basecalled_fastq = reshapeSamples(nanofilt.out)
		} 
 	    bc_fastq = reshapeSamples(basecalled_fastq)
		bc_fast5 = reshapeSamples(outbc.basecalled_fast5)
		bc_stats = reshapeSamples(outbc.basecalling_stats)

	emit:
    	basecalled_fast5 = bc_fast5
    	basecalled_fastq = bc_fastq
    	basecalled_stats = bc_stats

}

/*
*  Basecalling and Demultiplexing
*/

workflow flow2 {		
    take: 
    	fast5_4_analysis
    main:
		// IF DEMULTIPLEXING IS DEEPLEXICON	
    	if(params.demultiplexing == "deeplexicon") {
			outbc = GUPPY_BASECALL(fast5_4_analysis, params.flowcell, params.kit)
			demux = DEMULTIPLEX_DEEPLEXICON(fast5_4_analysis)
			fast5_res = outbc.basecalled_fast5
		
			// Optional demultiplex fast5 		
			if (params.demulti_fast5 == "ON" ) {
				basecalledbc = reshapeSamples(outbc.basecalled_fast5)
				alldemux = reshapeSamples(demux)
				fast5_res = extracting_demultiplexed_fast5(alldemux.groupTuple().join(basecalledbc.groupTuple()))
			}
		// Demultiplex fastq	
			demufq = extracting_demultiplexed_fastq(demux.join(outbc.basecalled_fastq))

		} else if (params.demultiplexing == "guppy") {
			// IF DEMULTIPLEXING IS GUPPY	
			outbc = GUPPY_BASECALL_DEMULTI (fast5_4_analysis, params.flowcell, params.kit)
			demufq = outbc.basecalled_fastq
			fast5_res = outbc.basecalled_fast5
		}
		reshapedDemufq = demufq.transpose().map{
			[it[1].name.replace(".fastq.gz", ""), it[1] ]
		}		
 		// Optional fastq filtering	
		if (params.filtering == "nanofilt") {
 			nanofilt = NANOFILT_FILTER(reshapedDemufq)
 			reshapedDemufq = nanofilt
		}

		reshapedDemufq.map {
 			def nano_ids = it[0].split("---")
 			def nano_dems = it[0].split("\\.")
 				["${nano_ids[0]}---${nano_dems[-1]}", it[1]]
	 	}.set{basecalled_fastq_res}

	emit:
    	basecalled_fast5 =  fast5_res
    	basecalled_fastq = basecalled_fastq_res
    	basecalled_stats = reshapeSamples(outbc.basecalling_stats)

		
}


def num = 0
workflow {
	fast5_per_folder.map{
		def folder_name = it[0]
		def buffer_files = it[1].flatten().collate(params.granularity)
		[folder_name, buffer_files]
	}.transpose().map{
	    num++ 
		[ "${it[0]}---${num}", it[1] ]
	}.set{ fast5_4_analysis }

	GET_WORKFLOWS(params.flowcell, params.kit).view()
	if (params.basecalling == "guppy" && params.demultiplexing == "NO" ) outf = flow1(fast5_4_analysis)
	else outf = flow2(fast5_4_analysis)

	def bc_fast5 = outf.basecalled_fast5
	def bc_fastq = outf.basecalled_fastq
	
	// Concatenate fastq files
	fastq_files = concatenateFastQFiles(bc_fastq.groupTuple())

	// Perform fastqc QC on fastq
	fastqc_files = FASTQC(fastq_files)

	// Perform MinIONQC on basecalling stats
	MinIONQC(outf.basecalled_stats.groupTuple())

	// Perform mapping on fastq files
	if (params.mapping == "graphmap") aln_reads = GRAPHMAP(fastq_files, reference)
	if (params.mapping == "graphmap2") aln_reads = GRAPHMAP2(fastq_files, reference)
	if (params.mapping == "minimap2") aln_reads = MINIMAP2(fastq_files, reference)
	
	// Perform SORTING and INDEXING on bam files
	sorted_alns = SAMTOOLS_SORT(aln_reads)
	aln_indexes = SAMTOOLS_INDEX(sorted_alns)

	// Perform bam2stats on sorted bams
	aln_stats = bam2stats(sorted_alns)
	
	// ADDING NanoPlot on sorted bams
	nanoplot_qcs = NANOPLOT_QC(sorted_alns)
	
	// ADDING COUNTING / ASSIGNMENT
	if (params.counting == "nanocount" && params.ref_type == "transcriptome") {
		read_counts = NANOCOUNT(sorted_alns)
		assignments = AssignReads(sorted_alns, "nanocount")
		stat_counts = countStats(assignments)
		stats_counts = joinCountsStats(stat_counts.map{ it[1]}.collect())
	}
	else if (params.counting == "htseq" && params.ref_type == "genome") {
		htseq_out = HTSEQ_COUNT(params.annotation, sorted_alns)
		read_counts = htseq_out.counts
		assignments = AssignReads(htseq_out.bam, "htseq")
		stat_counts = countStats(assignments)
		stats_counts = joinCountsStats(stat_counts.map{ it[1]}.collect())

	} else if (params.counting == "NO") {
	} else {
		println "ERROR ################################################################"
		println "${params.counting} is not compatible with ${params.ref_type}"
		println "htseq requires a genome as reference and an annotation in GTF"
		println "nanocount requires a transcriptome as a reference"		
		println "ERROR ################################################################"
		println "Exiting ..."
		System.exit(0)
	} 
	
	stats_counts.view()
	// ADDING COUNTING QC
	
	// ADDING REPORTING 

	// ADDING MULTIQC
	
	
}


workflow.onComplete {
    println "Pipeline BIOCORE@CRG Master of Pore completed!"
    println "Started at  $workflow.start" 
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

 


/*
* FUNCTIONS
*/

def reshapeSamples(inputChannel) {
    def reshapedChannel = inputChannel.map{
		def ids = it[0].split("---")
		["${ids[0]}", it[1]]
	}
	return(reshapedChannel)
}

def checkTools(tool_names, tool_lists) {
	println "----------------------CHECK TOOLS -----------------------------"
	tool_names.each{ key, value -> 
		if (value == "NO" ) {
			println "> ${key} will be skipped"
		} else {
			def combid = "${key}--${value}".toString()
			if (tool_lists.containsKey(combid)) {
				println "${key} : ${value}"	
			} else {
				println "ERROR ################################################################"
				println "${value} is not a valid program for ${key}"
				println "ERROR ################################################################"
				println "Exiting ..."
				System.exit(0)
			}
		}
	}
	println "--------------------------------------------------------------"
}