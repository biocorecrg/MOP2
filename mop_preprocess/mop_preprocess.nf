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
fastq                     : ${params.fastq}

reference                 : ${params.reference}
annotation                : ${params.annotation}

granularity				  : ${params.granularity}

ref_type                  : ${params.ref_type}
pars_tools				  : ${params.pars_tools}

output                    : ${params.output}
qualityqc                 : ${params.qualityqc}

GPU                       : ${params.GPU}

basecalling               : ${params.basecalling} 
demultiplexing            : ${params.demultiplexing} 
demulti_fast5		      : ${params.demulti_fast5}

filtering                 : ${params.filtering}
mapping                   : ${params.mapping}

counting                  : ${params.counting}

saveSpace   			  : ${params.saveSpace}

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
logo = file("$baseDir/img/logo_small.png")

def gpu				    = params.GPU

def tools = [:]
tools["basecalling"] = params.basecalling
tools["demultiplexing"] = params.demultiplexing
tools["mapping"] = params.mapping
tools["filtering"] = params.filtering
tools["counting"] = params.counting
//tools["variantcall"] = params.variantcall

// Output folders
outputFastq    = "${params.output}/fastq_files"
outputFast5    = "${params.output}/fast5_files"
outputQual     = "${params.output}/QC_files"
outputMultiQC  = "${params.output}/report"
outputMapping  = "${params.output}/alignment"
//outputCRAM     = "${params.output}/cram_files"
outputCounts   = "${params.output}/counts"
//outputVars     = "${params.output}/variants"
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
* This is default value in case guppy will be used for RNA demultiplexing
*/
params.barcodekit = ""

if (params.ref_type == "genome") {
	if (params.annotation != "") {
		annotation = file(params.annotation)
		if( !annotation.exists() ) exit 1, "Missing annotation file: ${params.annotation}!"
	}
}
 
def subworkflowsDir = "${baseDir}/../BioNextflow/subworkflows"
def local_modules = "${baseDir}/../local_modules"
def guppy_basecall_label = (params.GPU == 'ON' ? 'basecall_gpus' : 'basecall_cpus')
def deeplexi_basecall_label = (params.GPU == 'ON' ? 'demulti_gpus' : 'demulti_cpus')
def output_bc = (params.demulti_fast5 == 'ON' ? '' : outputFast5)

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

if (params.saveSpace == "YES") outmode = "move"
else outmode = "copy"

include { GET_WORKFLOWS; BASECALL as GUPPY_BASECALL; BASECALL_DEMULTI as GUPPY_BASECALL_DEMULTI } from "${subworkflowsDir}/basecalling/guppy" addParams(EXTRAPARS_BC: progPars["basecalling--guppy"], EXTRAPARS_DEM: progPars["demultiplexing--guppy"], LABEL: guppy_basecall_label, GPU_OPTION: gpu, MOP: "YES", OUTPUT: output_bc, OUTPUTMODE: outmode)
include { DEMULTIPLEX as DEMULTIPLEX_DEEPLEXICON } from "${subworkflowsDir}/demultiplexing/deeplexicon" addParams(EXTRAPARS: progPars["demultiplexing--deeplexicon"], LABEL:deeplexi_basecall_label, GPU_OPTION: gpu)
include { extracting_demultiplexed_fastq; extracting_demultiplexed_fast5_deeplexicon; extracting_demultiplexed_fast5_guppy} from "${local_modules}" addParams(LABEL: 'big_cpus')
include { FILTER as NANOFILT_FILTER} from "${subworkflowsDir}/trimming/nanofilt" addParams(EXTRAPARS: progPars["filtering--nanofilt"])
include { MAP as GRAPHMAP} from "${subworkflowsDir}/alignment/graphmap" addParams(EXTRAPARS: progPars["mapping--graphmap"])
include { MAP as GRAPHMAP2} from "${subworkflowsDir}/alignment/graphmap2" addParams(EXTRAPARS: progPars["mapping--graphmap2"])
include { MAP as MINIMAP2} from "${subworkflowsDir}/alignment/minimap2" addParams(EXTRAPARS: progPars["mapping--minimap2"])
include { FASTQCP as FASTQC} from "${subworkflowsDir}/qc/fastqc" addParams(LABEL: 'big_cpus')
include { SORT as SAMTOOLS_SORT; INDEX as SAMTOOLS_INDEX } from "${subworkflowsDir}/misc/samtools" addParams(LABEL: 'big_cpus', OUTPUT:outputMapping)
include { CAT as SAMTOOLS_CAT } from "${subworkflowsDir}/misc/samtools"
include { MOP_QC as NANOPLOT_QC } from "${subworkflowsDir}/qc/nanoplot" 
include { COUNT as NANOCOUNT } from "${subworkflowsDir}/read_count/nanocount" addParams(EXTRAPARS: progPars["counting--nanocount"], OUTPUT:outputCounts)
include { COUNT_AND_ANNO as HTSEQ_COUNT } from "${subworkflowsDir}/read_count/htseq" addParams(EXTRAPARS: progPars["counting--htseq"], OUTPUT:outputCounts)
include { REPORT as MULTIQC } from "${subworkflowsDir}/reporting/multiqc" addParams(EXTRAPARS: "-c ${config_report}", OUTPUT:outputMultiQC)
include { concatenateFastQFiles; MinIONQC; bam2stats; AssignReads; countStats; joinCountStats; joinAlnStats } from "${local_modules}"
include { cleanFile as fastqCleanFile; cleanFile as bamCleanFile; cleanFile as fast5CleanFile} from "${local_modules}"


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
 	    //bc_fastq = reshapeSamples(basecalled_fastq)
		bc_fast5 = reshapeSamples(outbc.basecalled_fast5)
		bc_stats = reshapeSamples(outbc.basecalling_stats)

	emit:
    	basecalled_fast5 = bc_fast5
    	basecalled_fastq = basecalled_fastq
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
				fast5_res = extracting_demultiplexed_fast5_deeplexicon(alldemux.groupTuple().join(basecalledbc.transpose().groupTuple()))
				
				// OPTIONAL CLEANING FASTQ5 FILES
				fast5CleanFile(basecalledbc.transpose().groupTuple(), fast5_res.map{it[1]}.collect(), ".fast5")
			}
		// Demultiplex fastq	
			demufq = extracting_demultiplexed_fastq(demux.join(outbc.basecalled_fastq))

		} else if (params.demultiplexing == "guppy") {
			// IF DEMULTIPLEXING IS GUPPY	
			outbc = GUPPY_BASECALL_DEMULTI (fast5_4_analysis, params.flowcell, params.kit)
			demufq = outbc.basecalled_fastq
			fast5_res = outbc.basecalled_fast5

			// Optional demultiplex fast5 		
			if (params.demulti_fast5 == "ON" ) {
				basecalledbc = reshapeSamples(outbc.basecalled_fast5)
				alldemux = reshapeSamples(outbc.basecalling_stats)								
				fast5_res = extracting_demultiplexed_fast5_guppy(alldemux.groupTuple().join(basecalledbc.transpose().groupTuple()))
				
				// OPTIONAL CLEANING FASTQ5 FILES
				fast5CleanFile(basecalledbc.transpose().groupTuple(), fast5_res.map{it[1]}.collect(), ".fast5")
			}


		}
		reshapedDemufq = demufq.transpose().map{
			[it[1].name.replace(".fastq.gz", ""), it[1] ]
		}		
 		// Optional fastq filtering	
		if (params.filtering == "nanofilt") {
 			nanofilt = NANOFILT_FILTER(reshapedDemufq)
 			reshapedDemufq = nanofilt
		}

	emit:
    	basecalled_fast5 =  fast5_res
    	//basecalled_fastq = basecalled_fastq_res
    	basecalled_fastq = reshapedDemufq
    	basecalled_stats = reshapeSamples(outbc.basecalling_stats)
	
}


workflow preprocess_flow {
    take:
    	bc_fast5
    	bc_fastq
    	basecalled_stats
    	
	main:	
	// Perform MinIONQC on basecalling stats
	basecall_qc = MinIONQC(basecalled_stats.groupTuple())

	// Perform mapping on fastq files
	if (params.mapping == "NO") {
		stats_aln = Channel.value()	
		sorted_alns = Channel.value()	
		nanoplot_qcs = Channel.value()	
	}
	else {
		switch(params.mapping) { 
   			case "graphmap": 
   			aln_reads = GRAPHMAP(bc_fastq, reference)
   			break
   			case "graphmap2": 
   			aln_reads = GRAPHMAP2(bc_fastq, reference)
   			case "minimap2": 
   			break
   			aln_reads = MINIMAP2(bc_fastq, reference)
   			break
   			default: 
			println "ERROR ################################################################"
			println "${params.mapping} is not a supported alignment"
			println "ERROR ################################################################"
			println "Exiting ..."
			System.exit(0)
			break

		}	 

		// Concatenate bamfiles
 	    if (params.demultiplexing == "NO" ) reshaped_aln_reads = reshapeSamples(aln_reads)
		else reshaped_aln_reads = reshapeDemuxSamples(aln_reads)

		jaln_reads = SAMTOOLS_CAT(reshaped_aln_reads.groupTuple())

		// Perform SORTING and INDEXING on bam files
		sorted_alns = SAMTOOLS_SORT(jaln_reads)
		aln_indexes = SAMTOOLS_INDEX(sorted_alns)

		// OPTIONAL CLEANING BAM FILES
		bamCleanFile(reshaped_aln_reads.groupTuple(), jaln_reads.map{it[1]}.collect(), ".bam")

		// Perform bam2stats on sorted bams
		aln_stats = bam2stats(sorted_alns)
		stats_aln = joinAlnStats(aln_stats.map{ it[1]}.collect())
	
		// Perform NanoPlot on sorted bams
		nanoplot_qcs = NANOPLOT_QC(sorted_alns)

	}	

	// Concatenate fastq files
	if (params.demultiplexing == "NO" ) reshaped_bc_fastq = reshapeSamples(bc_fastq)
	else reshaped_bc_fastq = reshapeDemuxSamples(bc_fastq)

	fastq_files = concatenateFastQFiles(reshaped_bc_fastq.groupTuple())

	// Perform fastqc QC on fastq
	fastqc_files = FASTQC(fastq_files)

	// OPTIONAL CLEANING FASTQC FILES
	fastqCleanFile(reshaped_bc_fastq.groupTuple(), fastq_files.map{it[1]}.collect().mix(fastqc_files.map{it[1]}.collect(), jaln_reads.map{it[1]}.collect()).collect(), ".gz")
	
	// OPTIONAL Perform COUNTING / ASSIGNMENT
	if (params.counting == "nanocount" && params.ref_type == "transcriptome") {
		read_counts = NANOCOUNT(sorted_alns)
		assignments = AssignReads(sorted_alns, "nanocount")
		stat_counts = countStats(assignments)
		stats_counts = joinCountStats(stat_counts.map{ it[1]}.collect())
	}
	else if (params.counting == "htseq" && params.ref_type == "genome") {
		htseq_out = HTSEQ_COUNT(params.annotation, sorted_alns)
		read_counts = htseq_out.counts
		assignments = AssignReads(htseq_out.bam, "htseq")
		stat_counts = countStats(assignments)
		stats_counts = joinCountStats(stat_counts.map{ it[1]}.collect())
	} else if (params.counting == "NO") {
		// Default empty channels for reporting
		stats_counts = Channel.value()
	} else {
		println "ERROR ################################################################"
		println "${params.counting} is not compatible with ${params.ref_type}"
		println "htseq requires a genome as reference and an annotation in GTF"
		println "nanocount requires a transcriptome as a reference"		
		println "ERROR ################################################################"
		println "Exiting ..."
		System.exit(0)
	} 
	
	fastqc_files.mix(basecall_qc).map{it[1]}.set{qcs}

	// Perform MULTIQC report
	MULTIQC(qcs.mix(stats_counts, stats_aln, Channel.from(logo)).collect())
	
}


workflow preprocess_simple {
    take:
    	bc_fastq
    	
	main:	
	// Perform mapping on fastq files
	if (params.mapping == "NO") {
		stats_aln = Channel.value()	
		sorted_alns = Channel.value()	
		nanoplot_qcs = Channel.value()	
	}
	else {
		switch(params.mapping) { 
   			case "graphmap": 
   			aln_reads = GRAPHMAP(bc_fastq, reference)
   			break
   			case "graphmap2": 
   			aln_reads = GRAPHMAP2(bc_fastq, reference)
   			case "minimap2": 
   			break
   			aln_reads = MINIMAP2(bc_fastq, reference)
   			break
   			default: 
			println "ERROR ################################################################"
			println "${params.mapping} is not a supported alignment"
			println "ERROR ################################################################"
			println "Exiting ..."
			System.exit(0)
			break
		}	 

		// Perform SORTING and INDEXING on bam files
		sorted_alns = SAMTOOLS_SORT(aln_reads)
		aln_indexes = SAMTOOLS_INDEX(sorted_alns)

		// Perform bam2stats on sorted bams
		aln_stats = bam2stats(sorted_alns)
		stats_aln = joinAlnStats(aln_stats.map{ it[1]}.collect())
	
		// Perform NanoPlot on sorted bams
		nanoplot_qcs = NANOPLOT_QC(sorted_alns)
	}	

	// Perform fastqc QC on fastq
	fastqc_files = FASTQC(bc_fastq)
	
	// OPTIONAL Perform COUNTING / ASSIGNMENT
	if (params.counting == "nanocount" && params.ref_type == "transcriptome") {
		read_counts = NANOCOUNT(sorted_alns)
		assignments = AssignReads(sorted_alns, "nanocount")
		stat_counts = countStats(assignments)
		stats_counts = joinCountStats(stat_counts.map{ it[1]}.collect())
	}
	else if (params.counting == "htseq" && params.ref_type == "genome") {
		htseq_out = HTSEQ_COUNT(params.annotation, sorted_alns)
		read_counts = htseq_out.counts
		assignments = AssignReads(htseq_out.bam, "htseq")
		stat_counts = countStats(assignments)
		stats_counts = joinCountStats(stat_counts.map{ it[1]}.collect())
	} else if (params.counting == "NO") {
		// Default empty channels for reporting
		stats_counts = Channel.value()
	} else {
		println "ERROR ################################################################"
		println "${params.counting} is not compatible with ${params.ref_type}"
		println "htseq requires a genome as reference and an annotation in GTF"
		println "nanocount requires a transcriptome as a reference"		
		println "ERROR ################################################################"
		println "Exiting ..."
		System.exit(0)
	} 
	
	//fastqc_files.mix(map{it[1]}).set{qcs}

	// Perform MULTIQC report
	//MULTIQC(qcs.mix(stats_counts, stats_aln, Channel.from(logo)).collect())
	
}


 workflow {
 	if (params.fast5 != "" && params.fastq == "") {

		Channel
			.fromPath( params.fast5)                                             
			.ifEmpty { error "Cannot find any file matching: ${params.fast5}" }
			.set {fast5_files}

		fast5_files.map { 
			def filepath = file(it)
			def file_parts = "${filepath}".tokenize("/")
			def folder_name  = filepath[-2]
			[folder_name, it]
		}.groupTuple().set{ fast5_per_folder}
	
		// Check tools
		checkTools(tools, progPars)

 		def num = 0
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
		def basecalled_stats = outf.basecalled_stats

		preprocess_flow(bc_fast5, bc_fastq, basecalled_stats)
		
 	} else if(params.fast5 == "" && params.fastq != "") {

 		// Check tools
		tools["basecalling"] = "NO"
		tools["demultiplexing"] = "NO"
		checkTools(tools, progPars)
		Channel.fromFilePairs( params.fastq , size: 1)
			.ifEmpty { error "Cannot find any file matching: ${params.fastq}" }
			.set {fastq_files}
		
		preprocess_simple(fastq_files)
	
		
 	} else {
			println "ERROR ################################################################"
			println "Please choose one between fast5 and fastq as input!!!" 
			println "ERROR ################################################################"
			println "Exiting ..."
			System.exit(0)
 		
 	}
 	
 
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

def reshapeDemuxSamples(inputChannel) {
	def reshapedChannel = inputChannel.map {
 		def ids = it[0].split("---")
 		def dems = it[0].split("\\.")
 			["${ids[0]}---${dems[-1]}", it[1]]
	 }
	 return(reshapedChannel)
}

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