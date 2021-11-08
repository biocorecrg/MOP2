// LOCAL MODULES FOR MOP

params.LABEL = ""
params.OUTPUT = ""
params.saveSpace = "NO"

// MODULES 
// MOP_PREPROCESS
process extracting_demultiplexed_fastq {
    label (params.LABEL)
    tag "${ idfile }"
			
	input:
	tuple val(idfile), path(demux), path(fastq) 
	
	
	output:
	tuple val(idfile), path ("*.fastq.gz")

	script:
	"""
		extract_sequence_from_fastq.py ${demux} ${fastq}
		for i in *.fastq; do gzip \$i; done
	"""
}

process extracting_demultiplexed_fast5_deeplexicon {
    label (params.LABEL)
	container 'lpryszcz/deeplexicon:1.2.0'
    tag "${ idfile }"
    if (params.saveSpace == "YES") publishDir(params.OUTPUT, mode:'move') 
    else publishDir(params.OUTPUT, mode:'copy')    
		
	input:
	tuple val(idfile), path("demux_*"), file("*")

	output:
	file("${idfile}-*")
	
	script:
	"""
	cat demux_* | grep -v ReadID >> dem.files
	awk '{print \$2 > \$3".list" }' dem.files
	for i in *.list; do mkdir ${idfile}---`basename \$i .list`; fast5_subset --input ./ --save_path ${idfile}---`basename \$i .list`/ --read_id_list \$i --batch_size 4000 -c vbz -t ${task.cpus}; done 
	rm *.list
	rm */filename_mapping.txt
	rm dem.files 
	"""
} 

process extracting_demultiplexed_fast5_guppy {
    tag "${ idfile }"
    label (params.LABEL)
    if (params.saveSpace == "YES") publishDir(params.OUTPUT, mode:'move') 
    else publishDir(params.OUTPUT, mode:'copy')    

    container "quay.io/biocontainers/ont-fast5-api:4.0.0--pyhdfd78af_0"

             
	input:
	tuple val(idfile), path("summaries_*"), file("*")
    
	output:
	file("${idfile}-*")

    script:
    """
      if [ -f "summaries_" ]; then
	  ln -s summaries_ final_summary.stats
	  else 
		  head -n 1 summaries_1 > final_summary.stats
	      for i in summaries_*; do grep -v "filename" \$i | awk -F"\t" -v id=${idfile}  '{OFS="\t"; \$19 = id"---"\$19; print \$0}'  >> final_summary.stats; done
	  fi

		demux_fast5 -c vbz -t ${task.cpus} --input ./ --save_path ./ --summary_file final_summary.stats 
		rm -fr barcode_arrangement
    """
}


/*
*  Clean files
*/
process cleanFile {
    tag "${id}"
    label (params.LABEL)
    
    input:
    tuple val(id), path(file_to_remove)
    val(file_to_wait1)
	val(extension)

	when: params.saveSpace == "YES"
    
    script:
    """
		for i in *${extension}; do rm \$(readlink -f \$i); done
    """
}


/*
*  Concatenate FastQ files
*/
process concatenateFastQFiles {
    tag "${idfile}"
    label (params.LABEL)
    publishDir(params.OUTPUT, mode:'copy') 

    input:
    tuple val(idfile), path(demultifq)

    output:
    tuple val(idfile), path("${idfile}.fq.gz") 
    

    script:
    """
		cat *.fastq.gz  >> ${idfile}.fq.gz
    """
}

/*
*  Perform QC on fast5 files
*/

process MinIONQC {
    tag "${folder_name}"
    label (params.LABEL)
    container 'biocorecrg/mopprepr:0.7'
    errorStrategy 'ignore'
    publishDir(params.OUTPUT, mode:'copy', pattern: '*.stats') 

    
    input:
    tuple val(folder_name), path("summaries_*") 

    output:
    tuple val(folder_name), path ("${folder_name}_QC"), emit: QC_folder 
    tuple val(folder_name), path ("*_summary.stats"), emit: stats 

    script:
    """
      if [ -f "summaries_" ]; then
	  ln -s summaries_ ${folder_name}_final_summary.stats
	  else 
		  head -n 1 summaries_1 > ${folder_name}_final_summary.stats
	      for i in summaries_*; do grep -v "filename" \$i >> ${folder_name}_final_summary.stats; done
	  fi
      MinIONQC.R -i ${folder_name}_final_summary.stats -o ${folder_name}_QC -q ${params.qualityqc} -p ${task.cpus}
    """
}

/*
*  Perform bam2stats QC 
*/
process bam2stats {
    label (params.LABEL)
    tag "${id}" 
   
    input:
    tuple val(id), path(bamfile)

    output:
    tuple val(id), path ("${id}.stat")
    
    script:
    """
    bam2stats.py ${bamfile} > ${id}.stat
    """
}

/*
*
*/

process AssignReads {
    tag "${id}"
    publishDir(params.OUTPUT, mode:'copy') 
    label (params.LABEL)

    input:
    tuple val(id), path(input)
	val(tool)

    output:
    tuple val(id), path ("${id}.assigned")
    
    script:
    if (tool == "nanocount")
	    """
		samtools view -F 256 ${input} |cut -f 1,3 > ${id}.assigned
    	"""
    else if(tool == "htseq")
    	"""
			samtools view ${input} | awk '{gsub(/XF:Z:/,"",\$NF); print \$1"\t"\$NF}' |grep -v '__' > ${id}.assigned
    	"""
    else 
        error "Invalid alignment mode: ${tool}"
}

/*
*
*/

process countStats {
    tag "${id}"
    label (params.LABEL)
   
    input:
    tuple val(id), path(input)

    output:
    tuple val(id), path ("${id}.count.stats")
    
    script:
	"""
		wc -l ${input} |sed s@.assigned@@g | awk '{print \$2"\t"\$1}' > ${id}.count.stats
    """
}

/*
*  Join AlnStats 
*/
process joinAlnStats {
    label (params.LABEL)
    tag "joining aln stats"
 
    input:
    file "alnqc_*" 

    output:
    path("alnQC_mqc.txt") 
    
    script:
    """
    echo '# id: alnQC
# plot_type: \'table\'
# section_name: \'Alignment QC\' ' > alnQC_mqc.txt
    cat alnqc_* | head -n 1| sed s@#@@g >> alnQC_mqc.txt
    cat alnqc_* | grep -v "#" >> alnQC_mqc.txt
    """
}

/*
*  Join Count Stats 
*/
process joinCountStats {
    tag "joining count stats"
    label (params.LABEL)
  
    input:
    file "stats_*" 

	output:
	path("counts_mqc.txt")
	
	script:
	"""
	echo '# id: Assigned reads
	# plot_type: \'table\'
	# section_name: Assigned counts 
	File name	\'Counts\' ' > counts_mqc.txt 
		cat stats_*  >> counts_mqc.txt 
		"""
} 

 process bam2Cram {
    tag "${idfile}"  
    
    publishDir(params.OUTPUT, mode:'copy') 
    label (params.LABEL)

    input:
    path(reference)
    val(subsampling_val)
    tuple val(idfile), path(aln), path(index)
    
    output:
    file("*.sorted.cram*") optional true 
    
    script:
    def downcmd = ""
    def input = aln
    def perc = 0
    def clean = ""
 	if (subsampling_val != "NO") {
		perc = subsampling_val/100
		downcmd = "samtools view -@ ${task.cpus} -bs ${perc} ${aln} > subsample.bam"
		input = "subsample.bam"
		clean = "rm subsample.bam"
	}
 	"""
 		${downcmd}
		samtools view  -@ ${task.cpus} -C ${input} -T ${reference} -o ${idfile}.sorted.cram
		samtools index -@ ${task.cpus} ${idfile}.sorted.cram
		${clean}
    """
}

// COMMON TO ALL
process checkRef {
    tag "Checking ${ reference }"
    label (params.LABEL)
 
    input:
    path(reference)
    
    output:
    path("reference.fa")
    
    script:
	"""
	if [ `echo ${reference} | grep ".gz"` ]; then 
   		zcat ${reference} > reference.fa
	else 
        ln -s ${reference} reference.fa
	fi
	"""	
}

// MOP_MOD and MOP_TAIL

process indexReference {
    label (params.LABEL)
    container 'biocorecrg/mopmod:0.6'
    tag "Indexing ${ reference }"

    input:
    path(reference)
    
    output:
    tuple path("*.dict"), path ("*.fai")
    
    script:
	"""
	\$PICARD CreateSequenceDictionary R=${reference} O=reference.dict
	samtools faidx ${reference}
	"""
}

/*
* 
*/

/*
* CALC MEAN PER POSITION
*/

process mean_per_pos {

    container 'biocorecrg/mopmod:0.6'
    label (params.LABEL)
    tag "${idsample}" 
	
    input:
    tuple val(idsample), path(event_align) 
    
    output:
    tuple val(idsample), path("*_perpos_median.parquete")


    script:
    
    """
	mean_per_pos.py -i ${event_align} -o `basename ${event_align} .fast5_event_align.tsv.gz` -s 500000
	#gzip *_processed_perpos_median.tsv
    """
}

/*
* CONCAT MEAN PER POS
*/
process concat_mean_per_pos {

    container 'biocorecrg/mopmod:0.6'
    label (params.LABEL)
    tag "${idsample}" 
	
    input:
    tuple val(idsample), path(event_align) 
    
    output:
    tuple val(idsample), path("*")


    script:
    """
	Merging_processed_nanopolish_data.py -i *.parquete -o ${idsample}
    """
}


/*
*
*/

process callVariants {
    tag "${sampleID}" 
    container 'biocorecrg/mopmod:0.6'
    label (params.LABEL)
	
    input:
    tuple val(sampleID), path(alnfile), path(reference), path(dict_index), path(faiidx) 

    output:
    tuple val(sampleID), path("${sampleID}.tsv")
   
    script:
	"""
	samtools view -h ${alnfile} -F 256 | \$SAM2TSV -R ${reference} | cut -f 3 --complement  > ${sampleID}.tsv
	"""
}

process makeEpinanoPlots {
	publishDir params.OUTPUT, mode: 'copy'
	container "biocorecrg/mopnanotail:0.3"
    label (params.LABEL)

    tag {"${sampleIDA}--${sampleIDB} ${mode}"}  
	
    input:
    tuple val(sampleIDA), val(sampleIDB), path(per_site_varA), path(per_site_varB) 
    val(mode)
    
    output:
    path("*.pdf")
       
    script:
	"""
	Rscript --vanilla ${baseDir}/bin/epinano_scatterplot.R ${per_site_varA} ${sampleIDA} ${per_site_varB} ${sampleIDB} ${mode}  
	"""
}

process multiToSingleFast5 {
    container 'biocorecrg/mopmod:0.6'
    label (params.LABEL)

    tag "${idsample}"  
	
    input:
    tuple val(idsample), path(fast5)
    
    output:
    tuple val(idsample), path("${idsample}-single")
       
    script:
	"""
    mkdir ${idsample}-single;
    multi_to_single_fast5 -i ./ -s ./ -t ${task.cpus}; 
    rm ./filename_mapping.txt; 
    mv ./*/*.fast5 ${idsample}-single;
	"""
}

/*
*
*/
process bedGraphToWig {
    container 'biocorecrg/mopmod:0.6'
    tag "${idsample}"  
	
    input:
    path(chromsizes)
    tuple val(idsample), path(bedgraph)
    
    output:
    tuple val(idsample), path("*.bw")
       
    script:
    def ofname = "${bedgraph.baseName}.wig"
	"""
	bedgraph2wig.pl --bedgraph ${bedgraph} --wig ${ofname} --step 1 --compact
	wigToBigWig ${ofname} ${chromsizes} ${bedgraph.baseName}.bw
	rm *.wig
	"""
}

/*
*/
process mergeTomboWigs {
    label (params.LABEL)
    tag "${combID}"  
	publishDir params.OUTPUT, pattern: "*_Tombo_Output.tsv",  mode: 'copy'
	container "biocorecrg/mopmod:0.6"

   input:
    val(strand)
    tuple val(combID), path(coverage), path(covcontrol), path(statistic)

	output:
	path("*_Tombo_Output.tsv.gz") optional true 
	
	script:
	"""
	Merge_Tombo.py ${statistic} ${covcontrol} ${coverage} ${combID}.${strand}
	gzip *_Tombo_Output.tsv
	"""
}
/*
*/
process wigToBigWig {
    label (params.LABEL)
    tag "${id}"  
	container "biocorecrg/mopmod:0.6"
    errorStrategy 'ignore'

   input:
    path(chromsizes)
    tuple val(id), path(bedgraph)

	output:
	tuple val(id), path("*.bw") optional true 
	
	script:
    def ofname = "${bedgraph.baseName}.bw"

	"""
	wigToBigWig ${bedgraph} ${chromsizes} ${ofname}
	"""
}


// MOP_TAIL

process collect_tailfindr_results {
	publishDir params.OUTPUT, pattern: "*_findr.csv",  mode: 'copy'
	tag "${ sampleID }"  
    label (params.LABEL)
	
	input:
	tuple val(sampleID), path("tailfin_*")
	
	output:
    tuple val(sampleID), path("${sampleID}.findr.len"), emit: length 
    tuple val(sampleID), file ("*_findr.csv"), emit: csv 

	script:
	"""
	cat tailfin_* | awk '!(NR>1 && /tail_start/)' >>  ${sampleID}_findr.csv
	awk -F"," '{if (\$5!="" && \$1!="read_id") print \$1"\t"\$5}' ${sampleID}_findr.csv > ${sampleID}.findr.len
	"""

}

/*
*/
process join_nanotail_results {
	container "biocorecrg/mopnanotail:0.2"
    label (params.LABEL)
    tag "joining nanotail results"

    publishDir params.OUTPUT,  mode: 'copy'
	tag { sampleID }  
	
	input:
	tuple val(sampleID), path(nanopol), path(tailfindr), path(genes)
	file(joinScript)
	
	output:
	file("${sampleID}_*")
	
	script:
	"""
	Rscript --vanilla ${joinScript} ${tailfindr} ${nanopol} ${genes} ${sampleID}
	"""
	
}


/*
* Filter bams
*/
process filter_bam {
	tag "${ sampleID }"
    label (params.LABEL)
	
	input:
	file(reference)
	tuple val(sampleID), path(alignment)

	output:
	tuple val(sampleID), path("${sampleID}_filt.bam")

	script:
	"""
    #to keep only mapped reads and remove secondary alignments 
    samtools view -@ {task.cpus} -bF 260 ${alignment} > ${sampleID}_filt.bam 
	"""
}

// MOP CONSENSUS
process indexFasta {
    label (params.LABEL)

    tag "${reference}" 
	
    input:
    path(reference)
    
    output:
    stdout   
       
    script:
	"""
	samtools faidx ${reference}
	cut -f 1,2 ${reference}.fai
	"""
}

process getChromInfo {
    label (params.LABEL)

    tag "${reference}" 
	
    input:
    path(reference)
    
    output:
    path("chrom.sizes")   
       
    script:
	"""
	samtools faidx ${reference}
	cut -f 1,2 ${reference}.fai > chrom.sizes
	"""
}


process nanoConsensus {
	publishDir "${params.OUTPUT}/${sampleIDs}-${chrName}", mode: 'copy'
	container "biocorecrg/mop_consensus:0.1"
    label (params.LABEL)
    errorStrategy 'ignore'

    tag "${sampleIDs} on ${chrName}"  
	
    input:
    path(nanoConScript)
    path(nanoScripts)
    path(reference)
    tuple val(sampleIDs), path(Epi_Sample), path(Epi_IVT), path(NP_Sample), path(NP_IVT), path(Tombo), path(Nanocomp), val(chrName), val(chrStart), val(chrEnd)
    
    output:
    path("*")
       
    script:
	"""
	Rscript --vanilla ${nanoConScript} -Epi_Sample ${Epi_Sample} \
	 -Epi_IVT ${Epi_IVT} \
	 -NP_Sample ${NP_Sample} \
	 -NP_IVT ${NP_IVT} \
	 -Tombo ${Tombo} \
	 -Nanocomp ${Nanocomp} \
	 -chr ${chrName} \
	 -ini_pos ${chrStart} -fin_pos ${chrEnd} \
	 -output ${sampleIDs} \
	 -fasta ${reference} \
	 --nanocomp_stat GMM_logit_pvalue 
	"""
}


/*
* COMMON FUNCTIONS
*/

// Create a channel for tool options
def getParameters(pars_tools_file) {
	def pars_tools = file(pars_tools_file)
	if( !pars_tools.exists() ) exit 1, "Missing tools options config: '$pars_tools'"

	def progPars = [:]
	def allLines  = pars_tools.readLines()

	for( line : allLines ) {
    	def list = line.split("\t")
    	if (list.length <3) {
			 error "ERROR!!! Tool option file has to be tab separated\n" 
		}
    	if (!(list[0] =~ /#/ )) {
			progPars["${list[0]}--${list[1]}"] = list[2].replace("\"", "").replace('$baseDir', "${baseDir}").replace('${baseDir}', "${baseDir}")
    	}  
	}	
	return(progPars)
}

// Create a channel for tool options
def parseFinalSummary(final_summary_file) {
	final_summary = file(final_summary_file)
	def outstring = ""
	if( final_summary.exists() ) {
		def allLines  = final_summary.readLines()

		for( line : allLines ) {
    		def list = line.split("=")
    		if (list[0] == "protocol") {
    			def vals = list[1].split(":")
    			outstring = "--flowcell ${vals[1]} --kit ${vals[2]}"
 		   	}  
		}	
	} else {
    	log.info '***No configuration file. You must specify kit and flowcell in the parameters!!***\n'
	}
	return(outstring)
}


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

def mapIDPairs (ids, values) {
	def combs = ids.combine(values, by:0).map{
		[it[1], it[0], it[2]]
	}.combine(values, by:0).map{
		[it[1], it[0], it[2],  it[3]]
	}
	return(combs)
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