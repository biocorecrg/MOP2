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

process extracting_demultiplexed_fastq {
	label 'basecall_cpus'
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
	label 'more_cpus'
	container 'lpryszcz/deeplexicon:1.2.0'
    tag "${ idfile }"
    if (params.saveSpace == "YES") publishDir(outputFast5, mode:'move') 
    else publishDir(outputFast5, mode:'copy')    
		
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
    tag { idfile }
    label (params.LABEL)
    if (params.saveSpace == "YES") publishDir(outputFast5, mode:'move') 
    else publishDir(outputFast5, mode:'copy')    

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
    //if (params.saveSpace == "YES") publishDir(outputFastq, mode:'move') 
    publishDir(outputFastq, mode:'copy') 

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
    label 'big_cpus'
    container 'biocorecrg/mopprepr:0.7'
    errorStrategy 'ignore'
    publishDir(outputQual, mode:'copy', pattern: '*.stats') 

    
    input:
    tuple val(folder_name), path("summaries_*") 

    output:
    tuple val(folder_name), path ("${folder_name}_QC")
    tuple val(folder_name), path ("*_summary.stats")

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
    publishDir(outputAssigned, mode:'copy') 

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
   
    input:
    file "stats_*" 

	output:
	path("counts_mqc.txt")
	
	script:
	"""
	echo '# id: Read counts
	# plot_type: \'table\'
	# section_name: Read counts 
	File name	\'Counts\' ' > counts_mqc.txt 
		cat stats_*  >> counts_mqc.txt 
		"""
} 

