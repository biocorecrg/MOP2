process extracting_demultiplexed_fastq {
	label 'basecall_cpus'
    tag { idfile }
	//publishDir outputFastq,  mode: 'copy'
			
	input:
	tuple val(idfile), file(demux), file(fastq) 
	
	
	output:
	tuple val(idfile), file ("*.fastq.gz")

	script:
	"""
		extract_sequence_from_fastq.py ${demux} ${fastq}
		for i in *.fastq; do gzip \$i; done
	"""
}

process extracting_demultiplexed_fast5 {
	label 'more_cpus'
	container 'lpryszcz/deeplexicon:1.2.0'

    tag { idfile }
	//publishDir outputFast5,  mode: 'copy'
		
	input:
	tuple val(idfile), file("demux_*"), file("*")

	output:
	file("*")
	
	script:
	"""
	cat demux_* | grep -v ReadID >> dem.files
	awk '{print \$2 > \$3".list" }' dem.files
	for i in *.list; do mkdir `basename \$i .list`; fast5_subset --input ./ --save_path `basename \$i .list`/ --read_id_list \$i --batch_size 4000 -t ${task.cpus}; done 
	rm *.list
	rm */filename_mapping.txt
	rm dem.files 
	"""
} 

/*
*  Concatenate FastQ files
*/
process concatenateFastQFiles {
    tag {idfile} 

	//publishDir outputFastq, pattern: "*.fq.gz",  mode: 'copy'

    input:
    tuple val(idfile), file(demultifq)

    output:
    tuple val(idfile), file("${idfile}.fq.gz") 
    

    script:
    """
		cat *.fastq.gz  >> ${idfile}.fq.gz
    """
}

