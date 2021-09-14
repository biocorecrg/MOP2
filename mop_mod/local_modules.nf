/*
* Output folders
*/
outputEpinanoFlow    = "${params.output}/epinano_flow"
outputTomboFlow      = "${params.output}/tombo_flow"

process checkRef {
    tag "Checking ${ reference }"
 
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

process indexReference {
	label 'more_cpus'
    container 'biocorecrg/mopmod1:0.2'
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

process callVariants {
    tag "${sampleID}" 
    container 'biocorecrg/mopmod1:0.2'
	
    input:
    tuple val(sampleID), path(alnfile) 
    tuple path(reference), path(dict_index), path(faiidx) 

    output:
    tuple val(sampleID), path("${sampleID}.tsv")
   
    script:
	"""
	samtools view -h ${alnfile} -F 256 | \$SAM2TSV -R ${reference} | cut -f 3 --complement  > ${sampleID}.tsv
	"""
}

epinanoScript = file("$baseDir/bin/epinano_scatterplot.R")

process makeEpinanoPlots {
	publishDir outputEpinanoFlow, mode: 'copy'
	container "biocorecrg/mopnanotail:0.3"

    tag {"${sampleIDA}--${sampleIDB} ${mode}"}  
	
    input:
    tuple val(sampleIDA), val(sampleIDB), path(per_site_varA), path(per_site_varB) 
    val(mode)
    
    output:
    path("*.pdf")
       
    script:
	"""
	Rscript --vanilla ${epinanoScript} ${per_site_varA} ${sampleIDA} ${per_site_varB} ${sampleIDB} ${mode}  
	"""
}

process multiToSingleFast5 {
    container 'biocorecrg/mopmod1:0.2'
	label 'more_cpus'

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
    container 'biocorecrg/mopmod1:0.2'
    tag "${idsample}"  
	
    input:
    tuple val(idsample), path(bedgraph)
    
    output:
    tuple val(idsample), path("*.wig")
       
    script:
    def ofname = "${bedgraph.baseName}.wig"
	"""
	bedgraph2wig.pl --bedgraph ${bedgraph} --wig ${ofname} --step 1 --compact
	"""
}

/*
*/
process mergeTomboWigs {
    label 'big_long_mem_cpus'
    tag {"${combID}"}  
	publishDir outputTomboFlow, pattern: "*_Tombo_Output.tsv",  mode: 'copy'
	container "biocorecrg/mopnanotail:0.2"

   input:
    val(strand)
    path(mergeTomboScript)
    tuple val(combID), path(coverage), path(covcontrol), path(statistic)

	output:
	path("*_Tombo_Output.tsv") optional true 
	
	script:
	"""
      Rscript --vanilla ${mergeTomboScript} \
      -stats_wig ${statistic} \
      -Cov_WT ${coverage} \
      -Cov_Control ${covcontrol} \
      -output ${combID}.${strand}
	"""
}
