/*
* Output folders
*/
outputEpinanoFlow    = "${params.output}/epinano_flow"

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
