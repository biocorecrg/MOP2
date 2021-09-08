#!/usr/bin/env nextflow


/* 
 * Define the pipeline parameters
 *
 */

// Pipeline version
version = '0.1'

params.help            = false
params.resume          = false

log.info """

╔╦╗┌─┐┌─┐┌┬┐┌─┐┬─┐  ┌─┐┌─┐  ╔═╗╔═╗╦═╗╔═╗╔═╗
║║║├─┤└─┐ │ ├┤ ├┬┘  │ │├┤   ╠═╝║ ║╠╦╝║╣ ╚═╗
╩ ╩┴ ┴└─┘ ┴ └─┘┴└─  └─┘└    ╩  ╚═╝╩╚═╚═╝╚═╝
                                                                                       
====================================================
BIOCORE@CRG NanoDirectRNA. Detection of modification and polyA length (RNA) - N F  ~  version ${version}
====================================================

*****************   Input files    *******************
input_path                              : ${params.input_path}
comparison                              : ${params.comparison}

********** reference has to be the genome *************
reference                               : ${params.reference}
output                                  : ${params.output}

coverage                                : ${params.coverage}

************* tombo and epinano params ****************
tombo_mode                              : ${params.tombo_mode}
tombo_lsc_opt                           : ${params.tombo_lsc_opt}
tombo_msc_opt                           : ${params.tombo_msc_opt}
epinano_opt                             : ${params.epinano_opt}

************* nanopolish and nanocompore params ***********
nanopolish_opt                          : ${params.nanopolish_opt}
nanocompore_opt                         : ${params.nanocompore_opt}


email                                   : ${params.email}
"""

// Help and avoiding typos
if (params.help) exit 1
if (params.resume) exit 1, "Are you making the classical --resume typo? Be careful!!!! ;)"

// check input files
reference = file(params.reference)
if( !reference.exists() ) exit 1, "Missing reference file: ${reference}!"
//config_report = file("$baseDir/config.yaml")
//if( !config_report.exists() ) exit 1, "Missing config.yaml file!"
logo = file("$baseDir/../docs/logo_small.png")

model_folder = file("$baseDir/models/")
if( !model_folder.exists() ) exit 1, "Missing folders with EpiNano's models!"
joinScript = file("$baseDir/bin/join.r")
mergeTomboScript = file("$baseDir/bin/Merge_Tombo_wigs_MoP_v2.R")

tombo_mode      = params.tombo_mode
tombo_msc_opt   = params.tombo_msc_opt
tombo_lsc_opt   = params.tombo_lsc_opt
epinano_opt     = params.epinano_opt

if( tombo_mode != "level_sample_compare" &&  tombo_mode != "model_sample_compare" && tombo_mode != "both") exit 1, "please specify a valid tombo_mode: level_sample_compare, model_sample_compare or both"

// Output folders
//outputtombo       = "${params.output}/Tombo_Level_Sample_Compare"
outputtombo       = "${params.output}/Tombo"
//outputtombo_comp  = "${params.output}/Tombo_Model_Sample_Compare"
outputEpinano     = "${params.output}/Epinano"
outputCombined    = "${params.output}/Comb_mod"
outputNanocompore = "${params.output}/Nanocomp"
outputNanopolish  = "${params.output}/Nanopolish"
//outputReport   = file("${outputMultiQC}/multiqc_report.html")

/*
* move old multiQCreport

if( outputReport.exists() ) {
  log.info "Moving old report to multiqc_report.html multiqc_report.html.old"
  outputReport.moveTo("${outputMultiQC}/multiqc_report.html.old")
}
*/

compfile = file(params.comparison)
if( !compfile.exists() ) exit 1, "Missing comparison file: ${compfile}. Specify path with --comparisons"

/*
 * Creates the channels with comparisons
 */
 Channel
    .from(compfile.readLines())
    .map { line ->
        list = line.split("\t")
        if (list[0]!= "") {
            def sampleID = list[0]
            def ctrlID = list[1]
            [ sampleID, ctrlID ]
        }
    }
    .into{ id_to_tombo_fast5; id_to_tombo_idx; id_to_epinano; id_for_resquiggling; id_for_index_np; id_for_resquiggling_np; id_for_nanocompore}

/*
 * Creates the channels with samples and KOs for EpiNano
 */
 Channel
    .from(compfile.readLines())
    .map { line ->
        list = line.split("\t")
        if (list[0]!= "") {
            list[0] 
        }
    }
    .set{ samples_for_epinano_filtering }

 Channel
    .from(compfile.readLines())
    .map { line ->
        list = line.split("\t")
        if (list[1]!= "") {
            list[1] 
        }
    }
    .set{ ko_for_epinano_filtering }

id_for_resquiggling.flatten().unique().map {
    [it, file("${params.input_path}/${it}/fast5_files/*.fast5")]
}.transpose().set{fast5_for_resquiggling}

id_for_index_np.flatten().unique().map {
    [it, file("${params.input_path}/${it}/fast5_files/")]
}.transpose().set{fast5_fold_for_np_indexing}

fast5_fold_for_np_indexing.map {
    [it[0],  file(it[1]), file("${params.input_path}/${it[0]}/fastq_files/*.fq.gz"), file("${params.input_path}/${it[0]}/QC_files/final_summary.stats")]
}.transpose().set{data_for_nanopolishindex}

id_for_resquiggling_np.flatten().unique().map {
    [it, file("${params.input_path}/${it}/fast5_files/*.fast5")]
}.transpose().set{fast5_for_resquiggling_np}

fast5_for_resquiggling_np.map {
    [it[0], file(it[1]), file("${params.input_path}/${it[0]}/alignment/*.bam"), file("${params.input_path}/${it[0]}/alignment/*.bai"), file("${params.input_path}/${it[0]}/fastq_files/*.fq.gz"), file("${params.input_path}/${it[0]}/QC_files/final_summary.stats")]
}.transpose().set{data_for_nop_eventalign}


id_to_epinano.flatten().unique().map {
    [it, file("${params.input_path}/${it}/alignment/*.bam")]
}.transpose().set{bams_for_variant_calling}

/*
* Perform resquiggling
*/

process resquiggling_with_tombo {
    tag {"${idsample}-${fast5.simpleName}"}  
    label 'big_mem_cpus'
	
    input:
    set idsample, file(fast5) from fast5_for_resquiggling
    file(reference)
    
    output:
    file ("*.resquiggle.failed.tsv") into failed_resquiggles
    set idsample, file ("${idsample}-${fast5.simpleName}") into singlefast5
    set idsample, file (".*.tombo.index") into tombo_indexes

    script:
    def infolder = "${idsample}-${fast5.simpleName}"
    
    """ 
	#from multifast5 to singlefast5
    mkdir ${infolder};
    multi_to_single_fast5 -i ./ -s ./ -t ${task.cpus}; 
    rm ./filename_mapping.txt; 
    mv ./*/*.fast5 ${infolder};
    # resquiggling
    tombo resquiggle ${infolder} ${reference} --rna --processes ${task.cpus} --overwrite --failed-reads-filename ${infolder}.resquiggle.failed.tsv 
    """
}


/*
* Perform resquiggling with nanopolish
*/

process indexing_with_nanopolish {
    label 'single_cpu'
    tag {"${idsample}"}  
	container "biocorecrg/mopnanopolish:0.2"
	
    input:
    set idsample, file(fast5_folder), file(fastq), file(seqsummary) from data_for_nanopolishindex
    file(reference)
    
    output:
    set idsample, file ("${fastq}.*") into nanopolish_idx

    script:
    
    """ 
      nanopolish index -s ${seqsummary} -d ${fast5_folder} ${fastq}
    """
}


/*
*/
process eventalign_with_nanopolish {
    label 'big_mem_cpus'
    tag {"${idsample}-${fast5_folder.simpleName}"}  
	container "biocorecrg/mopnanopolish:0.2"
	
    input:
    set idsample, file(fast5_folder), file(bam), file(bai), file(fastq), file(seqsummary), file("*") from data_for_nop_eventalign.combine(nanopolish_idx,by: 0)
    file(reference)
    
    output:
    set idsample, file ("*_event_align.tsv.gz") into np_eventalign_folders

    script:
    """ 
    mkdir fast5_files/
    cd fast5_files/; ln -s ../*.fast5 .; cd ../
    nanopolish eventalign -t ${task.cpus} --reads ${fastq} --bam ${bam} --genome ${reference} --samples --print-read-names --scale-events --samples | pigz -p ${task.cpus} > ${idsample}_${fast5_folder}_event_align.tsv.gz
    """
}


/*
*/

process cat_collapsed_nanopolish {
	publishDir outputNanopolish, pattern: "*_combined.eventalign.tsv.gz",  mode: 'copy'
    publishDir outputNanopolish, pattern: "*_processed_perpos_median.tsv.gz", mode: 'copy'
	container "biocorecrg/mopnanopolish:0.2"

	label 'big_long_mem_cpus'
	tag {"${idsample}"}  
	
    input:
    set idsample, file("event_align_*") from np_eventalign_folders.groupTuple() 
    
    output:
    set idsample, file("${idsample}_event_collapsed_align.tsv") into np_event_collapsed 
    file("${idsample}_combined.eventalign.tsv.gz") 
    file("${idsample}_processed_perpos_median.tsv.gz")


    script:
    """
	zcat event_align* | awk '!(/^contig/ && NR>1)' | tee   >(pigz -p ${task.cpus} -9 - > ${idsample}_combined.eventalign.tsv.gz) | NanopolishComp Eventalign_collapse -t ${task.cpus} -o ${idsample}_event_collapsed_align.tsv

	mean_per_pos.py -i ${idsample}_combined.eventalign.tsv.gz -o ${idsample} -s 500000
gzip ${idsample}_processed_perpos_median.tsv
         
    """
}



/*
* Group data together with indexes
*/
np_event_collapsed.groupTuple().into{grouped_event_A; grouped_event_B}
id_for_nanocompore.combine(grouped_event_A, by: 0).map {
	[ it[1], it[0], it[2] ]
}.combine(grouped_event_B, by: 0)map {
	[ "${it[1]}--${it[0]}", it[2], it[3]]
}.set{grouped_event_for_nanocompore}


/*
* detect modification 
*/

process calcNanoCompore {
    label 'big_long_mem_cpus'
    tag {"${combID}"}  
	publishDir outputNanocompore,  mode: 'copy'
    container "quay.io/biocontainers/nanocompore:1.0.0rc3.post1--py38_0"

    input:
    file(reference)
    set val(combID), file(tsv_file_A), file(tsv_file_B) from grouped_event_for_nanocompore
    
    output:
    file("${combID}_nanocompres")
    
    script:
	def folder_names = "${combID}".split("--")
	def folder_name_A = folder_names[0]
	def folder_name_B = folder_names[1]
	"""
nanocompore sampcomp --nthreads ${task.cpus}\
    --file_list1 ${tsv_file_A}/out_eventalign_collapse.tsv \
    --file_list2 ${tsv_file_B}/out_eventalign_collapse.tsv \
    --label1 ${folder_name_A} \
    --label2 ${folder_name_B} \
    --fasta ${reference} \
    --outpath ./${combID}_nanocompres/ \
    ${params.nanocompore_opt} --pvalue_thr 1 --logit --comparison_methods GMM,KS,MW,TT --overwrite 
	"""
}





/*
* Group data together with indexes
*/
singlefast5.groupTuple().into{grouped_single_fast5_A; grouped_single_fast5_B}
id_to_tombo_fast5.combine(grouped_single_fast5_A, by: 0).map {
	[ it[1], it[0], it[2] ]
}.combine(grouped_single_fast5_B, by: 0)map {
	[ "${it[1]}--${it[0]}", it[2], it[3]]
}.set{fast5_for_tombo_modifications}


tombo_indexes.groupTuple().into{grouped_indexes_A; grouped_indexes_B}
id_to_tombo_idx.combine(grouped_indexes_A, by: 0).map {
	[ it[1], it[0], it[2] ]
}.combine(grouped_indexes_B, by: 0)map {
	[ "${it[1]}--${it[0]}", it[2], it[3] ]
}.set{idx_for_tombo_modifications}

fast5_for_tombo_modifications.combine(idx_for_tombo_modifications, by: 0).into{data_for_tombo_modifications; data_for_tombo_modifications_compare; luca}

/*
* detect modification 
*/

process getModificationsWithTombo_LevelSampleCompare {
    label 'big_long_mem_cpus'
    tag {"${combID}"}  
         
    input:
    file(reference)
    //file(mergeTomboScript)
    set val(combID), file(fast5s_A), file(fast5s_B), file(index_A), file(index_B) from data_for_tombo_modifications

    when: 
    tombo_mode == "level_sample_compare" || tombo_mode == "both"
    
    output:
    //file ("*_Tombo_Output.tsv") into sign_tombo_regions
    set val("${combID}"), val("lsc_plus"), file("*.statistic.plus.wig"), file("*.coverage.sample.plus.wig"), file("*.coverage.control.plus.wig") optional true into wig_plus_tombo_lsc
    set val("${combID}"), val("lsc_minus"), file("*.statistic.minus.wig"), file("*.coverage.sample.minus.wig"), file("*.coverage.control.minus.wig") optional true into wig_minus_tombo_lsc

    script:
	def reference_cmd = unzipFile(reference, "reference.fa")
	def folder_names = "${combID}".split("--")
	def folder_name_A = folder_names[0]
	def folder_name_B = folder_names[1]
	"""
	${reference_cmd}
	mkdir ${folder_name_A} ${folder_name_B}
	mv ${fast5s_A} ${folder_name_A}
	mv ${index_A} ${folder_name_A}
	mv ${fast5s_B} ${folder_name_B}
	mv ${index_B} ${folder_name_B}

	tombo detect_modifications level_sample_compare \
	--fast5-basedirs ${folder_name_A}/* \
    --alternate-fast5-basedirs ${folder_name_B}/* \
    --processes ${task.cpus} ${params.tombo_lsc_opt} \
    --statistics-file-basename ${folder_name_A}_${folder_name_B}_model_sample_compare --store-p-value
    
	tombo text_output browser_files --fast5-basedirs ${folder_name_A}/* \
	--control-fast5-basedirs ${folder_name_B}/* \
	--browser-file-basename ${folder_name_A}_${folder_name_B} \
	--statistics-filename ${folder_name_A}_${folder_name_B}_model_sample_compare.tombo.stats \
	--file-types {'coverage','statistic'}

	bedgraph2wig.pl --bedgraph ${folder_name_A}_${folder_name_B}.coverage.sample.plus.bedgraph --wig ${folder_name_A}_${folder_name_B}.coverage.sample.plus.wig --step 1 --compact
	bedgraph2wig.pl --bedgraph ${folder_name_A}_${folder_name_B}.coverage.sample.minus.bedgraph --wig ${folder_name_A}_${folder_name_B}.coverage.sample.minus.wig --step 1 --compact
	bedgraph2wig.pl --bedgraph ${folder_name_A}_${folder_name_B}.coverage.control.plus.bedgraph --wig ${folder_name_A}_${folder_name_B}.coverage.control.plus.wig --step 1 --compact
	bedgraph2wig.pl --bedgraph ${folder_name_A}_${folder_name_B}.coverage.control.minus.bedgraph --wig ${folder_name_A}_${folder_name_B}.coverage.control.minus.wig --step 1 --compact

	rm reference.fa
	"""
}

/*
*/
process getModificationsWithTombo_ModelSampleCompare {
    label 'big_long_mem_cpus'
    tag {"${combID}"}  
            
    input:
    file(reference)
    //file(mergeTomboScript)
    set val(combID), file(fast5s_A), file(fast5s_B), file(index_A), file(index_B) from data_for_tombo_modifications_compare

	when: 
    tombo_mode == "both" || tombo_mode == "model_sample_compare"
    
    output:
    set val("${combID}"), val("msc_plus"), file("*.dampened_fraction_modified_reads.plus.wig"), file("*.coverage.sample.plus.wig"), file("*.coverage.control.plus.wig") optional true  into wig_plus_tombo_msc
    set val("${combID}"), val("msc_minus"), file("*.dampened_fraction_modified_reads.minus.wig"), file("*.coverage.sample.minus.wig"), file("*.coverage.control.minus.wig") optional true  into wig_minus_tombo_msc

    script:
	def reference_cmd = unzipFile(reference, "reference.fa")
	def folder_names = "${combID}".split("--")
	def folder_name_A = folder_names[0]
	def folder_name_B = folder_names[1]

	"""
	${reference_cmd}
	mkdir ${folder_name_A} ${folder_name_B}
	mv ${fast5s_A} ${folder_name_A}
	mv ${index_A} ${folder_name_A}
	mv ${fast5s_B} ${folder_name_B}
	mv ${index_B} ${folder_name_B}

       tombo detect_modifications model_sample_compare --fast5-basedirs ${folder_name_A}/* \
       --control-fast5-basedirs ${folder_name_B}/* \
       --processes ${task.cpus} ${params.tombo_msc_opt} \
       --statistics-file-basename ${folder_name_A}_${folder_name_B}_model_sample_compare
      
       tombo text_output browser_files --fast5-basedirs ${folder_name_A}/* \
       --control-fast5-basedirs ${folder_name_B}/* \
       --browser-file-basename ${folder_name_A}_${folder_name_B}.features \
       --statistics-filename ${folder_name_A}_${folder_name_B}_model_sample_compare.tombo.stats \
       --file-types 'dampened_fraction' 'coverage'
       
        bedgraph2wig.pl --bedgraph ${folder_name_A}_${folder_name_B}.features.coverage.sample.plus.bedgraph \
        --wig ${folder_name_A}_${folder_name_B}.features.coverage.sample.plus.wig --step 1 --compact

        bedgraph2wig.pl --bedgraph ${folder_name_A}_${folder_name_B}.features.coverage.control.plus.bedgraph \
        --wig ${folder_name_A}_${folder_name_B}.features.coverage.control.plus.wig --step 1 --compact

        bedgraph2wig.pl --bedgraph ${folder_name_A}_${folder_name_B}.features.coverage.sample.minus.bedgraph \
        --wig ${folder_name_A}_${folder_name_B}.features.coverage.sample.minus.wig --step 1 --compact

        bedgraph2wig.pl --bedgraph ${folder_name_A}_${folder_name_B}.features.coverage.control.minus.bedgraph \
        --wig ${folder_name_A}_${folder_name_B}.features.coverage.control.minus.wig --step 1 --compact

	rm reference.fa
	"""
}



/*
*/
process mergeTomboWigs {
    label 'big_long_mem_cpus'
    tag {"${combID}"}  
	publishDir outputtombo, pattern: "*_Tombo_Output.tsv",  mode: 'copy'
	container "biocorecrg/mopnanotail:0.2"

   input:
    file(mergeTomboScript)
    set val(combID), val(strand), file(statistic), file(coverage), file(covcontrol) from wig_plus_tombo_msc.mix(wig_minus_tombo_msc).mix(wig_plus_tombo_lsc).mix(wig_minus_tombo_lsc)

	output:
	file("*_Tombo_Output.tsv") optional true 
	
	script:
	"""
      Rscript --vanilla ${mergeTomboScript} \
      -stats_wig ${statistic} \
      -Cov_WT ${coverage} \
      -Cov_Control ${covcontrol} \
      -output ${combID}.${strand}
	"""

}

/*
* Perform preprocessing for Epinano
*/

process IndexReferenceForEpinano {

    input:
    file(reference)
    
    output:
    set file("reference.fa"), file("*.dict"), file ("*.fai") into indexes
    
    script:
	"""
	if [ `echo ${reference} | grep ".gz"` ]; then 
   		zcat ${reference} > reference.fa
	else 
        ln -s ${reference} reference.fa
	fi
	\$PICARD CreateSequenceDictionary R=reference.fa O=reference.dict
	samtools faidx reference.fa
	"""
}

/*
* Calling variants for Epinano
*/

process CallVariantsForEpinano {

    tag {"${sampleID}"}  
	
    input:
    set val(sampleID), file(alnfile) from bams_for_variant_calling
    set file(reference), file(dict_index), file(faiidx) from indexes

    output:
    set sampleID, file("${sampleID}.tsv") into variants_for_frequency
   
    script:
	"""
	samtools view -h ${alnfile} -F 256 | \$SAM2TSV -R ${reference} | cut -f 3 --complement  > ${sampleID}.tsv
	"""
}

/*
* 
*/

process calcVarFrequenciesForEpinano {
	publishDir outputEpinano,  pattern: "*.csv.gz", mode: 'copy'
	container "biocorecrg/mopepinano:0.1"

    tag {"${sampleID}"}  
    label 'big_mem_cpus'
	
    input:
    set val(sampleID), file(tsvfile) from variants_for_frequency
    
    output:
    set val(sampleID), file("*per_site_var.5mer.csv.gz") into per_site_vars
    file("*.csv.gz")
   
    script:
	"""
	TSV_to_Variants_Freq.py3 -f ${tsvfile} -t ${task.cpus}
	for i in *.csv; do gzip \$i; done
	"""
}

per_site_vars.map{
	[ it[0], it[1].splitText( by: 1000000, keepHeader:true, compress:true, file:true ) ]
}.transpose().set{per_site_vars_splitted}


/*
* functions
*/
// Get the name from the folder
def getFolderName(sample) {
   folder_info = sample.toString().tokenize("/")
   return folder_info[-2]
}




// make named pipe 
def unzipBash(filename) { 
    cmd = filename.toString()
    if (cmd[-3..-1] == ".gz") {
    	cmd = "<(zcat ${filename})"
    }
    return cmd
}

def unzipFile(zipfile, filename) {
    cmd = "ln -s ${zipfile} ${filename}"
    filestring = zipfile.toString()
    if (filestring[-3..-1] == ".gz") {
    	cmd = "zcat ${zipfile} > ${filename}"
    }
    return cmd	
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
        NanoMod module's execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        Error report: ${workflow.errorReport ?: '-'}
        """
        .stripIndent()

        sendMail(to: params.email, subject: "Master of Pore execution", body: msg)
    }
}
