title=RNA-Seq pipeline
date=2014-07-23
type=page
status=published
syntaxhighlighter=yes
~~~~~~

<div class="blg-summary example">
<h3><a href="javascript:void(0)">RNA-Seq pipeline</a></h3>

<p class="text-muted">
	The example below shows how put together a RNAseq pipeline with basic functionality. It maps a collection of read-pairs to a given 
	reference genome and outputs the respective transcript model.
</p>

<script type="syntaxhighlighter" class="brush: groovy">
<![CDATA[
#!/usr/bin/env nextflow

/*
 * Defines pipeline parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
  
/*
 * The reference genome file
 */
genome_file = file(params.genome) 
 
/*
 * Creates the `read_pairs` channel that emits for each read-pair a tuple containing 
 * three elements: the pair ID, the first read-pair file and the second read-pair file 
 */
Channel
    .fromFilePairs( params.reads )                                              
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }   
    .set { read_pairs } 
 
/*
 * Step 1. Builds the genome index required by the mapping process
 */
process buildIndex {
    input:
    file genome_file
     
    output:
    file 'genome.index*' into genome_index
       
    """
    bowtie2-build ${genome_file} genome.index
    """
}
 
/*
 * Step 2. Maps each read-pair by using Tophat2 mapper tool
 */
process mapping {     
    input:
    file genome_file
    file genome_index from genome_index.first()
    set pair_id, file(reads) from read_pairs
 
    output:
    set pair_id, "tophat_out/accepted_hits.bam" into bam
 
    """
    tophat2 genome.index ${reads}
    """
}
 
/*
 * Step 3. Assembles the transcript by using the "cufflinks" 
 * and publish the transcript output files into the `results` folder
 */
process makeTranscript {
    publishDir "results"
    
    input:
    set pair_id, bam_file from bam
     
    output:
    set pair_id, 'transcripts.gtf' into transcripts
 
    """
    cufflinks ${bam_file}
    """
}
 

]]>
</script>
</div>
    

### Try it in your computer 

In order to run this pipeline in your computer you will required: 

* Unix-like operating system 
* Java 7 (or higher)
* Docker 


Install Nextflow entering the following command in the shell terminal:

    $ curl -fsSL get.nextflow.io | bash


Then launch the pipeline execution using this command: 

    $ nextflow run rnatoy -with-docker 

It will automatically download the pipeline [Github repository](https://github.com/nextflow-io/rnatoy) 
and the associated Docker images, thus the first execution can requires few minutes to complete 
depending you network connection. 
