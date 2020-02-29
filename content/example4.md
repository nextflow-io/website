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
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.annot = "$baseDir/data/ggal/ggal_1_48850000_49020000.bed.gff"
params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = 'results'

/*
 * Create the `read_pairs_ch` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs_ch } 
 
/*
 * Step 1. Builds the genome index required by the mapping process
 */
process buildIndex {
    tag "$genome.baseName"
    
    input:
    path genome from params.genome
     
    output:
    path 'genome.index*' into index_ch
       
    """
    bowtie2-build --threads ${task.cpus} ${genome} genome.index
    """
}
 
/*
 * Step 2. Maps each read-pair by using Tophat2 mapper tool
 */
process mapping {
    tag "$pair_id"
     
    input:
    path genome from params.genome 
    path annot from params.annot
    path index from index_ch
    tuple val(pair_id), path(reads) from read_pairs_ch
 
    output:
    set pair_id, "accepted_hits.bam" into bam_ch
 
    """
    tophat2 -p ${task.cpus} --GTF $annot genome.index $reads
    mv tophat_out/accepted_hits.bam .
    """
}
  
/*
 * Step 3. Assembles the transcript by using the "cufflinks" tool
 */
process makeTranscript {
    tag "$pair_id"
    publishDir params.outdir, mode: 'copy'  
       
    input:
    path annot from params.annot
    tuple val(pair_id), path(bam_file) from bam_ch
     
    output:
    tuple val(pair_id), path('transcript_*.gtf')
 
    """
    cufflinks --no-update-check -q -p $task.cpus -G $annot $bam_file
    mv transcripts.gtf transcript_${pair_id}.gtf
    """
}

]]>
</script>
</div>
    

### Try it in your computer 

In order to run this pipeline in your computer you will required: 

* Unix-like operating system 
* Java 8 (or higher)
* Docker 


Install Nextflow entering the following command in the shell terminal:

    $ curl -fsSL get.nextflow.io | bash


Then launch the pipeline execution using this command: 

    $ nextflow run rnatoy -with-docker 

It will automatically download the pipeline [Github repository](https://github.com/nextflow-io/rnatoy) 
and the associated Docker images, thus the first execution can requires few minutes to complete 
depending you network connection. 
