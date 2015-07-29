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
 * Defines the pipeline inputs parameters (giving a default value for each for them) 
 * Each of the following parameters can be specified as command line options
 */
params.pair1 = "$baseDir/data/ggal/*_1.fq"
params.pair2 = "$baseDir/data/ggal/*_2.fq"
params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
 
/*
 * Creates a channel from the path where the read files are located, 
 * then transforms each file to a pair containing a name prefix and the file itself. 
 * Finally assign this channel to a variable named 'reads1'
 */
Channel
    .fromPath( params.pair1 )
    .map {  path -> tuple(path.baseName[0..-2], path) }
    .set { reads1 }
  
/*
 * as above for pair2 reads 
 */
Channel
    .fromPath( params.pair2 )
    .map {  path -> tuple(path.baseName[0..-2], path) }
	.set { reads2 }
	     
/*
 * Matches the pairs having the same 'prefix' emitted by the channels "read1" and "read2", 
 * for each match emits a new pair containing the expected read-pair files. 
 * Finally assign this channel to a variable named 'read_pairs'
 */
reads1
	.phase(reads2)
	.map { pair1, pair2 -> tuple(pair1[0], pair1[1], pair2[1]) }
	.set { read_pairs } 
	
/*
 * Creates a reference genome file object given the genome string parameter 
 */
genome_file = file(params.genome)
 
 
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
    set pair_id, file(read1), file(read2) from read_pairs
 
    output:
    set pair_id, "tophat_out/accepted_hits.bam" into bam
 
    """
    tophat2 genome.index ${read1} ${read2}
    """
}
 
 
/*
 * Step 3. Assembles the transcript by using the "cufflinks" 
 */
process makeTranscript {
    input:
    set pair_id, bam_file from bam
     
    output:
    set pair_id, 'transcripts.gtf' into transcripts
 
    """
    cufflinks ${bam_file}
    """
}
 
/*
 * Step 4. Collects the transcripts files and print them
 */
transcripts
  .collectFile { tuple("${it[0]}transcript", it[1]) }
  .println { "Transcript model: $it" }

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
