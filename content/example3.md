title=BLAST pipeline
date=2014-07-23
type=page
status=published
syntaxhighlighter=yes
~~~~~~

<div class="blg-summary example">
<h3><a href="javascript:void(0)">BLAST pipeline</a></h3>

<p class="text-muted">
	The example below splits a FASTA file into chunks and executes for each of them a BLAST query in a parallel manner.
	Then, all the sequences for the top hits are collected and merged to a single result file.
</p>

<script type="syntaxhighlighter" class="brush: groovy">
<![CDATA[
#!/usr/bin/env nextflow

/*
 * Defines the pipeline inputs parameters (giving a default value for each for them) 
 * Each of the following parameters can be specified as command line options
 */
params.query = "$baseDir/data/sample.fa"
params.db = "$baseDir/blast-db/pdb/tiny"
params.out = "result.txt"
params.chunkSize = 100 

db_name = file(params.db).name
db_dir = file(params.db).parent

/* 
 * Given the query parameter creates a channel emitting the query fasta file(s), 
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunkSize'.
 * Finally assign the result channel to the variable 'fasta_ch' 
 */
Channel
    .fromPath(params.query)
    .splitFasta(by: params.chunkSize, file:true)
    .set { fasta_ch }

/* 
 * Executes a BLAST job for each chunk emitted by the 'fasta_ch' channel 
 * and creates as output a channel named 'top_hits' emitting the resulting 
 * BLAST matches  
 */
process blast {
    input:
    path 'query.fa' from fasta_ch
    path db from db_dir

    output:
    file 'top_hits' into hits_ch

    """
    blastp -db $db/$db_name -query query.fa -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits
    """
}

/* 
 * Each time a file emitted by the 'top_hits' channel an extract job is executed 
 * producing a file containing the matching sequences 
 */
process extract {
    input:
    path 'top_hits' from hits_ch
    path db from db_dir

    output:
    file 'sequences' into sequences_ch

    """
    blastdbcmd -db $db/$db_name -entry_batch top_hits | head -n 10 > sequences
    """
}

/* 
 * Collects all the sequences files into a single file 
 * and prints the resulting file content when complete 
 */ 
sequences_ch
    .collectFile(name: params.out)
    .view { file -> "matching sequences:\n ${file.text}" }

]]>
</script>
</div>

### Try it in your computer 

In order to run this pipeline in your computer you will required: 

* Unix-like operating system 
* Java 8 (or higher)
* Docker 

Install Nextflow entering the following command in the shell terminal:

    $ curl -fsSL https://get.nextflow.io | bash


Then launch the pipeline execution using this command: 

    $ ./nextflow run blast-example -with-docker 

It will automatically download the pipeline [Github repository](https://github.com/nextflow-io/blast-example) 
and the associated Docker images, thus the first execution can requires few minutes to complete 
depending you network connection. 
