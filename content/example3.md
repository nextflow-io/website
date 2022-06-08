title=BLAST pipeline
date=2014-07-23
type=page
status=published
syntaxhighlighter=yes
~~~~~~

<div class="blg-summary example">
<h3><a href="javascript:void(0)">BLAST pipeline</a></h3>

<p class="text-muted">
    This example splits a FASTA file into chunks and executes a BLAST query for each chunk in parallel. Then, all the sequences for the top hits are collected and merged into a single result file.
</p>

<script type="syntaxhighlighter" class="brush: groovy">
<![CDATA[
#!/usr/bin/env nextflow

/*
 * Defines the pipeline input parameters (with a default value for each one).
 * Each of the following parameters can be specified as command line options.
 */
params.query = "$baseDir/data/sample.fa"
params.db = "$baseDir/blast-db/pdb/tiny"
params.out = "result.txt"
params.chunkSize = 100

db_name = file(params.db).name
db_dir = file(params.db).parent


workflow {
    /*
     * Create a channel emitting the given query fasta file(s).
     * Split the file into chunks containing as many sequences as defined by the parameter 'chunkSize'.
     * Finally, assign the resulting channel to the variable 'ch_fasta'
     */
    Channel
        .fromPath(params.query)
        .splitFasta(by: params.chunkSize, file:true)
        .set { ch_fasta }

    /*
     * Execute a BLAST job for each chunk emitted by the 'ch_fasta' channel
     * and emit the resulting BLAST matches.
     */
    ch_hits = blast(ch_fasta, db_dir)

    /*
     * Each time a file emitted by the 'blast' process, an extract job is executed,
     * producing a file containing the matching sequences.
     */
    ch_sequences = extract(ch_hits, db_dir)

    /*
     * Collect all the sequences files into a single file
     * and print the resulting file contents when complete.
     */
    ch_sequences
        .collectFile(name: params.out)
        .view { file -> "matching sequences:\n ${file.text}" }
}


process blast {
    input:
    path 'query.fa'
    path db

    output:
    path 'top_hits'

    """
    blastp -db $db/$db_name -query query.fa -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits
    """
}


process extract {
    input:
    path 'top_hits'
    path db

    output:
    path 'sequences'

    """
    blastdbcmd -db $db/$db_name -entry_batch top_hits | head -n 10 > sequences
    """
}
]]>
</script>
</div>

### Try it on your computer

To run this pipeline on your computer, you will need:

* Unix-like operating system
* Java 11 (or higher)
* Docker

Install Nextflow by entering the following command in the terminal:

    $ curl -fsSL https://get.nextflow.io | bash

Then launch the pipeline with this command:

    $ ./nextflow run blast-example -with-docker

It will automatically download the pipeline [Github repository](https://github.com/nextflow-io/blast-example) and the associated Docker images, thus the first execution may take a few minutes to complete depending on your network connection.

__NOTE__: To run this example with versions of Nextflow older than 22.04, you must include the `-dsl2` flag with `nextflow run`.
