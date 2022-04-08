title=RNA-Seq pipeline
date=2014-07-23
type=page
status=published
syntaxhighlighter=yes
~~~~~~

<div class="blg-summary example">
<h3><a href="javascript:void(0)">RNA-Seq pipeline</a></h3>

<p class="text-muted">
    This example shows how to put together a basic RNAseq pipeline. It maps a collection of read-pairs to a given reference genome and outputs the respective transcript model.
</p>

<script type="syntaxhighlighter" class="brush: groovy">
<![CDATA[
#!/usr/bin/env nextflow

/*
 * The following pipeline paramemters specify the refence genomes
 * and read pairs and can be provided as command line options
 */
params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.annot = "$baseDir/data/ggal/ggal_1_48850000_49020000.bed.gff"
params.genome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = 'results'

workflow {
    /*
     * Create the `ch_read_pairs` channel, which emits tuples containing three elements:
     * the pair ID, the first read-pair file and the second read-pair file.
     */
    Channel
        .fromFilePairs( params.reads )
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { ch_read_pairs }

    /*
     * Step 1. Build the genome index required by the mapping process
     */
    buildIndex(params.genome)

    /*
     * Step 2. Map each read-pair by using Tophat2 mapper tool
     */
    mapping(
        params.genome,
        params.annot,
        buildIndex.out.ch_index,
        ch_read_pairs)

    /*
     * Step 3. Assemble the transcript using the "cufflinks" tool
     */
    makeTranscript(params.annot, mapping.out.ch_bam)
}

process buildIndex {
    tag "${genome.baseName}"

    input:
    path genome

    output:
    path 'genome.index*', emit: ch_index

    """
    bowtie2-build --threads ${task.cpus} ${genome} genome.index
    """
}

process mapping {
    tag "${pair_id}"

    input:
    path genome
    path annot
    path index
    tuple val(pair_id), path(reads)

    output:
    set pair_id, "accepted_hits.bam", emit: ch_bam

    """
    tophat2 -p ${task.cpus} --GTF ${annot} genome.index ${reads}
    mv tophat_out/accepted_hits.bam .
    """
}

process makeTranscript {
    tag "${pair_id}"
    publishDir params.outdir, mode: 'copy'

    input:
    path annot
    tuple val(pair_id), path(bam_file)

    output:
    tuple val(pair_id), path('transcript_*.gtf')

    """
    cufflinks --no-update-check -q -p ${task.cpus} -G ${annot} ${bam_file}
    mv transcripts.gtf transcript_${pair_id}.gtf
    """
}
]]>
</script>
</div>


### Try it in your computer

To run this pipeline on your computer, you will need:

* Unix-like operating system
* Java 11 (or higher)
* Docker

Install Nextflow by entering the following command in the terminal:

    $ curl -fsSL get.nextflow.io | bash

Then launch the pipeline with this command:

    $ nextflow run rnatoy -with-docker

It will automatically download the pipeline [Github repository](https://github.com/nextflow-io/rnatoy) and the associated Docker images, thus the first execution may take a few minutes to complete depending on your network connection.
