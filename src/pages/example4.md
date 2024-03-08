---
title: RNA-Seq pipeline
layout: @layouts/Page.astro
---


<div class="blg-summary example">
<h3><a href="javascript:void(0)">RNA-Seq pipeline</a></h3>

<p class="text-muted">
    This example shows how to put together a basic RNA-Seq pipeline. It maps a collection of read-pairs to a given reference genome and outputs the respective transcript model.
</p>

<script type="syntaxhighlighter" class="brush: groovy">
<![CDATA[
#!/usr/bin/env nextflow

/*
 * The following pipeline parameters specify the reference genomes
 * and read pairs and can be provided as command line options
 */
params.reads = "$baseDir/data/ggal/ggal_gut_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = "results"

workflow {
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )

    INDEX(params.transcriptome)
    FASTQC(read_pairs_ch)
    QUANT(INDEX.out, read_pairs_ch)
}

process INDEX {
    tag "$transcriptome.simpleName"

    input:
    path transcriptome

    output:
    path 'index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}

process FASTQC {
    tag "FASTQC on $sample_id"
    publishDir params.outdir

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    fastqc.sh "$sample_id" "$reads"
    """
}

process QUANT {
    tag "$pair_id"
    publishDir params.outdir

    input:
    path index
    tuple val(pair_id), path(reads)

    output:
    path pair_id

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
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

    $ nextflow run rnaseq-nf -with-docker

It will automatically download the pipeline [GitHub repository](https://github.com/nextflow-io/rnaseq-nf) and the associated Docker images, thus the first execution may take a few minutes to complete depending on your network connection.

__NOTE__: To run this example with versions of Nextflow older than 22.04.0, you must include the `-dsl2` flag with `nextflow run`.
