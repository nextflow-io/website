---
title: RNA-Seq pipeline
layout: "@layouts/MarkdownPage.astro"
---

<div class="blg-summary example">
<h3>RNA-Seq pipeline</h3>

<p class="text-muted">
    This example shows a basic RNA-Seq pipeline.
</p>

```groovy
#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Input data
params.reads = "${workflow.projectDir}/data/ggal/ggal_gut_{1,2}.fq"

// Reference file
params.transcriptome = "${workflow.projectDir}/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"

// Output directory
params.outdir = "results"

/*
 * Index reference transcriptome file
 */
process INDEX {
    tag "$transcriptome.simpleName"
    container "community.wave.seqera.io/library/salmon:1.10.3--482593b6cd04c9b7"
    conda "bioconda::salmon=1.10.3"

    input:
    path transcriptome

    output:
    path 'index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}

/*
 * Generate FastQC reports
 */
process FASTQC {
    tag "FASTQC on $sample_id"
    publishDir params.outdir, mode:'copy'
    container "community.wave.seqera.io/library/fastqc:0.12.1--5cfd0f3cb6760c42"
    conda "bioconda::fastqc:0.12.1"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

/*
 * Quantify reads
 */
process QUANT {
    tag "$pair_id"
    publishDir params.outdir, mode:'copy'
    container "community.wave.seqera.io/library/salmon:1.10.3--482593b6cd04c9b7"
    conda "bioconda::salmon=1.10.3"

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

/*
 * Generate MultiQC report
 */
process MULTIQC {
  publishDir params.outdir, mode:'copy'
    container "community.wave.seqera.io/library/multiqc:1.24.1--789bc3917c8666da"
    conda "bioconda::multiqc:1.24.1"

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

workflow {

    // Paired reference data
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )

    // Index reference transcriptome file
    INDEX(params.transcriptome)

    // Generate FastQC reports
    FASTQC(read_pairs_ch)

    // Quantify reads
    QUANT(INDEX.out, read_pairs_ch)

    // Generate MultiQC report
    MULTIQC(QUANT.out.mix(FASTQC.out).collect())
}
```

</div>

### Synopsis

This example shows a basic Nextflow pipeline consisting of four processes. The `INDEX` process indexes a reference transcriptome file. The `FASTQC` process creates reports for the input fastq files. The `QUANT` process takes the indexed transcriptome and input fastq files and quantifies the reads. The `MULTIQC` process collects the output from the `QUANT` and `FASTQC` processes and generates a html report.

### Try it

This pipeline is available on the [nextflow-io/rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf/tree/example) GitHub repository.

An active internet connection and Docker are required for Nextflow to download the pipeline and the necessary Docker images to run the pipeline within containers. The data used by this pipeline is stored on the GitHub repository and will download automatically.

To try this pipeline:

1. Follow the [Nextflow installation guide](https://www.nextflow.io/docs/latest/install.html#install-nextflow) to install Nextflow.
2. Follow the [Docker installation guide](https://docs.docker.com/get-started/get-docker/) to install Docker.
3. Launch the `example` branch of the pipeline:

   nextflow run nextflow-io/rnaseq-nf -r example

**NOTE**: The main branch of the `rnaseq-nf` pipeline on GitHub is under active development and differs from the example shown above. The `rnaseq-nf` pipeline will use Docker to manage software dependencies by default.
