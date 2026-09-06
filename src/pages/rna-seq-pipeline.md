---
title: RNA-Seq pipeline
layout: "@layouts/ExampleLayout.astro"
---

<div class="blg-summary example">
<h2>RNA-Seq pipeline</h2>

<p class="">
    This pipeline shows an example of a basic for RNA-Seq analysis that performs quality control, transcript quantification, and result aggregation. The pipeline processes paired-end FASTQ files, generates quality control reports with FastQC, quantifies transcripts with Salmon, and produces a unified report with MultiQC.
</p>

```groovy
// Parameter inputs
params.reads = "${baseDir}/data/ggal/ggal_gut_{1,2}.fq"
params.transcriptome = "${baseDir}/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = "results"
params.multiqc = "${baseDir}/multiqc"


// Component imports
include { RNASEQ } from './modules/rnaseq'
include { MULTIQC } from './modules/multiqc'


// Workflow block
workflow {

    log.info """\
      R N A S E Q - N F   P I P E L I N E
      ===================================
      transcriptome: ${params.transcriptome}
      reads        : ${params.reads}
      outdir       : ${params.outdir}
    """.stripIndent()

    read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true, flat: true)

    (fastqc_ch, quant_ch) = RNASEQ(read_pairs_ch, params.transcriptome)

    multiqc_files_ch = fastqc_ch.mix(quant_ch).collect()

    MULTIQC(multiqc_files_ch, params.multiqc)

    workflow.onComplete = {
        log.info(
            workflow.success
                ? "\nDone! Open the following report in your browser --> ${params.outdir}/multiqc_report.html\n"
                : "Oops .. something went wrong"
        )
    }
}
```

</div>

### Synopsis

The pipeline uses two imported components:

  <p style="padding-left: 40px;">&#8226; <code>RNASEQ</code>: a subworkflow that contains three processes:</p>

  <p style="padding-left: 80px;">&#8226; <code>INDEX</code>: creates a Salmon index from the transcriptome (runs once)</p>

  <p style="padding-left: 80px;">&#8226; <code>FASTQC</code>: analyzes each sample in parallel</p>

  <p style="padding-left: 80px;">&#8226; <code>QUANT</code>: quantifies transcripts for each sample after indexing completes</p>

  <p style="padding-left: 40px;">&#8226; <code>MULTIQC</code>: aggregates all quality control and quantification outputs into a comprehensive HTML report</p>

The `workflow` block uses `channel.fromFilePairs` to create a channel of paired-end read files. It passes the reads and transcriptome to the `RNASEQ` subworkflow, then mixes the FastQC and quantification outputs and passes them to `MULTIQC`.

<br>

### Get started

To run this pipeline:

 <p style="padding-left: 40px;">1. <a href="https://docs.seqera.io/nextflow/install">Install Nextflow</a> (version 25.10 or later)</p>

 <p style="padding-left: 40px;">2. <a href="https://docs.docker.com/get-started/get-docker/"> Install Docker</a></p>

 <p style="padding-left: 40px;">3. Run the pipeline directly from its GitHub repository:</p>

<div style="padding-left: 60px; margin-top: 1rem;">

```
nextflow run nextflow-io/rnaseq-nf -profile docker
```

</div>

<br>

See the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) GitHub repository for all of the pipeline code and the [rnaseq-nf tutorial](https://nextflow.io/docs/stable/tutorials/rnaseq-nf.html) for a full pipeline description.
