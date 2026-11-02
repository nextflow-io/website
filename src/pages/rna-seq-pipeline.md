---
title: RNA-Seq pipeline
layout: "@layouts/ExampleLayout.astro"
---

<div class="blg-summary example">
<h2>RNA-Seq pipeline</h2>

<p class="">
    This pipeline shows an example of a basic RNA-Seq analysis that performs quality control, transcript quantification, and result aggregation. The pipeline processes paired-end FASTQ files, generates quality control reports with FastQC, quantifies transcripts with Salmon, and produces a unified report with MultiQC.
</p>

```nextflow
#!/usr/bin/env nextflow

/*
 * Proof of concept of a RNAseq pipeline implemented with Nextflow
 */

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.reads` specify on the run command line `--reads some_value`.
 */

params.reads = null
params.transcriptome = null
params.multiqc = "${projectDir}/multiqc"


// import modules
include { RNASEQ } from './modules/rnaseq'
include { MULTIQC } from './modules/multiqc'

/*
 * main script flow
 */
workflow {
    main:
    log.info """\
      R N A S E Q - N F   P I P E L I N E
      ===================================
      transcriptome: ${params.transcriptome}
      reads        : ${params.reads}
      outdir       : ${workflow.outputDir}
    """.stripIndent()

    // Create a channel of paired-end reads from the glob in `params.reads`.
    // `checkIfExists: true` fails fast if no files match; `flat: true` flattens
    // the pair into the tuple so each emission is `[sample_id, read1, read2]`,
    // e.g. ['ggal_gut', ggal_gut_1.fq, ggal_gut_2.fq].
    read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true, flat: true)

    // Run the RNASEQ sub-workflow over each read pair against the transcriptome.
    // It emits a channel of tuples with the following elements:
    //   - per-sample ID
    //   - per-sample FastQC report directories
    //   - per-sample Salmon quantification directories
    samples_ch = RNASEQ(read_pairs_ch, params.transcriptome)

    // Extract the FastQC and Salmon outputs into a single stream (`flatMap`) and
    // collect them into a list (`collect`) so that MultiQC sees all results at
    // once, e.g. [fastqc/ggal_gut, quant/ggal_gut, fastqc/ggal_liver, ...].
    multiqc_files_ch = samples_ch
        .flatMap { id, fastqc, quant -> [fastqc, quant] }
        .collect()

    // Run MultiQC on the collected result files using the config in `params.multiqc`
    // to produce a single consolidated HTML report (`multiqc_report.html`).
    multiqc_report = MULTIQC(multiqc_files_ch, params.multiqc)

    // Publish per-sample results and aggregated MultiQC report as workflow outputs
    // The directory paths for these outputs are defined in the `output` block below.
    publish:
    samples = samples_ch.map { id, fastqc, quant -> [id: id, fastqc: fastqc, quant: quant] }
    multiqc_report = multiqc_report
}

/*
 * workflow outputs
 *
 * The `output` block defines the directory structure of published outputs.
 * The base output directory can be set using the `-output-dir` CLI option
 * or `outputDir` config option. It defaults to `results`.
 * By default, all published files are saved directly in the output directory.
 * The `path` directive can be used to route files to specific paths within
 * the output directory.
 */
output {
    // Publish the per-sample FastQC and Salmon results into
    // the `fastqc` and `quant` subdirectories, respectively.
    samples {
        // Per-sample results are separated further into subdirectories
        // by ID (e.g. `fastqc/ggal_gut`, `fastqc/ggal_liver`, etc).
        path { sample ->
            sample.fastqc >> "fastqc/${sample.id}"
            sample.quant >> "quant/${sample.id}"
        }
        // Create an index file (samplesheet) of the published
        // samples in the output directory, e.g. `results/samples.csv`.
        index {
            path 'samples.csv'
            header true
        }
    }

    // Publish the aggregated MultiQC HTML report in the
    // output directory, e.g. `results/multiqc_report.html`
    multiqc_report {
    }
}
```

</div>

### Synopsis

The pipeline uses two imported components:

- `RNASEQ`: a subworkflow that contains three processes:
  - `INDEX`: creates a Salmon index from the transcriptome (runs once)
  - `FASTQC`: generates a quality control report for each sample in parallel
  - `QUANT`: quantifies transcripts for each sample after indexing completes
- `MULTIQC`: aggregates all quality control and quantification outputs into a single HTML report

The `workflow` block uses `channel.fromFilePairs` to create a channel of paired-end read files. It passes the reads and transcriptome to the `RNASEQ` subworkflow, which emits the FastQC and Salmon results for each sample. These results are collected and passed to `MULTIQC`.

The `publish:` section and the `output` block publish the per-sample results to the `fastqc` and `quant` subdirectories, along with a `samples.csv` index file and the MultiQC report, in the output directory (`results` by default).

### Get started

To run this pipeline:

1. [Install Nextflow](https://docs.seqera.io/nextflow/install) (version 25.10 or later).
2. [Install Docker](https://docs.docker.com/get-started/get-docker/).
3. Run the pipeline directly from its GitHub repository:

   ```bash
   nextflow run nextflow-io/rnaseq-nf -profile docker
   ```

When the run completes, open `results/multiqc_report.html` in your browser to view the report.

See the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) GitHub repository for all of the pipeline code and the [rnaseq-nf tutorial](https://docs.seqera.io/nextflow/tutorials/rnaseq-nf) for a full pipeline description.
