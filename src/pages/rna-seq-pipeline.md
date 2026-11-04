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

nextflow.enable.types = true

include { RNASEQ } from './modules/rnaseq'
include { MULTIQC } from './modules/multiqc'

include { Sample ; AlignedSample } from './modules/rnaseq'

params {
    // The input samplesheet of paired-end reads. Nextflow loads the
    // samplesheet into a channel of records, one record per row.
    reads: Channel<Sample>

    // The input transcriptome file
    transcriptome: Path

    // Directory containing multiqc configuration
    multiqc: Path = "${projectDir}/multiqc"
}

workflow {
    main:

    // Run the RNASEQ sub-workflow over each read pair against the transcriptome.
    // It emits a channel of records with the following fields:
    //   - id: per-sample ID
    //   - fastqc: per-sample FastQC directories
    //   - quant: per-sample Salmon quantification directories
    samples_ch = RNASEQ(params.reads, params.transcriptome)

    // Extract the FastQC and Salmon outputs into a single stream (`flatMap`) and
    // collect them into a list (`collect`) so that MultiQC sees all results at
    // once, e.g. [fastqc/ggal_gut, quant/ggal_gut, fastqc/ggal_liver, ...].
    multiqc_files_ch = samples_ch
        .flatMap { r -> [r.fastqc, r.quant] }
        .collect()

    // Run MultiQC on the collected result files using the config in `params.multiqc`
    // to produce a single consolidated HTML report (`multiqc_report.html`).
    multiqc_report = MULTIQC(multiqc_files_ch, params.multiqc)

    // Publish per-sample results and aggregated MultiQC report as workflow outputs
    // The directory paths for these outputs are defined in the `output` block below.
    publish:
    samples = samples_ch
    multiqc_report = multiqc_report
}

output {
    // Publish the per-sample FastQC and Salmon results into
    // the `fastqc` and `quant` subdirectories, respectively.
    samples: Channel<AlignedSample> {
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
    multiqc_report: Path {
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

The `params` block declares the pipeline inputs with their types. The `reads` parameter is a `Channel<Sample>`, so Nextflow loads the input samplesheet into a channel of `Sample` records, one per row. The `workflow` block passes the reads and transcriptome to the `RNASEQ` subworkflow, which emits a record with the FastQC and Salmon results for each sample. These results are collected and passed to `MULTIQC`.

The `publish:` section and the `output` block publish the per-sample results to the `fastqc` and `quant` subdirectories, along with a `samples.csv` index file and the MultiQC report, in the output directory (`results` by default).

### Get started

To run this pipeline:

1. [Install Nextflow](https://docs.seqera.io/nextflow/install) (version 26.10 or later).
2. [Install Docker](https://docs.docker.com/get-started/get-docker/).
3. Run the pipeline directly from its GitHub repository:

   ```bash
   nextflow run nextflow-io/rnaseq-nf -profile docker
   ```

When the run completes, open `results/multiqc_report.html` in your browser to view the report.

See the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) GitHub repository for all of the pipeline code and the [rnaseq-nf tutorial](https://docs.seqera.io/nextflow/tutorials/rnaseq-nf) for a full pipeline description.
