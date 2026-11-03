---
title: BLAST pipeline
layout: "@layouts/ExampleLayout.astro"
---

<div class="blg-summary example">
<h2>BLAST pipeline</h2>

<p class="">
    This pipeline splits a FASTA file into chunks and executes a BLAST query for each chunk in parallel. Then, all the sequences for the top hits are collected and merged into a single result file.
</p>

```nextflow
nextflow.enable.types = true

/*
 * Pipeline parameters. Each parameter can be overridden on the command line (e.g. `--query query.fa`)
 */
params {
    // Query sequences (FASTA)
    query: Path

    // BLAST database (directory and database name)
    db: String

    // Name of the output file with the matching sequences
    out: String = "result.txt"

    // Number of sequences per query chunk
    chunkSize: Integer = 100
}


// Entry workflow
workflow {
    main:
    db_dir = file(params.db).parent
    db_name = file(params.db).name

    // Split the query file into chunks of `chunkSize` sequences
    ch_fasta = channel.of(params.query)
        .flatMap { fasta -> fasta.splitFasta(by: params.chunkSize, file: true) }

    // Run BLAST on each chunk and emit the top hits
    ch_hits = blast(ch_fasta, db_dir, db_name)

    // Extract the sequences of the top hits from the database
    ch_sequences = extract(ch_hits, db_dir, db_name)

    // Merge the matching sequences into a single file and print it
    val_result = merge_sequences(ch_sequences.collect(), params.out)
    val_result.view { result -> "matching sequences:\n${result.text}" }

    publish:
    result = val_result
}


// Pipeline outputs (published to the `results` directory by default)
output {
    result {
    }
}


// Search each query chunk against the database and keep the IDs of the top 10 hits
process blast {
    input:
    query: Path
    db_dir: Path
    db_name: String

    output:
    file('top_hits')

    script:
    """
    blastp -db ${db_dir}/${db_name} -query ${query} -outfmt 6 > blast_result
    cat blast_result | head -n 10 | cut -f 2 > top_hits
    """
}


// Fetch the sequences of the top hits from the database
process extract {
    input:
    top_hits: Path
    db_dir: Path
    db_name: String

    output:
    file('sequences')

    script:
    """
    blastdbcmd -db ${db_dir}/${db_name} -entry_batch ${top_hits} | head -n 10 > sequences
    """
}


// Concatenate the matching sequences into a single file
process merge_sequences {
    input:
    sequences: Bag<Path>
    name: String

    stage:
    stageAs sequences, 'sequences_*'

    output:
    file(name)

    script:
    """
    cat sequences_* > ${name}
    """
}
```

</div>

### Synopsis

This pipeline defines three processes:

- `blast`: searches a chunk of query sequences against the BLAST database and writes the IDs of the top 10 hits to a file.
- `extract`: fetches the sequences of the top hits from the BLAST database.
- `merge_sequences`: concatenates the matching sequences into a single file.

The `workflow` block uses the `splitFasta` function with the `flatMap` operator to split the query file into chunks of `chunkSize` sequences, so that each chunk is processed by `blast` and `extract` in parallel. The `collect` operator gathers the sequences from every chunk so that `merge_sequences` runs once, and the `view` operator prints the result to the terminal.

The default values for `query` and `db` are set in the pipeline's `nextflow.config`.

The `output` block publishes the merged sequences to the output directory (`results` by default).

### Get started

To run this pipeline:

1. [Install Nextflow](https://docs.seqera.io/nextflow/install) (version 26.10 or later).
2. [Install Docker](https://docs.docker.com/get-started/get-docker/).
3. Run the pipeline directly from its GitHub repository:

   ```bash
   nextflow run nextflow-io/blast-example -profile docker
   ```

The matching sequences are printed to the console and saved to `results/result.txt`.

See the [blast-example](https://github.com/nextflow-io/blast-example) GitHub repository for all of the pipeline code.
