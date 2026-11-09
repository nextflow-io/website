---
title: Basic pipeline
layout: "@layouts/ExampleLayout.astro"
---

<div class="blg-summary example">
<h2>Basic pipeline</h2>

<p class="">
    This pipeline shows how to write a pipeline with two simple Bash processes. The first process splits a string into chunks, and the second process converts the lowercase letters in each chunk to uppercase. It also shows how to publish pipeline outputs to named directories.
</p>

```nextflow title="main.nf"
nextflow.enable.types = true

// Pipeline parameters
params {
    str: String = "Hello world!"
}

// Entry workflow
workflow {
    main:
    // Create a channel of the pipeline input
    ch_str = channel.of(params.str)

    // Split string into chunks and emit each chunk separately
    ch_chunks = split(ch_str).flatMap { chunks -> chunks }

    // Convert lowercase letters to uppercase letters
    ch_upper = convert_to_upper(ch_chunks)

    publish:
    lower = ch_chunks
    upper = ch_upper
}

// Pipeline outputs
output {
    lower {
        path 'lower'
    }
    upper {
        path 'upper'
    }
}

// split process
process split {
    input:
    x: String

    output:
    files('chunk_*')

    script:
    """
    printf '${x}' | split -b 6 - chunk_
    """
}

// convert_to_upper process
process convert_to_upper {
    tag "$y"

    input:
    y: Path

    output:
    file("upper_${y}")

    script:
    """
    cat $y | tr '[a-z]' '[A-Z]' > upper_${y}
    """
}
```

</div>

### Synopsis

The script enables [static typing](https://docs.seqera.io/nextflow/static-typing) with the `nextflow.enable.types` feature flag.

This script defines two processes:

- `split`: takes a string input, splits it into 6-character chunks, and writes the chunks to files with the prefix `chunk_`.
- `convert_to_upper`: takes an input file, transforms its contents to uppercase letters, and writes the uppercase text to a file with the prefix `upper_`.

The script also defines a `params` block, a `workflow` block, and an `output` block:

- The `params` block declares the pipeline parameters.
- The `workflow` block defines how processes are connected via channels (`main:`) and which channels are published as outputs (`publish:`). In this example, the `split` process emits a collection of files, and the `flatMap` operator emits each file separately to the `convert_to_upper` process.
- The `output` block defines where each output is published. In this example, the outputs from both processes are published to subdirectories (`lower` and `upper`) in the output directory (`results` by default).

### Get started

To run this pipeline:

1. [Install Nextflow](https://docs.seqera.io/nextflow/install) (version 26.10 or later).
2. Create a new file named `main.nf` in your current directory.
3. Copy and save the above script to your new file.
4. Run your pipeline:

   ```bash
   nextflow run main.nf
   ```

See [Your first script](https://docs.seqera.io/nextflow/your-first-script) for more information about this pipeline.
