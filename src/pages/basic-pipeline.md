---
title: Your first script
layout: "@layouts/ExampleLayout.astro"
---

<div class="blg-summary example">
<h2>Your first script</h2>

<p class="" >
    Your first script shows how to write a pipeline with two simple Bash processes. The first process splits a string into chunks, and the second process converts the lowercase letters in each chunk to uppercase. It also demonstrates how to publish pipeline outputs to named directories.
</p>

```groovy
// Default parameter input
params.str = "Hello world!"

// split process
process split {
    input:
    val x

    output:
    path 'chunk_*'

    script:
    """
    printf '${x}' | split -b 6 - chunk_
    """
}

// convert_to_upper process
process convert_to_upper {
    tag "$y"

    input:
    path y

    output:
    path 'upper_*'

    script:
    """
    cat $y | tr '[a-z]' '[A-Z]' > upper_${y}
    """
}

// Workflow block
workflow {
    main:
    ch_str = channel.of(params.str)                  // Create a channel using parameter input
    ch_chunks = split(ch_str)                        // Split string into chunks and create a named channel
    ch_upper = convert_to_upper(ch_chunks.flatten()) // Convert lowercase letters to uppercase letters

    publish:
    lower = ch_chunks.flatten()
    upper = ch_upper
}

output {
    lower {
        path 'lower'
    }
    upper {
        path 'upper'
    }
}
```

</div>

### Script synopsis

This script defines two processes:

- `split`: takes a string as input, splits it into 6-byte chunks, and writes the chunks to files with the prefix `chunk_`

- `convert_to_upper`: takes files as input, transforms their contents to uppercase letters, and writes the uppercase strings to files with the prefix `upper_`

The `split` output is emitted as a single element containing all chunk files. The `flatten` operator splits this combined element so that each file is treated as a sole element and processed independently by `convert_to_upper`.

The `workflow` block is organized into two sections:

- `main:`: defines the workflow logic and how processes are connected via channels

- `publish:`: declares which channels should be published as workflow outputs

The `output` block (outside the workflow) defines where and how each output should be published. In this example, the outputs from both processes are published in subdirectories (`lower` and `upper`) in the default results output directory (`params.outdir`).

### Try it

To run this pipeline:

1. Install Nextflow
   - See [Installation](https://docs.seqera.io/nextflow/install) for more information
2. Create a new file named `main.nf` in your current directory
3. Copy and save the above script to your new file
4. Run your pipeline using the following command:

   ```
   nextflow run main.nf
   ```

See [Your first script](https://docs.seqera.io/nextflow/your-first-script) for more information about this pipeline.
