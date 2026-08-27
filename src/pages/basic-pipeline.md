---
title: Basic pipeline
layout: "@layouts/ExampleLayout.astro"
---

<div class="blg-summary example">
<h2>Basic pipeline</h2>

<p class="" >
    This pipeline shows how to write a pipeline with two simple Bash processes. The first process splits a string into chunks, and the second process converts the lowercase letters in each chunk to uppercase. It also shows how to publish pipeline outputs to named directories.
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
    ch_str = channel.of(params.str)
    ch_chunks = split(ch_str)
    ch_upper = convert_to_upper(ch_chunks.flatten())

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

### Synopsis

This pipeline defines two processes:

  <p style="padding-left: 40px;">&#8226; <code>split</code>: takes a string as input, splits it into 6-byte chunks, and writes the chunks to files with the prefix <code>chunk_</code>

  <p style="padding-left: 40px;">&#8226; <code>convert_to_upper</code>: takes files as input, transforms their contents to uppercase letters, and writes the uppercase strings to files with the prefix <code>upper_</code>

The `split` output is emitted as a single element containing all chunk files. The `flatten` operator splits this combined element so that each file is treated as a sole element and processed independently by `convert_to_upper`.

The `workflow` block is organized into two sections:

  <p style="padding-left: 40px;">&#8226; <code>main:</code> defines the workflow logic and how processes are connected via channels</p>

  <p style="padding-left: 40px;">&#8226; <code>publish:</code> declares which channels should be published as workflow outputs</p>

The `output` block (outside the workflow) defines where and how each output should be published. In this example, the outputs from both processes are published in subdirectories (`lower` and `upper`) in the default results output directory (`params.outdir`).

<br>

### Get started

To run this pipeline:

 <p style="padding-left: 40px;">1. <a href="https://docs.seqera.io/nextflow/install">Install Nextflow</a> (version 25.10 or later)</p>

 <p style="padding-left: 40px;">2. Create a new file named <code>main.nf</code> in your current directory</p>

 <p style="padding-left: 40px;">3. Copy and save the above script to your new file</p>

 <p style="padding-left: 40px;">4. Run your pipeline:</p>

<div style="padding-left: 60px; margin-top: 1rem;">

```
nextflow run main.nf
```

</div>

<br>

See [Your first script](https://docs.seqera.io/nextflow/your-first-script) for more information about this pipeline.
