---
title: Mixing scripting languages
layout: "@layouts/ExampleLayout.astro"
---

<div class="blg-summary example">
<h2>Mixing scripting languages</h2>

<p class="">
    This pipeline shows how to use different scripting languages in a single Nextflow pipeline. The first process uses an R script to generate random number pairs, and the second process uses a Python script to calculate the average of each coordinate.
</p>

```nextflow
nextflow.enable.types = true

// Pipeline parameters
params {
    range: Integer = 100
}

// Entry workflow
workflow {
    pairs = rTask(params.range)
    result = pyTask(pairs)
    result.view()
}

// R process
process rTask {
    input:
    range: Integer

    output:
    stdout()

    script:
    """
    #!/usr/bin/env Rscript
    x <- runif(10, max = ${range})
    y <- runif(10, max = ${range})
    writeLines(paste(x, y, sep = ", "))
    """
}

// Python process
process pyTask {
    input:
    pairs: String

    stage:
    stdin pairs

    output:
    stdout()

    script:
    """
    #!/usr/bin/env python3
    import sys

    x = 0
    y = 0
    lines = 0
    for line in sys.stdin:
        items = line.strip().split(",")
        x += float(items[0])
        y += float(items[1])
        lines += 1

    print("avg: %s - %s" % (x / lines, y / lines))
    """
}
```

</div>

### Synopsis

This pipeline defines two processes:

- `rTask`: executes an R script that generates 10 pairs of random numbers within a configurable range and prints them to standard output.
- `pyTask`: executes a Python script that reads the pairs from standard input and calculates the average of each coordinate.

Each process uses a shebang declaration at the start of its script block to specify the scripting language. Nextflow detects the shebang and executes the script using the appropriate interpreter.

The `workflow` block passes the `range` parameter to `rTask`, passes the output of `rTask` to `pyTask`, and then uses the `view` operator to print the result to the terminal.

### Get started

To run this pipeline:

1. [Install Nextflow](https://docs.seqera.io/nextflow/install) (version 26.10 or later).
2. Make sure that R (`Rscript`) and Python 3 (`python3`) are installed.
3. Create a new file named `main.nf` in your current directory.
4. Copy and save the above script to your new file.
5. Run your pipeline:

   ```bash
   nextflow run main.nf
   ```

See [Scripts à la carte](https://docs.seqera.io/nextflow/process#scripts-%C3%A0-la-carte) for more information about using multiple scripting languages in Nextflow.
