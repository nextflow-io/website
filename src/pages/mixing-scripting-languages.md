---
title: Mixed script pipeline
layout: "@layouts/ExampleLayout.astro"
---

<div class="blg-summary example">
<h2>Mixed script pipeline</h2>

<p class="" >
    This pipeline shows how to use different scripting languages in a single Nextflow pipeline. The first process uses a Perl script to generate random number pairs, and the second process uses a Python script to calculate the average of each pair.
</p>

```groovy
// Default parameter input
params.range = 100

// Perl process
process perlTask {
    output:
    stdout

    script:
    """
    #!/usr/bin/env perl
    use strict;
    use warnings;

    my \$count;
    my \$range = ${params.range};
    for (\$count = 0; \$count < 10; \$count++) {
        print rand(\$range) . ', ' . rand(\$range) . "\\n";
    }
    """
}

// Python process
process pyTask {
    input:
    stdin

    output:
    stdout

    script:
    """
    #!/usr/bin/env python
    import sys

    x = 0
    y = 0
    lines = 0
    for line in sys.stdin:
        items = line.strip().split(",")
        x += float(items[0])
        y += float(items[1])
        lines += 1

    print("avg: %s - %s" % ( x/lines, y/lines ))
    """
}

// Workflow block
workflow {
    perlTask | pyTask | view
}
```

</div>

### Synopsis

This pipeline defines two processes:

  <p style="padding-left: 40px;">&#8226; <code>perlTask</code>: executes a Perl script that generates 10 pairs of random numbers within a configurable range and prints them to standard output</p>

  <p style="padding-left: 40px;">&#8226; <code>pyTask</code>: executes a Python script that reads the pairs from standard input and calculates the average of each coordinate</p>

Each process uses a shebang declaration at the start of its script block to specify the scripting language. Nextflow detects the shebang and executes the script using the appropriate interpreter.

The `workflow` block pipes the output of `perlTask` directly into `pyTask`, and then passes the result to the `view` operator to print it to the terminal.

<br>

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

See [Scripts à la carte](https://docs.seqera.io/nextflow/process#scripts-%C3%A0-la-carte) for more information about using multiple scripting languages in Nextflow.
